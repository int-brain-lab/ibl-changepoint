function [params,data,refitted_flag] = fit_model(model_name,data,Nopts,vbmc_flag,refit_flags,opt_init,save_flag)
%FIT_MODEL Fit model MODEL_NAME to dataset DATA.

% # restarts for fitting procedures (MLE and variational inference)
if nargin < 3 || isempty(Nopts); Nopts = [10,5]; end
if numel(Nopts) > 1; Nvbmc = Nopts(2); else; Nvbmc = ceil(Nopts(1)/2); end

% Get approximate posteriors with Variational Bayesian Monte Carlo
if nargin < 4 || isempty(vbmc_flag); vbmc_flag = false; end

% Force refits even if fit already exists
if nargin < 5 || isempty(refit_flags); refit_flags = false; end

% Starting point(s) for the first fit
if nargin < 6; opt_init = []; end
if isstruct(opt_init); opt_init = {opt_init}; end
if iscell(opt_init)
    for iOpt = 1:numel(opt_init)
        temp(iOpt,:) = opt_init{iOpt}.theta;
    end
    opt_init = temp;
end

% Save fits
if nargin < 7 || isempty(save_flag); save_flag = false; end

refitted_flag = false;  % Check if model has been refitted

% Load data file if passed as string
if ischar(data); data = read_data_from_csv(data); end

% Load fit if model was already fitted to this dataset
params = load_model_fit(data.fullname,model_name);

%% Maximum likelihood fit
if isempty(params) || refit_flags(1)
    % Initialize model/params struct
    params = params_new(model_name,data);
    params.mle_fits = [];    % Reset MLE fits struct
end

bounds = setup_params([],params);       % Get parameter bounds

% Find starting point from other models
x0_base = [];
if ~isempty(params.model_startingfit)
    params1 = load_model_fit(data.fullname,params.model_startingfit);
    if ~isempty(params1)
        x0_base = NaN(1,numel(params.names));
        %----------------------------------------------------------
        % This should be assigned based on matching parameter names!
        %----------------------------------------------------------
        x0_base(1:numel(params1.theta)) = params1.theta;  % TO BE FIXED
        fprintf('Reading starting point from model %s.\n', params.model_startingfit);
    end
end

% Fits all models but psychometric curves with BADS    
bads_flag = ~contains(params.model_name,{'psychofun'});

% Separate optimization runs
for iOpt = 1:Nopts(1)
    
    % Skip optimization run if it already exists
    if ~isempty(params.mle_fits) && numel(params.mle_fits.nll) >= iOpt
        continue;
    end
    
    if size(opt_init,1) >= iOpt     % Use provided starting points
        x0 = opt_init(iOpt,:);
    else
        if iOpt == 1
            x0 = bounds.x0;
        else
            x0 = rand(1,numel(bounds.PLB)).*(bounds.PUB-bounds.PLB) + bounds.PLB;
        end
        if ~isempty(x0_base) && iOpt <= ceil(Nopts(1)/3)
            x0(~isnan(x0_base)) = x0_base(~isnan(x0_base));
        end
    end

    if bads_flag
        MaxFunEvals = 2e3;
        badopts = bads('defaults');
        badopts.MaxFunEvals = MaxFunEvals;        
        [x,fval] = bads(@(x_)nllfun(x_,params,data),...
            x0,bounds.LB,bounds.UB,bounds.PLB,bounds.PUB,[],badopts);
    else
        fminopts.Display = 'iter';
        [x,fval] = fmincon(@(x_)nllfun(x_,params,data),...
            x0,[],[],[],[],bounds.LB,bounds.UB,[],fminopts);            
    end
    
    params.mle_fits.x0(iOpt,:) = x0;
    params.mle_fits.x(iOpt,:) = x;
    params.mle_fits.nll(iOpt) = fval;
    
    refitted_flag = true;
    params.mle_fits.x_best = [];
    params.mle_fits.nll_best = [];
    
    if save_flag; save_model_fit(data.fullname,params); end
end

params.mle_fits.x
params.mle_fits.nll
[nll_best,idx_best] = min(params.mle_fits.nll);
x_best = params.mle_fits.x(idx_best,:);
params = setup_params(x_best,params);

params.mle_fits.x_best = x_best;
params.mle_fits.nll_best = nll_best;

if refitted_flag && save_flag; save_model_fit(data.fullname,params); end

if vbmc_flag
    
    if ~isfield(params,'vbmc_fits') || isempty(params.vbmc_fits) || refit_flags(min(2,end))
        params.vbmc_fits = [];
    end
        
    bounds = setup_params([],params);       % Get parameter bounds    
    x0 = params.theta;

    vbmc_opts = vbmc('defaults');
    % vbmc_opts.Plot = 'on';
    vbmc_opts.NSgpMaxMain = 0;
    vbmc_opts.RetryMaxFunEvals = vbmc_opts.MaxFunEvals;

%         w.SearchCacheFrac = 0.1; w.HPDSearchFrac = 0.9; w.HeavyTailSearchFrac = 0; w.MVNSearchFrac = 0; w.SearchAcqFcn = @vbmc_acqpropregt; w.StopWarmupThresh = 0.1; w.SearchCMAESVPInit = false;
%         vbmc_opts.WarmupOptions = w; vbmc_opts.TolStableWarmup = 5; vbmc_opts.FastWarmup = true; vbmc_opts.NSgpMaxWarmup = 8;

    % Assume smoothed trapezoidal prior over finite box
    logprior = @(x) log(msplinetrapezpdf(x,bounds.LB,bounds.PLB,bounds.PUB,bounds.UB));
    
    for iOpt = 1:Nvbmc        
        % Skip variational optimization run if it already exists
        if ~isempty(params.vbmc_fits) && numel(params.vbmc_fits.vps) >= iOpt ...
                && ~isempty(params.vbmc_fits.vps{iOpt})
            continue;
        end
        
        [vp,~,~,~,output] = ...
            vbmc(@(x_) -nllfun(x_,params,data)+logprior(x_), ...
            x0,bounds.LB,bounds.UB,bounds.PLB,bounds.PUB,vbmc_opts);

        % Store results
        params.vbmc_fits.vps{iOpt} = vp;
        params.vbmc_fits.outputs{iOpt} = output;
        
        params.vbmc_fits.diagnostics = [];
        refitted_flag = true;
        if save_flag; save_model_fit(data.fullname,params); end
    end

    if refitted_flag
        [exitflag,best,idx_best,stats] = vbmc_diagnostics(params.vbmc_fits.vps);
        params.vbmc_fits.diagnostics.exitflag = exitflag;
        params.vbmc_fits.diagnostics.best = best;
        params.vbmc_fits.diagnostics.idx_best = idx_best;
        params.vbmc_fits.diagnostics.stats = stats;
        
        if save_flag; save_model_fit(data.fullname,params); end
    end

    
end

end
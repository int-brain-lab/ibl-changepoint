function [params,data,refitted_flag] = fit_model(model_name,data,Nopts,vbmc_flag,refit_flags,opt_init)
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

refitted_flag = false;  % Check if model has been refitted

% Load data file if passed as string
if ischar(data); data = read_data_from_csv(data); end

% Load fit if model was already fitted to this dataset
params = load_model_fit(data.name,model_name);

%% Maximum likelihood fit
if isempty(params) || refit_flags(1)
    params = params_new(model_name,data);   % Initialize model/params struct
    bounds = setup_params([],params);       % Get parameter bounds

    MaxFunEvals = 2e3;
    x = []; nll = [];
    badopts = bads('defaults');
    badopts.MaxFunEvals = MaxFunEvals;

    % Find starting point from other models
    x0_base = [];
    if ~isempty(params.model_startingfit)
        params1 = load_model_fit(data.name,params.model_startingfit);
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

    for iOpt = 1:Nopts(1)
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
            [x(iOpt,:),nll(iOpt)] = bads(@(x_)nllfun(x_,params,data),...
                x0,bounds.LB,bounds.UB,bounds.PLB,bounds.PUB,[],badopts);
        else
            fminopts.Display = 'iter';
            [x(iOpt,:),nll(iOpt)] = fmincon(@(x_)nllfun(x_,params,data),...
                x0,[],[],[],[],bounds.LB,bounds.UB,[],fminopts);            
        end
    end

    x
    nll
    [nll_best,idx_best] = min(nll);
    x_best = x(idx_best,:);
    params = setup_params(x_best,params);
    
    refitted_flag = true;
end

if vbmc_flag
    if ~isfield(params,'vbmc_fit') || isempty(params.vbmc_fit) || refit_flags(min(2,end))   
        if isfield(params,'vbmc_fit'); params = rmfield(params,'vbmc_fit'); end
        
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
        
        vbmc_fit = [];
        for iOpt = 1:Nvbmc
            [vp,~,~,~,output] = ...
                vbmc(@(x_) -nllfun(x_,params,data)+logprior(x_), ...
                x0,bounds.LB,bounds.UB,bounds.PLB,bounds.PUB,vbmc_opts);

            % Store results
            vbmc_fit.vps{iOpt} = vp;
            vbmc_fit.outputs{iOpt} = output;
        end
        
        [exitflag,best,idx_best,stats] = vbmc_diagnostics(vbmc_fit.vps);
        vbmc_fit.diagnostics.exitflag = exitflag;
        vbmc_fit.diagnostics.best = best;
        vbmc_fit.diagnostics.idx_best = idx_best;
        vbmc_fit.diagnostics.stats = stats;
        
        params.vbmc_fit = vbmc_fit;
        refitted_flag = true;
    end
end

end
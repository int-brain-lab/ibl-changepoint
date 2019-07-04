function [params,data,refitted_flag] = fit_model(model_name,data,Nopts,vbmc_flag,refit_flags,opt_init,save_flag,empirical_list)
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

% List of datasets for empirical Bayes prior
if nargin < 8
    empirical_list = [];
end
empirical_bayes = iscell(empirical_list) || (~isempty(empirical_list) && ~(empirical_list == false));

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
        
    nvars = numel(bounds.PLB);  % Number of dimensions
    
    if size(opt_init,1) >= iOpt     % Use provided starting points
        x0 = opt_init(iOpt,:);
    else
        if iOpt == 1
            x0 = bounds.x0;
        else
            x0 = rand(1,nvars).*(bounds.PUB-bounds.PLB) + bounds.PLB;
        end
        if ~isempty(x0_base) && iOpt <= ceil(Nopts(1)/3)
            x0(~isnan(x0_base)) = x0_base(~isnan(x0_base));
        end
    end
    
    % Evaluate nLL on a bunch of quasirandom points, start from best
    Ninit = 5*nvars;    
    P = sobolset(nvars);
    P = scramble(P,'MatousekAffineOwen');
    xx = bsxfun(@plus,bsxfun(@times,net(P,Ninit),bounds.PUB-bounds.PLB),bounds.PLB);    
    x0_list = [x0; xx];

    nll0 = zeros(size(x0_list,1),1);
    for i0 = 1:size(x0_list,1)
        nll0(i0) = nllfun(x0_list(i0,:),params,data);        
    end
    [~,idx0] = min(nll0);
    x0 = x0_list(idx0,:);

    % Run optimization from best initial point
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
    params = compute_mle_stats(params);
    params.mle_fits.ndata = size(data.tab,1);
    
    if save_flag; save_model_fit(data.fullname,params); end
end

params.mle_fits.x
params.mle_fits.nll

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
    % vbmc_opts.NoiseSize = 0.1;
    % vbmc_opts.gpNoiseFun = [1 0 1];
    vbmc_opts.SGDStepSize = 5e-4;           % Very narrow posteriors
    vbmc_opts.Bandwidth = 0;             % Deal with high-frequency noise
    vbmc_opts.UncertaintyHandling = true;
    
%         w.SearchCacheFrac = 0.1; w.HPDSearchFrac = 0.9; w.HeavyTailSearchFrac = 0; w.MVNSearchFrac = 0; w.SearchAcqFcn = @vbmc_acqpropregt; w.StopWarmupThresh = 0.1; w.SearchCMAESVPInit = false;
%         vbmc_opts.WarmupOptions = w; vbmc_opts.TolStableWarmup = 5; vbmc_opts.FastWarmup = true; vbmc_opts.NSgpMaxWarmup = 8;

    
    % Compute empirical Bayes (truncated) Student's t prior over parameters
    if empirical_bayes        
        if ~iscell(empirical_list)
            empirical_list = get_mice_list([],[]);
        end
        theta = [];
        for iMouse = 1:numel(empirical_list)
            params1 = load_model_fit(empirical_list{iMouse},model_name);
            if ~isempty(params1); theta = [theta; params1.theta]; end
        end
        prior_mean = mean(theta);
        prior_std = std(theta);
        df = 3;     % Use Student's t with nu = 3, milder evidence
        logprior = @(x) log(mtrunctpdf(x,bounds.LB,prior_mean-prior_std,prior_mean+prior_std,bounds.UB,df));
        [prior_mean; prior_std]
        
        width = bounds.UB - bounds.LB;
        bounds.PLB = max(prior_mean - 2*prior_std,bounds.LB+0.025*width);
        bounds.PUB = min(prior_mean + 2*prior_std,bounds.UB-0.025*width);        
    else
        % Assume smoothed trapezoidal prior over finite box
        logprior = @(x) log(msplinetrapezpdf(x,bounds.LB,bounds.PLB,bounds.PUB,bounds.UB));    
        % logprior = @(x) log(mtruncgausspdf(x,bounds.LB,bounds.PLB,bounds.PUB,bounds.UB));
    end
    
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
        
        refitted_flag = true;
        params = compute_vbmc_stats(params);
        
        if save_flag; save_model_fit(data.fullname,params); end
    end

    
end

end

%--------------------------------------------------------------------------
function params = compute_mle_stats(params)

[nll_best,idx_best] = min(params.mle_fits.nll);
x_best = params.mle_fits.x(idx_best,:);
params = setup_params(x_best,params);

params.mle_fits.x_best = x_best;
params.mle_fits.nll_best = nll_best;

end

%--------------------------------------------------------------------------
function params = compute_vbmc_stats(params)

[exitflag,best,idx_best,stats] = vbmc_diagnostics(params.vbmc_fits.vps);
params.vbmc_fits.diagnostics.exitflag = exitflag;
params.vbmc_fits.diagnostics.best = best;
params.vbmc_fits.diagnostics.idx_best = idx_best;
params.vbmc_fits.diagnostics.stats = stats;

end
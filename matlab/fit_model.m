function [params,data,refitted_flag] = fit_model(model_name,data,Nopts,methods_flags,refit_flags,opt_init,save_flag,empirical_list)
%FIT_MODEL Fit model MODEL_NAME to dataset DATA.

% # restarts for fitting procedures (MLE and variational inference)
if nargin < 3 || isempty(Nopts); Nopts = [10,5]; end
if numel(Nopts) > 1; Nvbmc = Nopts(2); else; Nvbmc = ceil(Nopts(1)/2); end

% Which advanced fitting methods should be applied?
if nargin < 4 || isempty(methods_flags); methods_flags = false; end

Nmethods = 3;
if numel(methods_flags) < Nmethods; methods_flags = [methods_flags, false(1,Nmethods-numel(methods_flags))]; end

hmmfit_flag = methods_flags(1); % Compute state transitions with Hidden Markov model fitting
vbmc_flag = methods_flags(2);   % Get approximate posteriors with Variational Bayesian Monte Carlo
mcmc_flag = methods_flags(3);   % Get approximate posteriors via MCMC

% Force refits even if fit already exists (default is FALSE)
if nargin < 5 || isempty(refit_flags); refit_flags = false; end
if isscalar(refit_flags); refit_flags = refit_flags*ones(1,4); end
% Flags are: 1 for MLE, 2 for HMM, 3 for VBMC, 4 for MCMC

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
nvars = numel(bounds.PLB);              % Number of dimensions

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
            
    x0 = NaN(1,nvars);
    if size(opt_init,1) >= iOpt     % Use provided starting points
        x0 = opt_init(iOpt,:);
    else
        if iOpt == 1; x0 = bounds.x0; end
        if ~isempty(x0_base) && iOpt <= ceil(Nopts(1)/3)
            x0(~isnan(x0_base)) = x0_base(~isnan(x0_base));
        end
    end
    
    % Fit model via optimization
    [x,fval,x0] = optimize_model(params,data,x0);
        
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

%% Cross-validated maximum-likelihood fit

% Build only-even and only-odd sessions datasets
sessions = unique(data.tab(:,2));
tab = data.tab;
rmsess = sessions(1:2:end);    
for iSess = 1:numel(rmsess)
    tab(tab(:,2) == rmsess(iSess),:) = [];
end
data_fold{1} = format_data(tab,data.filename,data.fullname);

tab = data.tab;
rmsess = sessions(2:2:end);    
for iSess = 1:numel(rmsess)
    tab(tab(:,2) == rmsess(iSess),:) = [];
end
data_fold{2} = format_data(tab,data.filename,data.fullname);

if ~isfield(params.mle_fits,'cv') || refit_flags(1) || refitted_flag
    params.mle_fits.cv = [];
end

% Skip cross-validation run if there is only one session
if numel(sessions) > 1
    for iFold = 1:2    
        if iFold == 1; iTest = 2; else; iTest = 1; end

        % Skip cross-validation run if it already exists
        if ~isempty(params.mle_fits.cv) && numel(params.mle_fits.cv.nll_test) >= iFold
            continue;
        end

        params_train = params_new(model_name,data_fold{iFold});

        % Cross-validation via optimization
        x0 = [params.mle_fits.x0; params.mle_fits.x];
        [x,fval,x0] = optimize_model(params_train,data_fold{iFold},x0);    

        % Evaluate log likelihood on other fold
        params_test = params_new(model_name,data_fold{iTest});
        fun = @(x_) sum(nllfun(x_,params_test,data_fold{iTest}));

        params.mle_fits.cv.x0{iFold} = x0;
        params.mle_fits.cv.x{iFold} = x;
        params.mle_fits.cv.nll_train{iFold} = fval;

        fval_test = fun(x);
        params.mle_fits.cv.nll_test(iFold) = fval_test;
        refitted_flag = true;

        if save_flag; save_model_fit(data.fullname,params); end
    end
end

%% Hidden Markov model fit
if hmmfit_flag
    
    if ~isfield(params,'hmm_fits') || isempty(params.hmm_fits) || refit_flags(2)
        params.hmm_fits = [];
    end    
    
    K_vec = 1:5;  % Run HMM fit with these components
            
    x0 = params.theta;
    hmmopts = hmmfit('defaults');
    hmmopts.MaxIter = 10;
    hmmopts.TolFun = 0.01;
    
    for k = 1:numel(K_vec)
        K = K_vec(k);
        
        if isempty(params.hmm_fits) || numel(params.hmm_fits.hmm) < K ...
                || isempty(params.hmm_fits.hmm{K})
            
            if K == 1
                [hmm,ll,exitflag,output] = compute_hmm1(params.mle_fits);
                hmm.ell(1,:) = exp(-nllfun(hmm.mu,params,data));
            else
                [hmm,ll,exitflag,output] = hmmfit(@(x_) -nllfun(x_,params,data),[],K,x0, ...
                    bounds.LB,bounds.UB,bounds.PLB,bounds.PUB,hmmopts);
            end

            params.hmm_fits.hmm{K} = hmm;
            params.hmm_fits.nll(K) = -ll;
            params.hmm_fits.exitflag(K) = exitflag;
            params.hmm_fits.output{K} = output;

            params = compute_hmm_stats(params);
            
            if save_flag; save_model_fit(data.fullname,params); end            
        end   
        
        x0 = params.hmm_fits.hmm{K}.mu;    % Use current solution as next starting point
        
    end
    
end

%% Variational Bayesian Monte Carlo (VBMC)
if vbmc_flag
    
    if ~isfield(params,'vbmc_fits') || isempty(params.vbmc_fits) || refit_flags(3)
        params.vbmc_fits = [];
    end
        
    bounds = setup_params([],params);       % Get parameter bounds    
    x0 = params.theta;

    vbmc_opts = vbmc('defaults');
    % vbmc_opts.Plot = 'on';
    vbmc_opts.NSgpMaxMain = 0;
    vbmc_opts.ActiveSampleVPUpdate = true;
    vbmc_opts.ActiveSampleGPUpdate = true;    
%    vbmc_opts.SearchAcqFcn = @acqmireg_vbmc;
    vbmc_opts.WarpRotoScaling = 1;
    vbmc_opts.RetryMaxFunEvals = vbmc_opts.MaxFunEvals;
    % vbmc_opts.NoiseSize = 0.1;
    % vbmc_opts.gpNoiseFun = [1 0 1];
    vbmc_opts.Bandwidth = 0;             % Deal with high-frequency noise
    % vbmc_opts.UncertaintyHandling = true;
    % vbmc_opts.gpOutwarpFun = @outwarp_negpowc1;
        
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
        if ~isempty(params.vbmc_fits) && numel(params.vbmc_fits.vp) >= iOpt ...
                && ~isempty(params.vbmc_fits.vp{iOpt})
            continue;
        end
        
        [vp,~,~,~,output] = ...
            vbmc(@(x_) -sum(nllfun(x_,params,data))+logprior(x_), ...
            x0,bounds.LB,bounds.UB,bounds.PLB,bounds.PUB,vbmc_opts);

        % Store results
        params.vbmc_fits.vp{iOpt} = vp;
        params.vbmc_fits.output{iOpt} = output;
        
        refitted_flag = true;
        params = compute_vbmc_stats(params);
        
        if save_flag; save_model_fit(data.fullname,params); end
    end    
end

%% Markov Chain Monte Carlo (via slice sampling)
if mcmc_flag
    
    if ~isfield(params,'mcmc_fits') || isempty(params.mcmc_fits) || refit_flags(4)
        params.mcmc_fits = [];
    end
    
    bounds = setup_params([],params);       % Get parameter bounds    
    x0 = params.theta;
    mcmc_fun = @(x_) -sum(nllfun(x_,params,data)); % + logprior(x_);    
    
    Widths = bounds.PUB - bounds.PLB;
    Nsamples = 1e3;
    mcmc_opts.Thin = 3;
    mcmc_opts.Burnin = Nsamples;
    mcmc_opts.Display = 'notify';
    mcmc_opts.Diagnostics = true;
    
    for iOpt = 1:Nvbmc        
        % Skip MCMC run if it already exists
        if ~isempty(params.mcmc_fits) && numel(params.mcmc_fits.samples) >= iOpt ...
                && ~isempty(params.mcmc_fits.samples{iOpt})
            continue;
        end
        
        [samples,fvals,exitflag,output] = ...
            slicesamplebnd(mcmc_fun,x0,Nsamples,Widths,bounds.LB,bounds.UB,mcmc_opts);        

        % Store results
        params.mcmc_fits.samples{iOpt} = samples;
        params.mcmc_fits.output{iOpt} = output;
        
        refitted_flag = true;        
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
function params = compute_hmm_stats(params)

Kmax = numel(params.hmm_fits.hmm);
params.hmm_fits.K = 1:Kmax;
for k = 1:Kmax
    params.hmm_fits.nparams(k) = params.hmm_fits.hmm{k}.nparams;
end
Ntrials = params.mle_fits.ndata;
params.hmm_fits.bic = 2*params.hmm_fits.nll + 0.5*params.hmm_fits.nparams*log(Ntrials);
end

%--------------------------------------------------------------------------
function [hmm,ll,exitflag,output] = compute_hmm1(mle_fits)

% Create first component from MLE fits
hmm = [];
hmm.D = size(mle_fits.x_best,2);
hmm.K = 1;
hmm.mu = mle_fits.x_best;
hmm.pi = 1;
hmm.A = 1;
hmm.lZ = -mle_fits.nll_best;
hmm.ell = [];
hmm.Pk = ones(1,mle_fits.ndata);
hmm.nparams = numel(hmm.mu);

ll = -mle_fits.nll_best;
exitflag = 1;
output = [];

end

%--------------------------------------------------------------------------
function params = compute_vbmc_stats(params)

[exitflag,best,idx_best,stats] = vbmc_diagnostics(params.vbmc_fits.vp);
params.vbmc_fits.diagnostics.exitflag = exitflag;
params.vbmc_fits.diagnostics.best = best;
params.vbmc_fits.diagnostics.idx_best = idx_best;
params.vbmc_fits.diagnostics.stats = stats;

end
function [params,data] = fit_model(model_name,data,Nopts,vbmc_flag,refit_flag,opt_init)
%FIT_MODEL Fit model MODEL_NAME to dataset DATA.

% # restarts for optimization procedure (maximum-likelihood estimation)
if nargin < 3 || isempty(Nopts); Nopts = 10; end

% Get approximate posteriors with Variational Bayesian Monte Carlo
if nargin < 4 || isempty(vbmc_flag); vbmc_flag = false; end

% Force refit even if fit already exists
if nargin < 5 || isempty(refit_flag); refit_flag = false; end

% Starting point(s) for the first fit
if nargin < 6; opt_init = []; end
if isstruct(opt_init); opt_init = {opt_init}; end
if iscell(opt_init)
    for iOpt = 1:numel(opt_init)
        temp(iOpt,:) = opt_init{iOpt}.theta;
    end
    opt_init = temp;
end

if ischar(data)
    data = read_data_from_csv(data);   % Load data
end

% Load existing fits
mypath = fileparts(mfilename('fullpath'));
fits_path = [mypath filesep 'fits'];
addpath(fits_path);
addpath([mypath filesep 'utils']);

matfilename = [data.name '_fits.mat'];
if exist(matfilename,'file')
    fprintf('Found existing file ''%s'', loading previous fits.\n', matfilename);
    for iTry = 1:10
        try
            load(matfilename);
            break;
        catch
            fprintf('Try %d: Load failed, possible I/O error. Waiting a few second before retrying...\n', iTry);
            pause(5 + 5*rand());
        end
    end
else
    modelfits.data = data;
    modelfits.params = [];
end

% Check if model was already fitted to this dataset
params = [];
idx_params = [];
if ~isempty(modelfits.params)
    for iFit = 1:numel(modelfits.params)
        pp = modelfits.params{iFit}; 
        if strcmp(pp.model_name,model_name)
            params = pp;
            idx_params = iFit;
        end
    end
end

%% Maximum likelihood fit
if isempty(params) || refit_flag
    params = params_new(model_name,data);   % Initialize model/params struct
    bounds = setup_params([],params);       % Get parameter bounds

    MaxFunEvals = 2e3;
    x = []; nll = [];
    badopts = bads('defaults');
    badopts.MaxFunEvals = MaxFunEvals;

    % Find starting point from other models
    x0_base = [];
    if ~isempty(params.model_startingfit) && ~isempty(modelfits)
        for iFit = 1:numel(modelfits.params)
            pp = modelfits.params{iFit};
            if strcmp(pp.model_name,params.model_startingfit)
                x0_base = NaN(1,numel(params.names));
                %----------------------------------------------------------
                % This should be assigned based on matching parameter names!
                %----------------------------------------------------------
                x0_base(1:numel(pp.theta)) = pp.theta;  % TO BE FIXED
                fprintf('Reading starting point from model %s.\n', params.model_startingfit);
            end
        end
    end
        
    % Fit the complex change-point models with BADS, otherwise use fmincon
    bads_flag = contains(params.model_name,'changepoint');

    for iOpt = 1:Nopts
        if size(opt_init,1) >= iOpt     % Use provided starting points
            x0 = opt_init(iOpt,:);
        else
            if iOpt == 1
                x0 = bounds.x0;
            else
                x0 = rand(1,numel(bounds.PLB)).*(bounds.PUB-bounds.PLB) + bounds.PLB;
            end
            if ~isempty(x0_base)
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
end

if vbmc_flag
    if ~isfield(params,'vbmc_fit') || isempty(params.vbmc_fit) || refit_flag    
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
        
        vp_array = [];
        for iOpt = 1:ceil(Nopts/2)
            [vp,elbo,elbo_sd,exitflag,output] = ...
                vbmc(@(x_) -nllfun(x_,params,data)+logprior(x_), ...
                x0,bounds.LB,bounds.UB,bounds.PLB,bounds.PUB,vbmc_opts);

            % Store results
            vbmc_fit.vp = vp;
            vbmc_fit.elbo = elbo;
            vbmc_fit.elbo_sd = elbo_sd;
            vbmc_fit.exitflag = exitflag;
            vbmc_fit.output = output;

            params.vbmc_fit(iOpt) = vbmc_fit;
            vp_array{iOpt} = vp;
        end
        
        vbmc_diagnostics(vp_array);
    end
end

end
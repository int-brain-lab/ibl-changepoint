% Plot parameters from fitted changepoint model

vp = [];
mice_list = get_mice_list('charles');
% model_name = 'changepoint_nakarushton_runlength_probs_lapse';
% model_name = 'changepoint_nakarushton_runlength_probs_lapse';
if ~exist('model_name','var') || isempty(model_name)
    model_name = 'changepoint_contrastnoise_runlength_probs';
end
% model_name = 'exponential_contrastnoise';
theta = []; theta_mle = [];
xpdf = [];  ypdf = [];

for iMouse = 1:numel(mice_list)
    fprintf('%s\n',mice_list{iMouse});
    
    params = load_model_fit(mice_list{iMouse},model_name);
    
    if isempty(params)
        theta(iMouse,:) = NaN(1,size(theta,2));        
    else
        if isfield(params,'vbmc_fits') && ~isempty(params.vbmc_fits)
            if isempty(params.vbmc_fits.diagnostics.best)
                fprintf('vbmc did not converge...\n');
                
                stats = params.vbmc_fits.diagnostics.stats;
                elbo = stats.elbo
                elbo_sd = stats.elbo_sd
                beta_lcb = 10;  % Extreme loss aversion
                
                [~,idx] = max(elbo - beta_lcb*elbo_sd);
                vp{iMouse} = params.vbmc_fits.vp{idx};
                
            else
                vp{iMouse} = params.vbmc_fits.diagnostics.best.vp;
            end
            
            exitflag(iMouse) = params.vbmc_fits.diagnostics.exitflag;

            X = get_posterior_samples(params,1e5);
            logflag = false(1,size(X,2));

            switch lower(model_name)
                case 'exponential_contrastnoise'
                    logflag([1 2 5]) = true;
                case 'changepoint_contrastnoise_runlength_probs'
                    logflag([1 2 5 6]) = true;
            end

            bounds = setup_params([],params);
            LB = bounds.LB;
            UB = bounds.UB;

            LB(logflag) = exp(LB(logflag));
            UB(logflag) = exp(UB(logflag));
            X(:,logflag) = exp(X(:,logflag));

            for iParam = 1:size(X,2)
                [~,ypdf{iMouse}(iParam,:),xpdf{iMouse}(iParam,:)] = kde1d(X(:,iParam),2^10,LB(iParam),UB(iParam));
            end            
            
            theta(iMouse,:) = mean(X,1);
            theta(iMouse,logflag) = log(theta(iMouse,logflag));
        end
            
        if isempty(params.mle_fits.nll_best)
            [nLL,idx] = min(params.mle_fits.nll);
            params.theta = params.mle_fits.x(idx,:);
        end
        theta_mle(iMouse,:) = params.theta;
    end
end
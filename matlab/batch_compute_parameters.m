% Plot parameters from fitted changepoint model

vp = [];
mice_list = get_mice_list(5);
% model_name = 'changepoint_nakarushton_runlength_probs_lapse';
% model_name = 'changepoint_nakarushton_runlength_probs_lapse';
if ~exist('model_name','var') || isempty(model_name)
    model_name = 'changepoint_contrastnoise_runlength_probs';
end
% model_name = 'exponential_contrastnoise';
theta = [];
xpdf = [];  ypdf = [];

for iMouse = 1:numel(mice_list)
    fprintf('%s\n',mice_list{iMouse});
    
    params = load_model_fit(mice_list{iMouse},model_name);
    
    if isempty(params)
        theta(iMouse,:) = NaN(1,size(theta,2));        
    else
        if isfield(params,'vbmc_fits') && ~isempty(params.vbmc_fits)
            if isempty(params.vbmc_fits.diagnostics.best)
                fprintf('vbmc did not converge, skipping...\n');
            else            
                vp{iMouse} = params.vbmc_fits.diagnostics.best.vp;
                exitflag(iMouse) = params.vbmc_fits.diagnostics.exitflag;

                X = get_posterior_samples(params,1e5);
                logflag = false(1,size(X,2));
                
                switch lower(model_name)
                    case 'changepoint_contrastnoise_runlength_probs'
                        logflag([5 6]) = true;
                end
                
                bounds = setup_params([],params);
                LB = bounds.LB;
                UB = bounds.UB;
                
                LB(logflag) = exp(LB(logflag));
                UB(logflag) = exp(UB(logflag));
                X(:,logflag) = exp(X(:,logflag));

                for iParam = 1:size(X,2)
                    [~,ypdf{iMouse}(iParam,:),xpdf{iMouse}(iParam,:)] = kde(X(:,iParam),2^10,LB(iParam),UB(iParam));
                end            
            end
            % theta(iMouse,:) = mean(X,1);
        end
            
        if isempty(params.mle_fits.nll_best)
            [nLL,idx] = min(params.mle_fits.nll);
            params.theta = params.mle_fits.x(idx,:);
        end
        theta(iMouse,:) = params.theta;
    end
end
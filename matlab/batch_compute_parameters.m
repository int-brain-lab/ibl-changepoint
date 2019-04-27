% Plot parameters from fitted changepoint model

vp = [];
% mice_list = get_mice_list();
model_name = 'changepoint_doublenoise_loadnoiseendtrain_runlength_probs';
theta = [];

for iMouse = 1:numel(mice_list)
    fprintf('%s\n',mice_list{iMouse});
    
    params = load_model_fit(mice_list{iMouse},model_name);
    
    if isempty(params)
        theta(iMouse,:) = NaN(1,size(theta,2));        
    else
        vp{iMouse} = params.vbmc_fit.diagnostics.best.vp;
        exitflag(iMouse) = params.vbmc_fit.diagnostics.exitflag;

        X = get_posterior_samples(params,1e6); % cornerplot(X)
        theta(iMouse,:) = mean(X,1);
    end
end
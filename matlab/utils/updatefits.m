function updatefits()
%UPDATEFITS Update fits to new version/format. Do not use.

fitfiles = dir('*_fits.mat');
for iFile = 1:numel(fitfiles)
    load(fitfiles(iFile).name); 
    params = modelfits.params;
    data = modelfits.data;
    
    name = fitfiles(iFile).name;
    idx = strfind(name,'_fits.mat');
    data.name = name(1:idx-1);
    
    for iFit = 1:numel(params)
        pp = params{iFit};
        if isfield(pp,'changeprob_model')
            if strcmp(pp.changeprob_model,'bayesian')
                pp.changeprob_model = 'changepoint_runlength';
            end
            
            pp.model_name = pp.changeprob_model;            
            pp = rmfield(pp,'changeprob_model');
            
            temp = params_new(pp.model_name,data);
            
            pp.model_desc = temp.model_desc;
            pp.model_startingfit = temp.model_startingfit;            
        end
        
        if isfield(pp,'runlength_prior')
            pp.runlength_prior = ['@(t) exp(-t/' num2str(pp.runlength_tau,'%.8f') ')'];
        end
        
        if strcmp(func2str(pp.model_nLLfun),'ChangeProb_bocpd_nll_v2')
            pp.model_nLLfun = @changepoint_bayesian_nll;
        end
        
        params{iFit} = pp;
        pp
    end    
    
    modelfits.params = params;
    modelfits.data = data;
    save(fitfiles(iFile).name,'modelfits');
end

end
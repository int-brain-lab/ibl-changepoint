function tab = collect_model_comparison(data_list,model_list)

if nargin < 1; data_list = []; end
if nargin < 2; model_list = []; end

if isempty(data_list); data_list = get_mice_list(); end
if isempty(model_list); model_list = {'omniscient_fixedprior_doublenoise','changepoint_doublenoise','changepoint_doublenoise_runlength_probs'}; end

Ndata = numel(data_list);
Nmodels = numel(model_list);

tab.data = data_list;
tab.models = model_list;

tab.loglike = NaN(Ndata,Nmodels);
tab.aic = NaN(Ndata,Nmodels);
tab.bic = NaN(Ndata,Nmodels);
tab.elbo = NaN(Ndata,Nmodels);
tab.exitflag = NaN(Ndata,Nmodels);


for iData = 1:Ndata
    fprintf('%s ',data_list{iData}); if mod(iData,10) == 0; fprintf('\n'); end
    
    for iModel = 1:Nmodels
        params = load_model_fit(data_list{iData},model_list{iModel});
        if isempty(params); continue; end
        
        Nparams = numel(params.theta);
        
        data1 = read_data_from_csv(data_list{iData});
        Ntrials = size(data1.tab,1);
        
        nLL = nllfun([],params,data1);
        tab.loglike(iData,iModel) = -nLL;
        tab.aic(iData,iModel) = 2*nLL + 2*Nparams;
        tab.bic(iData,iModel) = 2*nLL + log(Ntrials)*Nparams;
        
        if isfield(params,'vbmc_fit')
            if ~isempty(params.vbmc_fit.diagnostics.best)
                tab.elbo(iData,iModel) = params.vbmc_fit.diagnostics.best.elbo;
            end
            tab.exitflag(iData,iModel) = params.vbmc_fit.diagnostics.exitflag;
        end
    end
end

end


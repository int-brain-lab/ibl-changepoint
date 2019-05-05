function tab = collect_model_comparison(data_list,model_list)

if nargin < 1; data_list = []; end
if nargin < 2; model_list = []; end

if isempty(data_list); data_list = get_mice_list(); end
if isempty(model_list)
%    model_list = {'omniscient_fixedprior_doublenoise','changepoint_doublenoise','changepoint_doublenoise_runlength_probs'}; 
    model_list = {'omniscient_fixedprior_nakarushton','changepoint_nakarushton','changepoint_nakarushton_runlength_probs', ...
        'omniscient_fixedprior_nakarushton_lapse','changepoint_nakarushton_lapse','changepoint_nakarushton_runlength_probs_lapse'};
%        'omniscient_fixedprior_doublenoise','changepoint_doublenoise','changepoint_doublenoise_runlength_probs'}; 
end

Ndata = numel(data_list);
Nmodels = numel(model_list);

tab.data = data_list;
tab.models = model_list;

tab_fields = {'nopts','loglike','aic','bic','nvps','elbo','exitflag'};

for iField = 1:numel(tab_fields)
    tab.(tab_fields{iField}) = NaN(Ndata,Nmodels);
end

for iData = 1:Ndata
    fprintf('%s ',data_list{iData}); if mod(iData,10) == 0; fprintf('\n'); end
    
    for iModel = 1:Nmodels
        params = load_model_fit(data_list{iData},model_list{iModel});
        if isempty(params); continue; end
                
        data1 = read_data_from_csv(data_list{iData});
        Ntrials = size(data1.tab,1);
                        
        nLL = params.mle_fits.nll_best;
        if isempty(nLL)
            [nLL,idx] = min(params.mle_fits.nll);
            params.theta = params.mle_fits.x(idx,:);
        end
                
        Nparams = numel(params.theta);
        
        tab.nopts(iData,iModel) = size(params.mle_fits.x,1);
        tab.loglike(iData,iModel) = -nLL;
        tab.aic(iData,iModel) = 2*nLL + 2*Nparams;
        tab.bic(iData,iModel) = 2*nLL + log(Ntrials)*Nparams;
        
        if isfield(params,'vbmc_fits')
            if ~isempty(params.vbmc_fits.diagnostics.best)
                tab.elbo(iData,iModel) = params.vbmc_fits.diagnostics.best.elbo;
            end
            tab.exitflag(iData,iModel) = params.vbmc_fits.diagnostics.exitflag;
            tab.nvps(iData,iModel) = numel(params.vbmc_fits.vps);
        end
    end
end

end


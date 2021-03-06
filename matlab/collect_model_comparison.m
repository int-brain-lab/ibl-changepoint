function tab = collect_model_comparison(data_list,model_list)

if nargin < 1; data_list = []; end
if nargin < 2; model_list = []; end

if isempty(data_list); data_list = get_mice_list(5); end
if isempty(model_list)
%     model_list = { ...
%         'omniscient_fixedprior_nakarushton','changepoint_nakarushton','changepoint_nakarushton_runlength_probs', ...
%         'exponential_nakarushton', 'exponential_hyperprobs_nakarushton', ...
%         'omniscient_fixedprior_nakarushton_lapse','changepoint_nakarushton_lapse','changepoint_nakarushton_runlength_probs_lapse', ...
%         'exponential_nakarushton_lapse', 'exponential_hyperprobs_nakarushton_lapse' ...
%         'omniscient_fixedprior_nakarushton_marginalize_approx','changepoint_nakarushton_marginalize_approx','changepoint_nakarushton_runlength_probs_marginalize_approx', ...
%         'omniscient_fixedprior_nakarushton_marginalize_approx_lapse','changepoint_nakarushton_marginalize_approx_lapse','changepoint_nakarushton_runlength_probs_marginalize_approx_lapse', ...
%         };
%    model_list = { ...
%        'changepoint_contrastnoise',
%        'changepoint_contrastnoise_runlength_probs', 'exponential_contrastnoise', ...
%		'omniscient_contrastnoise_fixedfreeprior', 'exponential_contrastnoise_nobias'...    
%    }
    model_list = { ...
        'changepoint_contrastnoise', 'changepoint_contrastnoise_runlength_probs', 'exponential_contrastnoise', ...
		'omniscient_contrastnoise_fixedfreeprior' ...    
    }
end

Ndata = numel(data_list);
Nmodels = numel(model_list);

tab.data = data_list;
tab.models = model_list;

tab_fields = {'ntrials','nopts','loglike','cvll','aic','bic','nvps','elbo','exitflag'};

for iField = 1:numel(tab_fields)
    tab.(tab_fields{iField}) = NaN(Ndata,Nmodels);
end

for iData = 1:Ndata
    fprintf('%s ',data_list{iData}); if mod(iData,10) == 0; fprintf('\n'); end
    
    for iModel = 1:Nmodels
        params = load_model_fit(data_list{iData},model_list{iModel});
        if isempty(params); continue; end
                
        if isfield(params.mle_fits,'ndata')
            Ntrials = params.mle_fits.ndata;
        else
            data1 = read_data_from_csv(data_list{iData});
            Ntrials = size(data1.tab,1);
        end
        
        nLL = params.mle_fits.nll_best;
        if isempty(nLL)
            [nLL,idx] = min(params.mle_fits.nll);
            params.theta = params.mle_fits.x(idx,:);
        end
                
        Nparams = numel(params.theta);
        
        tab.ntrials(iData,iModel) = params.mle_fits.ndata;
        tab.nopts(iData,iModel) = size(params.mle_fits.x,1);
        tab.loglike(iData,iModel) = -nLL;
        tab.aic(iData,iModel) = 2*nLL + 2*Nparams;
        tab.bic(iData,iModel) = 2*nLL + log(Ntrials)*Nparams;

        if isfield(params.mle_fits,'cv') && ~isempty(params.mle_fits.cv)
            tab.cvll(iData,iModel) = -sum(params.mle_fits.cv.nll_test);            
        end
        
        if isfield(params,'vbmc_fits')
            if ~isempty(params.vbmc_fits.diagnostics.best)
                tab.elbo(iData,iModel) = params.vbmc_fits.diagnostics.best.elbo;
            end
            tab.exitflag(iData,iModel) = params.vbmc_fits.diagnostics.exitflag;
            tab.nvps(iData,iModel) = numel(params.vbmc_fits.vp);
        end
    end
end

end


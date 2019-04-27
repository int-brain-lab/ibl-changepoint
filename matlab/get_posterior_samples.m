function X = get_posterior_samples(params,N)
%GET_POSTERIOR_SAMPLES Return samples from (approximate) posterior.

X = [];

if isfield(params,'vbmc_fit') && ~isempty(params.vbmc_fit)
    vp = params.vbmc_fit.diagnostics.best.vp;    
    X = vbmc_rnd(vp,N);
end


end
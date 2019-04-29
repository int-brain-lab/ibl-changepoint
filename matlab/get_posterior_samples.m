function X = get_posterior_samples(params,N)
%GET_POSTERIOR_SAMPLES Return samples from (approximate) posterior.

X = [];

if isfield(params,'vbmc_fits') && ~isempty(params.vbmc_fits)
    vp = params.vbmc_fits.diagnostics.best.vp;    
    X = vbmc_rnd(vp,N);
end


end
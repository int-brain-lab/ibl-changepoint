function X = get_posterior_samples(params,N,beta_lcb)
%GET_POSTERIOR_SAMPLES Return samples from (approximate) posterior.

if nargin < 3 || isempty(beta_lcb); beta_lcb = 3; end

X = [];

if isfield(params,'vbmc_fit') && ~isempty(params.vbmc_fit)
    vp = vbmc_fit.diagnostics.best.vp;    
    X = vbmc_rnd(vp,N);
end


end
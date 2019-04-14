function X = get_posterior_samples(params,N,beta_lcb)
%GET_POSTERIOR_SAMPLES Return samples from (approximate) posterior.

if nargin < 3 || isempty(beta_lcb); beta_lcb = 3; end

X = [];

if isfield(params,'vbmc_fit') && ~isempty(params.vbmc_fit)
    vbmc_fit = params.vbmc_fit;
    
    for iFit = 1:numel(vbmc_fit)
        elbo(iFit) = vbmc_fit(iFit).elbo;
        elbo_sd(iFit) = vbmc_fit(iFit).elbo_sd;
        exitflag(iFit) = vbmc_fit(iFit).exitflag;
    end

    idx_ok = exitflag == 1;
    
    if sum(idx_ok) == 0
        idx_ok = true(size(idx_ok)); 
        warning('No solution has converged, using potentially unstable solution.');
    end
    
    elcbo = elbo - beta_lcb*elbo_sd;
    elcbo(~idx_ok) = -Inf;
    [~,best_idx] = max(elcbo);
    
    vp = vbmc_fit(best_idx).vp;    
    X = vbmc_rnd(vp,N);
end


end
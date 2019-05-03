function [nLL,output] = omniscient_bayesian_nll(params,data)

% Compute prior
if isfield(params,'fixed_prior') && ~isempty(params.fixed_prior)
    priorL = params.fixed_prior*ones(size(data.p_true));
else
    priorL = data.p_true;
end

% Pre-compute response probability as a function of signed contrast level
% and log prior odds for speed
pgrid = sort(unique(priorL));
pgrid = [pgrid(1)-sqrt(eps),pgrid(:)',pgrid(end)+sqrt(eps)];
[params.PChatL_grid,params.lp_odds_grid] = ...
    precompute_sdt(params,data,pgrid);

% Compute negative log likelihood and probability of responding L
[nLL,PChatL] = sdt_nll(params,data,priorL);

if nargout > 1
    output.resp_model = PChatL;    
end


end
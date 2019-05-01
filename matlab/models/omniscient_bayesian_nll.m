function [nLL,output] = omniscient_bayesian_nll(params,data,mu,sigma)

% Potentially different SIGMA for left and right stimuli, copy if only one
% is passed
if size(sigma,2) == 1; sigma = repmat(sigma,[1,2]); end

% Compute prior
if isfield(params,'fixed_prior') && ~isempty(params.fixed_prior)
    priorL = params.fixed_prior*ones(size(data.p_true));
else
    priorL = data.p_true;
end

% Compute negative log likelihood and probability of responding L
[nLL,PChatL] = sdt_nll(data.S,mu,sigma,priorL,data.resp_obs,params);

if nargout > 1
    output.resp_model = PChatL;    
end


end
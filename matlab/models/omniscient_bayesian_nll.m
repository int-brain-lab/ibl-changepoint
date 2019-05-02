function [nLL,output] = omniscient_bayesian_nll(params,data)

% Compute prior
if isfield(params,'fixed_prior') && ~isempty(params.fixed_prior)
    priorL = params.fixed_prior*ones(size(data.p_true));
else
    priorL = data.p_true;
end

% Compute negative log likelihood and probability of responding L
[nLL,PChatL] = sdt_nll(params,data,priorL);

if nargout > 1
    output.resp_model = PChatL;    
end


end
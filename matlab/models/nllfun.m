function [nLL,output] = nllfun(theta,params,data)
%NLLFUN Negative log likelihood function

if isempty(theta)
    theta = params.theta;
end

Nparams = size(theta,1);
nLL = NaN(1,Nparams);

for iParam = 1:Nparams
    % Assign parameter vector to parameter struct
    params1 = setup_params(theta(iParam,:),params);
        
    if nargout == 1
        nLL(iParam) = params1.model_nLLfun(params1,data);
    else
        [nLL(iParam),output(iParam)] = params1.model_nLLfun(params1,data);
    end
end

end
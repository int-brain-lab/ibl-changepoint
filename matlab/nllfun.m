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

    % Compute SIGMA values for each contrast level
    if isfield(params1,'sigma_poly')
        sigma_vec = polyval(params1.sigma_poly,log(data.contrasts_vec));
        sigma_vec(data.contrasts_vec == 0) = Inf;
        sigma_vec = min(max(sigma_vec,1),360);
        sigma = sigma_vec(data.contrasts_idx);

        if isfield(params1,'attention_factor') && params1.attention_factor ~= 1
            sigma = repmat(sigma,[1 2]);
            sigma(:,1) = sigma(:,1)*params1.attention_factor;
            sigma(:,2) = sigma(:,2)/params1.attention_factor;
        end
    elseif isfield(params1,'precision')
       sigma_vec = 1./sqrt(params1.precision(1)*(data.contrasts_vec.^params1.precision_power(1)));
       sigma_vec(data.contrasts_vec == 0) = Inf;
       sigma_vec = min(max(sigma_vec,1),360);
       sigma(:,1) = sigma_vec(data.contrasts_idx);
       
       sigma_vec = 1./sqrt(params1.precision(2)*(data.contrasts_vec.^params1.precision_power(2)));
       sigma_vec(data.contrasts_vec == 0) = Inf;
       sigma_vec = min(max(sigma_vec,1),360);
       sigma(:,2) = sigma_vec(data.contrasts_idx);
    else
        sigma = [];
    end

    %----------------------------------------------------------------------
    % Temporary patch for added parameter
    %----------------------------------------------------------------------
    if ~isfield(params1,'lapse_bias')
        params1.lapse_bias = 0.5;
    end
    
    if nargout == 1
        nLL(iParam) = params1.model_nLLfun(params1,data,sigma);
    else
        [nLL(iParam),output(iParam)] = params1.model_nLLfun(params1,data,sigma);            
    end
end

end

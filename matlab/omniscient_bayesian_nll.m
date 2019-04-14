function [nLL,output] = omniscient_bayesian_nll(params,data,sigma)

MIN_P = 1e-4;   % Minimum lapse/error probability

%% Assign observer model parameters

lambda = max(MIN_P,params.lambda);     % Minimum lapse to avoid numerical trouble
lapse_bias = params.lapse_bias;
softmax_eta = params.softmax_eta;
softmax_bias = params.softmax_bias;

% Potentially different SIGMA for left and right stimuli, copy if only one
% is passed
if size(sigma,2) == 1; sigma = repmat(sigma,[1,2]); end

% Get measurement noise grid and pdf
[X,W] = get_noisy_measurements(data.S,sigma);

%% Compute log posterior for each measurement

loglikeL = bsxfun(@plus,-0.5*bsxfun(@rdivide,X-data.mu(1),sigma(:,1)).^2, -log(sigma(:,1)));
loglikeR = bsxfun(@plus,-0.5*bsxfun(@rdivide,X-data.mu(2),sigma(:,2)).^2, -log(sigma(:,2)));
logprior_odds = log(data.p_true./(1-data.p_true));

% Decision variable (>0 for L, <0 for R)
dhat = bsxfun(@plus, loglikeL - loglikeR, logprior_odds);

% Compute negative log likelihood and probability of responding L
[nLL,PChatL] = get_responses_nLL(dhat,softmax_eta,softmax_bias,lambda,lapse_bias,W,data.resp_obs);

if nargout > 1
    output.resp_model = PChatL;    
end


end
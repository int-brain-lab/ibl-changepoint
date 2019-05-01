function [nLL,PChatL] = sdt_nll(S,mu,sigma,priorL,resp_obs,params,nx)
%SDT_NLL Compute negative log likelihood with signal-detection theory model.

% Measurement grid size
if nargin < 7 || isempty(nx); nx = 1001; end

%% 1. Get noisy measurement grid and pdf

MAXSD = 5;  % When integrating a Gaussian go up to this distance
S_cat = (S >= 0) + 1;   % Category is 1 for L and 2 for R

% Extract 1st MU/SIGMA for L stimulus, 2nd MU/SIGMA for R stimulus
mu_vec = mu(sub2ind(size(mu),(1:size(mu,1))',S_cat));
sigma_vec = sigma(sub2ind(size(sigma),(1:size(sigma,1))',S_cat));

X = bsxfun(@plus, mu_vec, bsxfun(@times, sigma_vec, MAXSD*linspace(-1,1,nx)));
W = normpdf(linspace(-MAXSD,MAXSD,nx));
W = W./qtrapz(W);

%% 2. Compute decision variable according to noise model

loglikeL = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,mu(:,1)),sigma(:,1)).^2, -log(sigma(:,1)));
loglikeR = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,mu(:,2)),sigma(:,2)).^2, -log(sigma(:,2)));

logprior_odds = log(priorL./(1-priorL));

% Decision variable (>0 for L, <0 for R)
dhat = bsxfun(@plus, loglikeL - loglikeR, logprior_odds);

% Ntrials = size(X,1);
% iter = 0;
% idx_max = 0;
% stride = 2e4;
% 
% while idx_max < Ntrials
%     idx_min = 1+stride*iter;
%     idx_max = min(stride*(iter+1),Ntrials);
%     idx = (idx_min:idx_max)';
%     
%     loglikeL = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X(idx,:),mu(idx,1)),sigma(idx,1)).^2, -log(sigma(idx,1)));
%     loglikeR = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X(idx,:),mu(idx,2)),sigma(idx,2)).^2, -log(sigma(idx,2)));
% 
%     logprior_odds = log(priorL(idx)./(1-priorL(idx)));
% 
%     % Decision variable (>0 for L, <0 for R)
%     dhat(idx,:) = bsxfun(@plus, loglikeL - loglikeR, logprior_odds);
%     
%     iter = iter + 1;
% end

%% 3. Compute negative log likelihood of responses

MIN_P = 1e-4;   % Minimum lapse/error probability

lapse_rate = max(MIN_P,params.lapse_rate);     % Minimum lapse to avoid numerical trouble
lapse_bias = params.lapse_bias;
softmax_eta = params.softmax_eta;
softmax_bias = params.softmax_bias;

% Compute probability given decision variable DHAT
PCx_soft = 1./(1+exp(-softmax_eta*(dhat + softmax_bias)));

% Marginalize over noisy measurements and add lapse
PChatL = qtrapz(bsxfun(@times,PCx_soft,W),2);
PChatL = lapse_rate*lapse_bias + (1-lapse_rate)*PChatL;

% Log probability of observer responses (ignore no-response trials)
log_PChat = log(PChatL).*(resp_obs == -1) + log(1-PChatL).*(resp_obs == 1);

% Correct for NaNs or Infs
log_PChat(~isfinite(log_PChat)) = log(MIN_P);

% Sum negative log likelihood
nLL = -nansum(log_PChat);

end
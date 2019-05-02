function [nLL,PChatL] = sdt_nll(params,data,priorL,nx)
%SDT_NLL Compute negative log likelihood with signal-detection theory model.

% Measurement grid size
if nargin < 4 || isempty(nx); nx = 201; end


%% 1. Signal-detection theory noise model setup

contrasts_vec = data.contrasts_vec;

if isfield(params,'nakarushton_response_min')
    % Naka-Rushton sensory noise model
    response_min = params.nakarushton_response_min;
    response_delta = params.nakarushton_response_delta;
    n = params.nakarushton_n;
    c50 = params.nakarushton_c50;
    neff_left = params.nakarushton_neff_left;
    neff_right = params.nakarushton_neff_right;
    
    % Naka-Rushton contrast response curve (same for left and right)
    r_vec = response_min + response_delta./(1 + (c50./contrasts_vec).^n);
    
    mu_vec(1,:) = r_vec - response_min; % Stimulus left
    mu_vec(2,:) = response_min - r_vec; % Stimulus right
    
    sigma_vec(1,:) = sqrt(r_vec/neff_left + response_min/neff_right);
    sigma_vec(2,:) = sqrt(response_min/neff_left + r_vec/neff_right);

else
    mu_vec = repmat(data.mu(:),[1,numel(contrasts_vec)]);

    % Compute SIGMA values for each contrast level    
    if isfield(params,'sigma_poly')
        sigma_vec(1,:) = poly2sigmavec(params.sigma_poly,contrasts_vec);
        sigma_vec = repmat(sigma_vec,[2,1]);
        if isfield(params,'attention_factor') && params.attention_factor ~= 1
            sigma_vec(1,:) = sigma_vec(1,:)*params.attention_factor;
            sigma_vec(2,:) = sigma_vec(2,:)/params.attention_factor;
        end
    elseif isfield(params,'sigma_poly_left')
        sigma_vec(1,:) = poly2sigmavec(params.sigma_poly_left,contrasts_vec);
        sigma_vec(2,:) = poly2sigmavec(params.sigma_poly_right,contrasts_vec);
    else
        sigma = [];
    end    
end

mu(:,1) = mu_vec(1,data.contrasts_idx);
mu(:,2) = mu_vec(2,data.contrasts_idx);

sigma(:,1) = sigma_vec(1,data.contrasts_idx);
sigma(:,2) = sigma_vec(2,data.contrasts_idx);

%% 2. Get noisy measurement grid and pdf

MAXSD = 5;  % When integrating a Gaussian go up to this distance
S_cat = (data.S >= 0) + 1;   % Category is 1 for L and 2 for R

% Extract 1st MU/SIGMA for L stimulus, 2nd MU/SIGMA for R stimulus
mu_t = mu(sub2ind(size(mu),(1:size(mu,1))',S_cat));
sigma_t = sigma(sub2ind(size(sigma),(1:size(sigma,1))',S_cat));

X = bsxfun(@plus, mu_t, bsxfun(@times, sigma_t, MAXSD*linspace(-1,1,nx)));
W = normpdf(linspace(-MAXSD,MAXSD,nx));
W = W./qtrapz(W);

%% 3. Compute decision variable according to noise model

if isfield(params,'marginalize_contrasts') && params.marginalize_contrasts
    % Marginalize over non-zero contrasts (assumes uniform prior over contrasts)
    
    % Ignore zero-contrast
    muL_vec3(1,1,:) = mu_vec(1,2:end);
    sigmaL_vec3(1,1,:) = sigma_vec(1,2:end);
    muR_vec3(1,1,:) = mu_vec(2,2:end);
    sigmaR_vec3(1,1,:) = sigma_vec(2,2:end);
    
    % Sum over likelihoods (with numerical savviness to avoid overflows)
    loglikeL_contrast = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,muL_vec3),sigmaL_vec3).^2, -log(sigmaL_vec3));
    maxL = max(loglikeL_contrast,[],3);
    loglikeL = log(sum(exp(bsxfun(@minus,loglikeL_contrast,maxL)),3)) + maxL;
    
    loglikeR_contrast = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,muR_vec3),sigmaR_vec3).^2, -log(sigmaR_vec3));
    maxR = max(loglikeR_contrast,[],3);
    loglikeR = log(sum(exp(bsxfun(@minus,loglikeR_contrast,maxR)),3)) + maxL;
else    
    loglikeL = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,mu(:,1)),sigma(:,1)).^2, -log(sigma(:,1)));
    loglikeR = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,mu(:,2)),sigma(:,2)).^2, -log(sigma(:,2)));
end

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

%% 4. Compute negative log likelihood of responses

MIN_P = 1e-4;   % Minimum lapse/error probability

lapse_rate = max(MIN_P,params.lapse_rate);     % Minimum lapse to avoid numerical trouble
lapse_bias = params.lapse_bias;
softmax_eta = params.softmax_eta;
softmax_bias = params.softmax_bias;

% Compute probability given decision variable DHAT
PCx_soft = 1./(1+exp(-softmax_eta*(dhat + softmax_bias)));
PCx_soft(dhat == 0) = 0.5;

% Marginalize over noisy measurements and add lapse
PChatL = qtrapz(bsxfun(@times,PCx_soft,W),2);
PChatL = lapse_rate*lapse_bias + (1-lapse_rate)*PChatL;

% Log probability of observer responses (ignore no-response trials)
log_PChat = log(PChatL).*(data.resp_obs == -1) + log(1-PChatL).*(data.resp_obs == 1);

% Correct for NaNs or Infs
log_PChat(~isfinite(log_PChat)) = log(MIN_P);

% Sum negative log likelihood
nLL = -sum(log_PChat);

end

%--------------------------------------------------------------------------
function sigma_vec = poly2sigmavec(sigma_poly,contrasts_vec)
%POLY2SIGMAVEC Evaluates polynomial to SIGMA values

sigma_vec = polyval(sigma_poly,log(contrasts_vec));
sigma_vec(contrasts_vec == 0) = Inf;
sigma_vec(~isfinite(sigma_vec)) = Inf;
sigma_vec = min(max(sigma_vec,1),360);
        
end

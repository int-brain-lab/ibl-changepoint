function [PChatL,logprior_odds,mu_sc,sigma_sc] = precompute_sdt(params,data,pgrid,nx)
%PRECOMPUTE_SDT Precompute response probability as a function of log prior 
% odds and log prior odds for the signal-detection theory noise model

if nargin < 4; nx = []; end

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

% MU and SIGMA of Gaussians in the SDT model by signed contrast level
mu_sc = [fliplr(mu_vec(1,2:end)),mu_vec(2,:)];
sigma_sc = [fliplr(sigma_vec(1,2:end)),sigma_vec(2,:)];


%% 2. Get noisy measurement grid and pdf

if isempty(nx); nx = 5001; end

if 1
    MAXSD = 10;
    X = bsxfun(@plus, mu_sc(:), bsxfun(@times, sigma_sc(:), MAXSD*linspace(-1,1,nx)));
    W = normpdf(linspace(-MAXSD,MAXSD,nx));
    W = W./qtrapz(W);
else
    MAXSD = 6;
    Nc = numel(mu_sc);  nx = ceil(nx/Nc);
    X = bsxfun(@plus, mu_sc(:), bsxfun(@times, sigma_sc(:), MAXSD*linspace(-1,1,nx)));
    X = repmat(sort(X(:)'),[Nc,1]);
    
    W = normpdf(X,mu_sc(:),sigma_sc(:));
    dx = 0.5*[diff(X(1,:)),0] + 0.5*[0,diff(X(1,:))];
    W = W.*dx;    
    W = W./sum(W,2);
    
end

%% 3. Compute decision variable according to noise model

Ncontrasts = numel(data.contrasts_vec);

if isfield(params,'marginalize_contrasts') && params.marginalize_contrasts
    % Marginalize over non-zero contrasts (assumes uniform prior over contrasts)

    % Ignore zero-contrast
    muL_vec3(1,1,:) = mu_sc(1:Ncontrasts-1);
    sigmaL_vec3(1,1,:) = sigma_vec(1:Ncontrasts-1);
    muR_vec3(1,1,:) = mu_sc(Ncontrasts+1:end);
    sigmaR_vec3(1,1,:) = sigma_sc(Ncontrasts+1:end);

    % Sum over likelihoods (with numerical savviness to avoid overflows)
    loglikeL_contrast = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,muL_vec3),sigmaL_vec3).^2, -log(sigmaL_vec3));
    maxL = max(loglikeL_contrast,[],3);
    loglikeL = log(sum(exp(bsxfun(@minus,loglikeL_contrast,maxL)),3)) + maxL;

    loglikeR_contrast = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,muR_vec3),sigmaR_vec3).^2, -log(sigmaR_vec3));
    maxR = max(loglikeR_contrast,[],3);
    loglikeR = log(sum(exp(bsxfun(@minus,loglikeR_contrast,maxR)),3)) + maxL;
else    
    muL(:,1) = [mu_sc(1:Ncontrasts),fliplr(mu_sc(1:Ncontrasts-1))];
    sigmaL(:,1) = [sigma_sc(1:Ncontrasts),fliplr(sigma_sc(1:Ncontrasts-1))];
    muR(:,1) = [fliplr(mu_sc(Ncontrasts+1:end)),mu_sc(Ncontrasts:end)];
    sigmaR(:,1) = [fliplr(sigma_sc(1:Ncontrasts)),sigma_sc(Ncontrasts:end)];

    loglikeL = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,muL),sigmaL).^2, -log(sigmaL));
    loglikeR = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,muR),sigmaR).^2, -log(sigmaR));
end    

% Compute log prior odds for a grid of provided P(Left) values
logprior_odds(1,1,:) = log(pgrid./(1-pgrid));

% Decision variable (3-D table) for contrast level, measurement, and log prior odds
dhat = bsxfun(@plus, loglikeL - loglikeR, logprior_odds);

%% 4. Compute negative log likelihood of responses

softmax_eta = params.softmax_eta;
softmax_bias = params.softmax_bias;

% Compute probability given decision variable DHAT
PCx_soft = 1./(1+exp(-softmax_eta*(dhat + softmax_bias)));
PCx_soft(dhat == 0) = 0.5;

% Marginalize over noisy measurements
PChatL(:,:) = qtrapz(bsxfun(@times,PCx_soft,W),2);

% PCHATL is a 2-D matrix representing P(resp = Left) for contrast level 
% (rows) times log prior odds (columns)

logprior_odds = logprior_odds(:)';

end

%--------------------------------------------------------------------------
function sigma_vec = poly2sigmavec(sigma_poly,contrasts_vec)
%POLY2SIGMAVEC Evaluates polynomial to SIGMA values

sigma_vec = polyval(sigma_poly,log(contrasts_vec));
sigma_vec(contrasts_vec == 0) = Inf;
sigma_vec(~isfinite(sigma_vec)) = Inf;
sigma_vec = min(max(sigma_vec,1),360);
        
end

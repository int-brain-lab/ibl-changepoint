function [PChatL,logprior_odds,mu_sc,sigma_sc] = precompute_sdt(params,data,pgrid,nx)
%PRECOMPUTE_SDT Precompute response probability as a function of log prior 
% odds and log prior odds for the signal-detection theory noise model

if nargin < 4; nx = []; end

if all(isnan(data.resp_obs))
    % Skip pre-computation of the log-likelihood
    PChatL = NaN;
    logprior_odds = NaN;
    mu_sc = NaN;
    sigma_sc = NaN;
    return;
end

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

elseif isfield(params,'contrast_sigma')
    
    Nc = numel(contrasts_vec);
    
    mu_vec(1,:) = -contrasts_vec;
    mu_vec(2,:) = contrasts_vec;
    
    sigma_vec(1,:) = params.contrast_sigma(1)*ones(1,Nc);
    sigma_vec(2,:) = params.contrast_sigma(2)*ones(1,Nc);
    
    % Adjust zero-contrast to be equal
    sigma_vec(mu_vec == 0) = sqrt(0.5*(params.contrast_sigma(1)^2 + params.contrast_sigma(2)^2));
    
    nf_vec = normcdf(1,mu_vec,sigma_vec) - normcdf(-1,mu_vec,sigma_vec);
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

if isempty(nx); nx = 2001; end

if isfield(params,'contrast_sigma')
    
    X = linspace(-1,1,nx);
    L = X(end)-X(1);
    dx = X(2)-X(1);
    
    contrast_eps = [params.contrast_epsilon(1)*ones(1,Nc),params.contrast_epsilon(2)*ones(1,Nc)];
    contrast_eps([mu_vec(1,:),mu_vec(2,:)]' == 0) = mean(params.contrast_epsilon);
    
    W = bsxfun(@rdivide, bsxfun_normpdf(X,[mu_vec(1,:),mu_vec(2,:)]',[sigma_vec(1,:),sigma_vec(2,:)]'), [nf_vec(1,:),nf_vec(2,:)]');
    W = bsxfun(@plus, bsxfun(@times, W, 1-contrast_eps'), contrast_eps'/L)*dx;
    
    
else

    if 1
        randomize = false;
        MAXSD = 10;
        xx = MAXSD*linspace(-1,1,nx);
        if randomize
            shift = 2*(rand()-0.5)/(nx-1)*MAXSD;
            xx = xx - shift;
        end
        X = bsxfun(@plus, mu_sc(:), bsxfun(@times, sigma_sc(:), xx));
        W = normpdf(xx);
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
end

%% 3. Compute decision variable according to noise model

Ncontrasts = numel(data.contrasts_vec);

marginalize_contrasts = false;
approx_marginalization = false;
if isfield(params,'marginalize_contrasts'); marginalize_contrasts = params.marginalize_contrasts; end
if isfield(params,'marginalize_approx'); approx_marginalization = params.marginalize_approx; end

if isfield(params,'contrast_sigma')
    
    muL_vec3(1,1,:) = mu_vec(1,:);
    sigmaL_vec3(1,1,:) = sigma_vec(1,:);
    nfL_vec3(1,1,:) = nf_vec(1,:);
    muR_vec3(1,1,:) = mu_vec(2,:);
    sigmaR_vec3(1,1,:) = sigma_vec(2,:);
    nfR_vec3(1,1,:) = nf_vec(2,:);
        
    loglikeL = log(mean(bsxfun(@rdivide, bsxfun_normpdf(X,muL_vec3,sigmaL_vec3),nfL_vec3),3)*(1 - contrast_eps(1)) + contrast_eps(1)/L);
    loglikeR = log(mean(bsxfun(@rdivide, bsxfun_normpdf(X,muR_vec3,sigmaR_vec3),nfR_vec3),3)*(1 - contrast_eps(end)) + contrast_eps(end)/L);
    
elseif marginalize_contrasts || approx_marginalization
    % Marginalize over non-zero contrasts (assumes uniform prior over contrasts)
    
    % Ignore zero-contrast
    muL_vec3(1,1,:) = mu_sc(1:Ncontrasts-1);
    sigmaL_vec3(1,1,:) = sigma_vec(1:Ncontrasts-1);
    muR_vec3(1,1,:) = mu_sc(Ncontrasts+1:end);
    sigmaR_vec3(1,1,:) = sigma_sc(Ncontrasts+1:end);

    if approx_marginalization
        % Approximate marginalized likelihood with single Gaussian
        mu_approx = mean(muL_vec3);
        sigma2_approx = mean(muL_vec3.^2 + sigmaL_vec3.^2) - mu_approx^2;
        loglikeL = -0.5*(X-mu_approx).^2./sigma2_approx - 0.5*log(2*pi*sigma2_approx);
        
        mu_approx = mean(muR_vec3);
        sigma2_approx = mean(muR_vec3.^2 + sigmaR_vec3.^2) - mu_approx^2;
        loglikeR = -0.5*(X-mu_approx).^2./sigma2_approx - 0.5*log(2*pi*sigma2_approx);        
    else
        % Sum over likelihoods (with numerical savviness to avoid overflows)
        loglikeL_contrast = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,muL_vec3),sigmaL_vec3).^2, -log(sigmaL_vec3));
        maxL = max(loglikeL_contrast,[],3);
        loglikeL = log(sum(exp(bsxfun(@minus,loglikeL_contrast,maxL)),3)) + maxL;
        
        loglikeR_contrast = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,muR_vec3),sigmaR_vec3).^2, -log(sigmaR_vec3));
        maxR = max(loglikeR_contrast,[],3);
        loglikeR = log(sum(exp(bsxfun(@minus,loglikeR_contrast,maxR)),3)) + maxL;
    end

else    
    muL(:,1) = [mu_sc(1:Ncontrasts),fliplr(mu_sc(1:Ncontrasts-1))];
    sigmaL(:,1) = [sigma_sc(1:Ncontrasts),fliplr(sigma_sc(1:Ncontrasts-1))];
    muR(:,1) = [fliplr(mu_sc(Ncontrasts+1:end)),mu_sc(Ncontrasts:end)];
    sigmaR(:,1) = [fliplr(sigma_sc(Ncontrasts+1:end)),sigma_sc(Ncontrasts:end)];

    loglikeL = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,muL),sigmaL).^2, -log(sigmaL));
    loglikeR = bsxfun(@plus,-0.5*bsxfun(@rdivide,bsxfun(@minus,X,muR),sigmaR).^2, -log(sigmaR));
end    

% Compute log prior odds for a grid of provided P(Left) values
logprior_odds(1,:) = log(pgrid./(1-pgrid));

[Nc,Nx] = size(loglikeL);
Np = numel(logprior_odds);

loglikediff_t = (loglikeL - loglikeR)';

% Compute probability of responding Left
if isfield(params,'contrast_sigma')
    PChatL(:,:) = compute_pchatL(Nc,Nx,Np,loglikediff_t,logprior_odds,...
        params.softmax_eta,params.softmax_bias,W);    
else
    PChatL(:,:) = compute_pchatL_mex(Nc,Nx,Np,loglikediff_t,logprior_odds,...
        params.softmax_eta,params.softmax_bias,W);
end

% PCHATL is a 2-D matrix representing P(resp = Left) for contrast level 
% (rows) times log prior odds (columns)


end

%--------------------------------------------------------------------------
function sigma_vec = poly2sigmavec(sigma_poly,contrasts_vec)
%POLY2SIGMAVEC Evaluates polynomial to SIGMA values

sigma_vec = polyval(sigma_poly,log(contrasts_vec));
sigma_vec(contrasts_vec == 0) = Inf;
sigma_vec(~isfinite(sigma_vec)) = Inf;
sigma_vec = min(max(sigma_vec,1),360);
        
end

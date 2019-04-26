function [X,W] = get_noisy_measurements(S,sigma,nx)
%GET_NOISY_MEASUREMENTS Return grid of noisy measurements and their pdf.

% Measurement grid size
if nargin < 3 || isempty(nx); nx = 1001; end

MAXSD = 5;  % When integrating a Gaussian go up to this distance
S_cat = (S >= 0) + 1;   % Category is 1 for L and 2 for R

% Extract 1st SIGMA for L stimulus, 2nd SIGMA for R stimulus
sigma_vec = sigma(sub2ind(size(sigma),(1:size(sigma,1))',S_cat));
X = bsxfun(@plus, S, bsxfun(@times, sigma_vec, MAXSD*linspace(-1,1,nx)));
W = normpdf(linspace(-MAXSD,MAXSD,nx));
W = W./qtrapz(W);

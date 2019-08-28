function y = bsxfun_normpdf(x,mu,sigma)
%BSXFUN_NORMPDF Vectorized normal probability density function (pdf).
%   Y = BSXFUN_NORMPDF(X,MU,SIGMA) returns the pdf of the normal 
%   distribution with mean MU and standard deviation SIGMA, evaluated at 
%   the values in X. Dimensions of X, MU, and SIGMA must either match, or 
%   be equal to one. Computation of the pdf is performed with singleton
%   expansion enabled via BSXFUN. The size of Y is the size of the input 
%   arguments (expanded to non-singleton dimensions).
%
%   All elements of SIGMA are assumed to be non-negative (no checks).
%
%   See also BSXFUN, BSXFUN_NORMCDF, NORMPDF.

%   Author: Luigi Acerbi
%   Release date: 15/07/2015

if nargin<3
    error('bmp:bsxfun_normpdf:TooFewInputs','Input argument X, MU or SIGMA are undefined.');
end

try
    if isscalar(mu)
        y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, x - mu, sigma).^2), sigma)/sqrt(2*pi);
    elseif isscalar(sigma)
        y = exp(-0.5*(bsxfun(@minus, x, mu)/sigma).^2)/(sigma*sqrt(2*pi));
    else
        y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma).^2), sigma)/sqrt(2*pi);
    end
catch
    error('bmp:bsxfun_normpdf:InputSizeMismatch',...
          'Non-singleton dimensions must match in size.');
end

function y = mtruncgausspdf(x,a,b,c,d)
%MGAUSSPDF Multivariate truncated Gaussian probability density function.

% a < b < c < d

y = zeros(size(x));

mu = 0.5*(b+c);
sigma = 0.5*(c-b);

% Normalization factor
lnf = -log(normcdf(d,mu,sigma) - normcdf(a,mu,sigma));

for ii = 1:size(x,2)
    y(:,ii) = -0.5*((x(:,ii)-mu(ii))./sigma(ii)).^2 - 0.5*log(2*pi*sigma(ii)^2) + lnf(ii);
end

y = exp(sum(y,2));
function y = mtrunctpdf(x,a,b,c,d,df)
%MTRUNCTPDF Multi-univariate truncated Student's t probability density function.

if nargin < 6 || isempty(df); df = 7; end

% a < b < c < d

y = zeros(size(x));

mu = 0.5*(b+c);
sigma = 0.5*(c-b);

% Normalization factor
lnf = gammaln((df+1)/2) - gammaln(df/2) - 0.5*log(df*pi) - log(sigma) ...
    -log(tcdf((d-mu)./sigma,df) - tcdf((a-mu)./sigma,df));

for ii = 1:size(x,2)
    y(:,ii) = -(df+1)/2*log1p( ((x(:,ii)-mu(ii))./sigma(ii)).^2/df ) + lnf(ii);
end

y = exp(sum(y,2));
function y = mquadpdf(x,a,b,c,d)
%MQUADPDF Multivariate truncated quadratic probability density function.

% a < b < c < d

y = zeros(size(x));

for ii = 1:size(x,2)
    nf = c - b + (b - a + d - c)*2/3;
    
    idx = x(:,ii) >= a(ii) & x(:,ii) < b(ii);
    y(idx,ii) = (1 - ((x(idx,ii) - b(ii)) / (b(ii) - a(ii))).^2) / nf(ii);
    
    idx = x(:,ii) >= b(ii) & x(:,ii) < c(ii);
    y(idx,ii) = 1 / nf(ii);
    
    idx = x(:,ii) >= c(ii) & x(:,ii) < d(ii);
    y(idx,ii) = (1 - ((x(idx,ii) - c(ii))/(d(ii) - c(ii))).^2) / nf(ii);    
end

y = max(prod(y,2),0);








function rldist = runlengthdist(tau,p)

N = 1e6;

% Min/max block length
rlmin = 20;
rlmax = 100;

maxn = 500;
counts = zeros(1,maxn);

for iter = 1:50    
    idx = true(N,1);
    x = zeros(N,1);
    
    while any(idx)
        Nactive = sum(idx);
        x(idx) = x(idx) + min(rlmax,max(rlmin,round(exprnd(tau,[Nactive,1]))));
        idx(idx) = rand(Nactive,1) < p;    
    end
    
    n = histcounts(x,'BinLimits',[0.5,maxn+0.5],'BinMethod','integer');    
    counts = counts + n;
end

rldist = counts / sum(counts);


end
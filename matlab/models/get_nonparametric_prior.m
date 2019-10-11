function [pvec,xx] = get_nonparametric_prior(lnp_hyp,xx)
%GET_NONPARAMETRIC_PRIOR Calculate hyperprior over probabilities.

MIN_P = 1e-6;

if nargin < 2 || isempty(xx)
    np = 2001;
    xx = linspace(MIN_P,1-MIN_P,np);
end

[Ns,nhyp] = size(lnp_hyp);
p = [lnp_hyp(:,1:nhyp/2),zeros(Ns,1),lnp_hyp(:,nhyp/2+1:end)];
px = min(max(MIN_P,linspace(0,1,size(p,2))),1-MIN_P);    
lnpvec = interp1(px,p',xx,'pchip')';

lnZ = max(lnpvec,[],2);
dx = xx(2) - xx(1);
pvec = exp(bsxfun(@minus,lnpvec,lnZ));
pvec = bsxfun(@rdivide,pvec,sum(pvec,2)*dx);

end

function plot_hyperprior(mice_name,model_name)

if nargin < 1 || isempty(mice_name); mice_name = get_mice_list(); end
if nargin < 2 || isempty(model_name); model_name = 'exponential_hyperprobs_nakarushton'; end

fontsize = 16;

MIN_P = 1e-6;
np = 2001;
xx = linspace(MIN_P,1-MIN_P,np);
dx = xx(2) - xx(1);

pvec = zeros(numel(mice_name),np);

nrows = 5;
ncols = 5;

for iMouse = 1:numel(mice_name)
    
    subplot(nrows,ncols,iMouse);
    
    params = load_model_fit(mice_name{iMouse},model_name);
    
    lnp_hyp = params.lnp_hyp;    
    p = [log(MIN_P) lnp_hyp 0 fliplr(lnp_hyp) log(MIN_P)];
    px = min(max(MIN_P,linspace(0,1,numel(p))),1-MIN_P);    
    lnpvec = interp1(px,p,xx,'pchip');
    lnZ = max(lnpvec);
    pvec(iMouse,:) = exp(lnpvec - lnZ)/sum(exp(lnpvec - lnZ)*dx);
    
    plot(xx,pvec(iMouse,:),'k-','LineWidth',1); hold on;    
    
    name = mice_name{iMouse};
    name(name == '_') = '-';
    prettyplot(name,fontsize);
    
    if iMouse ~= 21; xlabel(''); ylabel(''); end
end

subplot(nrows,ncols,iMouse+1);
plot(xx,mean(pvec),'-k','LineWidth',3);
prettyplot('Average mouse',fontsize);

end

function prettyplot(titlestr,fontsize)

set(gca,'TickDir','out');
box off;
xlim([0 1]);
ylim([0 15]);
title(titlestr);
xlabel('P(Left)','Fontsize',fontsize);
ylabel('Probability','Fontsize',fontsize);
set(gcf,'Color','w');


end

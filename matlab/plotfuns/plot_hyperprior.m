function plot_hyperprior(mice_name,model_name,data_suffix)

if nargin < 1 || isempty(mice_name); mice_name = get_mice_list(4); end
if nargin < 2 || isempty(model_name); model_name = 'exponential_contrastnoise_hyperprobs'; end
if nargin < 3; data_suffix = []; end

fontsize = 16;

MIN_P = 1e-6;
np = 2001;
xx = linspace(MIN_P,1-MIN_P,np);
dx = xx(2) - xx(1);

pvec = zeros(numel(mice_name),np);

nrows = 7;
ncols = 8;

for iMouse = 1:numel(mice_name)
    
    subplot(nrows,ncols,iMouse);
    
    data_name = mice_name{iMouse};
    if ~isempty(data_suffix); data_name = [data_name '_' data_suffix]; end
    params = load_model_fit(data_name,model_name);
    
    nhyp = numel(params.lnp_hyp);
    lnp_hyp = params.lnp_hyp;    
    p = [lnp_hyp(1:nhyp/2),0,lnp_hyp(nhyp/2+1:end)];
    px = min(max(MIN_P,linspace(0,1,numel(p))),1-MIN_P);    
    xx = linspace(MIN_P,1-MIN_P,np);
    lnpvec = interp1(px,p,xx,'pchip');
    
    lnZ = max(lnpvec);
    pvec(iMouse,:) = exp(lnpvec - lnZ)/sum(exp(lnpvec - lnZ)*dx);
    
    plot(xx,pvec(iMouse,:),'k-','LineWidth',1); hold on;    
    
    name = mice_name{iMouse};
    name(name == '_') = '-';
    prettyplot(name,fontsize);
    
    if iMouse ~= (nrows*(ncols-1)); xlabel(''); ylabel(''); end
end

subplot(nrows,ncols,iMouse+1);
plot(xx,mean(pvec),'-k','LineWidth',3);
prettyplot('Average mouse',fontsize);

end

function prettyplot(titlestr,fontsize)

set(gca,'TickDir','out');
box off;
xlim([0 1]);
ylim([0 10]);
text(0.05,0.85,titlestr,'FontSize',fontsize*0.75,'Units','normalized');
%title(titlestr);
xlabel('P(Left)','Fontsize',fontsize);
ylabel('Probability','Fontsize',fontsize);
set(gcf,'Color','w');


end

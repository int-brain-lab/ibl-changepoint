function plot_hyperprior(mice_name,model_name,data_suffix,mcmc_flag)

if nargin < 1 || isempty(mice_name); mice_name = get_mice_list(4); end
if nargin < 2 || isempty(model_name); model_name = 'exponential_contrastnoise_hyperprobs'; end
if nargin < 3; data_suffix = []; end
if nargin < 4 || isempty(mcmc_flag); mcmc_flag = false; end

fontsize = 16;

pvec = []; xx = [];
N = numel(mice_name)+1;

nrows = floor(sqrt(N));
ncols = ceil(N/nrows);
MinTrials = 300;

for iMouse = 1:numel(mice_name)
    
    subplot(nrows,ncols,iMouse);
    
    data_name = mice_name{iMouse};
    if ~isempty(data_suffix); data_name = [data_name '_' data_suffix]; end
    params = load_model_fit(data_name,model_name);
    
    if mcmc_flag
        nhyp = size(params.lnp_hyp,2);
        if isfield(params,'mcmc_fits')
            samples = params.mcmc_fits.samples{1}(:,end-nhyp+1:end);
            [pvec1,xx] = get_nonparametric_prior(samples,xx);
            pvec1 = mean(pvec1,1);
        else
            pvec1 = NaN(1,size(pvec,2));
        end
    else
        [pvec1,xx] = get_nonparametric_prior(params.lnp_hyp,xx);
    end
    
    if isempty(pvec); pvec = zeros(numel(mice_name),size(pvec1,2)); end
    pvec(iMouse,:) = pvec1;
    
    if params.mle_fits.ndata < MinTrials
        pvec(iMouse,:) = NaN;
    end    
    
    plot(xx,pvec(iMouse,:),'k-','LineWidth',1); hold on;    
    
    name = mice_name{iMouse};
    name(name == '_') = '-';
    prettyplot(name,fontsize);
    
    if iMouse ~= (nrows*(ncols-1)); xlabel(''); ylabel(''); end
end

subplot(nrows,ncols,iMouse+1);
plot(xx,nanmean(pvec),'-k','LineWidth',3);
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

%PLOT_SFN2019_MODEL Plot results of model comparison
function plot_sfn2019_model(tab,metric)

switch metric
    case 'cvll'
        mcm = tab.cvll;
    case 'bic'
        mcm = -tab.bic;
end

mcm = mcm(all(isfinite(mcm),2),:);  % Remove incomplete fits

% mcm = tab.cvll ./ tab.ntrials;
[~,best_idx] = max(nanmean(mcm,1));
mcm = mcm(:,best_idx) - mcm;

xx = 1:4;
yy = nanmean(mcm,1);

col = [70,162,221]/255;

bar(xx,yy,0.6,'FaceColor',col,'EdgeColor','none'); hold on;
yerr = stderr(mcm);
errorbar(xx,yy,yerr,'k','LineStyle','none','LineWidth',3,'CapSize',0);

set(gcf,'Color','w');
set(gca,'TickDir','out');
box off;


set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',0:500:1e4);
set(gca,'FontSize',24);
set(gcf,'position',[489,344,800,420]);



end
function plot_runlength_posterior(data,params,titlestr)
%PLOT_TRIALS Running plot of trial probabilities.

if nargin < 3; titlestr = []; end

axesfontsize = 14;
fontsize = 18;


[~,output] = nllfun([],params,data);

trials = 1:numel(data.p_true);

post = sum(output.fullpost,3);
surf(post','edgecolor','none');
view([0 90]);
colormap(1-gray);

set(gcf,'Color','w');
set(gca,'TickDir','out','FontSize',axesfontsize);
box off;
xlabel('Trial','FontSize',fontsize);
ylabel('Run length','FontSize',fontsize);
set(gca,'XTick',50:50:max(trials),'YTick',[1 50 100]);
xlim([1,max(trials)]);
ylim([0,100]);
if ~isempty(titlestr); title(titlestr,'FontSize',fontsize); end


end
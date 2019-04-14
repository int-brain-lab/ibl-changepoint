function plot_trials(data,params,titlestr)
%PLOT_TRIALS Running plot of trial probabilities.

if nargin < 3; titlestr = []; end

fontsize = 14;

[~,output] = nllfun([],params,data);

trials = 1:numel(data.p_true);

h(1) = plot(trials,data.p_true,'Color',[0 0 0],'LineWidth',3);
hold on;
h(2) = plot(trials,output.p_estimate,'Color',[0.2 0.2 0.2],'LineWidth',1);
h(3) = plot(trials,output.resp_model,'Color',[0 0.8 0.9],'LineWidth',1);

idx = data.C == 1;
h(4) = scatter(trials(idx),double(data.resp_obs(idx)==1),'ok','MarkerFaceColor',[0.7 0 0.7],'MarkerEdgeColor',[1 1 1]);
h(5) = scatter(trials(~idx),double(data.resp_obs(~idx)==1),'ok','MarkerFaceColor',[0.5 0.5 0],'MarkerEdgeColor',[1 1 1]);

set(gcf,'Color','w');
set(gca,'TickDir','out');
box off;
xlabel('Trial','FontSize',fontsize);
ylabel('Probability Left','FontSize',fontsize);
set(gca,'XTick',50:50:max(trials),'YTick',0:0.2:1);
xlim([1,max(trials)]);
ylim([-0.05,1.05]);

hl = legend(h,'True probability','Estimated probability','Choice probability','Responses (L stim)','Responses (R stim)');
set(hl,'Box','off','FontSize',fontsize,'Location','NorthEastOutside');

end
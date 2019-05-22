function plot_trials(data,params,titlestr,plotchoice_flag,plotbycontrast_flag,fontsizes)
%PLOT_TRIALS Running plot of trial probabilities.

if nargin < 3; titlestr = []; end
if nargin < 4 || isempty(plotchoice_flag); plotchoice_flag = true; end
if nargin < 5 || isempty(plotbycontrast_flag); plotbycontrast_flag = false; end
if nargin < 6 || isempty(fontsizes); fontsizes = [18 14]; end

fontsize = fontsizes(1);
axesfontsize = fontsizes(2);

[~,output] = nllfun([],params,data);

trials = 1:numel(data.p_true);

h(1) = plot(trials,data.p_true,'Color',[0 0 0],'LineWidth',3);
ltext{1} = 'True probability';

hold on;
h(end+1) = plot(trials,output.p_estimate,'Color',[83 154 206]/255,'LineWidth',3);
ltext{end+1} = 'Estimated probability';

if plotchoice_flag
    h(end+1) = plot(trials,output.resp_model,'Color',[0 0.8 0.9],'LineWidth',1);
    ltext{end+1} = 'Choice probability';
end

if plotbycontrast_flag
    idx = data.C == 1;
    
    h(end+1) = scatter(trials(idx),ones(size(trials(idx))),'ok','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
    ltext{end+1} = 'Stimuli';
    scatter(trials(~idx),zeros(size(trials(~idx))),'ok','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
    % ltext{end+1} = 'Right stimuli';
    
%     ccols = linspace(0.8,0,numel(data.contrasts_vec));    
%     for iContrast = 1:numel(data.contrasts_vec)-1        
%         idx1 = data.C == 1 & data.contrasts == data.contrasts_vec(iContrast);
%         idx0 = data.C ~= 1 & data.contrasts == data.contrasts_vec(iContrast);
%         scatter(trials(idx1),ones(size(trials(idx1))),'ok','MarkerFaceColor',ccols(iContrast)*[1 1 1],'MarkerEdgeColor',[1 1 1]);
%         scatter(trials(idx0),zeros(size(trials(idx0))),'ok','MarkerFaceColor',ccols(iContrast)*[1 1 1],'MarkerEdgeColor',[1 1 1]);
%     end
    
    
else
    idx = data.C == 1;
    h(end+1) = scatter(trials(idx),double(data.resp_obs(idx)==1),'ok','MarkerFaceColor',[0.7 0 0.7],'MarkerEdgeColor',[1 1 1]);
    ltext{end+1} = 'Responses (L stim)';
    h(end+1) = scatter(trials(~idx),double(data.resp_obs(~idx)==1),'ok','MarkerFaceColor',[0.5 0.5 0],'MarkerEdgeColor',[1 1 1]);
    ltext{end+1} = 'Responses (R stim)';
end

set(gcf,'Color','w');
set(gca,'TickDir','out','FontSize',axesfontsize);
box off;
xlabel('Trial','FontSize',fontsize);
ylabel('P(Left)','FontSize',fontsize);
set(gca,'XTick',50:50:max(trials),'YTick',0:0.2:1);
xlim([1,max(trials)]);
ylim([-0.05,1.05]);
if ~isempty(titlestr); title(titlestr,'FontSize',fontsize); end

hl = legend(h,ltext{:});
set(hl,'Box','off','FontSize',fontsize,'Location','NorthEastOutside');

end
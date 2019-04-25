function plot_changepoints(data,params,N,titlestr,fitinfo_flag)
%PLOT_CHANGEPOINTS Plot performance around changepoints.

if nargin < 2; params = []; end
if nargin < 3 || isempty(N); N = 20; end
if nargin < 4; titlestr = []; end
if nargin < 5 || isempty(fitinfo_flag); fitinfo_flag = true; end

fontsize = 14;


% List of all change points
cps = find(abs(abs(diff(data.p_true)) - 0.6) < eps)+1;

Ncontrasts = numel(data.contrasts_vec);

% Create table for each changepoint * surrounding trials * contrast * match
tab_match = zeros(size(cps,1),2*N,Ncontrasts,2);
tab_counts = tab_match;

p_true = data.p_true;
plot_data_flag = true;

for iC = 1:size(cps,1)
    offset = -N+1:N;
    idx_list = cps(iC)+offset(1):cps(iC)+offset(end);
    
    % Check that we are inside the trials
    ok_trials = idx_list > 1 & idx_list < size(data.tab,1);
    offset = offset(ok_trials); idx_list = idx_list(ok_trials);
    
    % Check that we stay within the same session as the change point
    ok_trials = data.tab(idx_list,2) == data.tab(cps(iC),2);
    offset = offset(ok_trials); idx_list = idx_list(ok_trials);
    
    % Loop over (valid) surrounding trials
    for iTrial = 1:numel(idx_list)
        idx = idx_list(iTrial);
        
        p_block = p_true(cps(iC));
        
        tab_index = offset(iTrial) + N;
        
        if (p_block > 0.5 && data.S(idx) < 0) || (p_block < 0.5 && data.S(idx) > 0) ...
                || (data.contrasts(idx) == 0)
            match_column = 1;
        else
            match_column = 2;
        end
                
        tab_counts(iC,tab_index,data.contrasts_idx(idx),match_column) = tab_counts(tab_index) + 1;
        % Response matches block type
        if (data.resp_obs(idx) > 0 && p_block < 0.5) || (data.resp_obs(idx) < 0 && p_block > 0.5)
            tab_match(iC,tab_index,data.contrasts_idx(idx),match_column) = tab_match(tab_index) + 1;          
        end        
    end    
    
end

contrast_color = linspace(0,0.8,size(tab_match,3))'*[1 1 1];



for iContrast = 1:size(tab_match,3)
    match = sum(tab_match(:,:,iContrast,1),1);
    total = sum(tab_counts(:,:,iContrast,1),1);
    %match = sum(sum(tab_match(:,:,iContrast,:),1),4);
    %total = sum(sum(tab_counts(:,:,iContrast,:),1),4);
    offset = -N+1:N;    
    idx = total > 0;
    
    match = match(idx); total = total(idx); offset = offset(idx);
    
    p = match ./ total;
    s = sqrt(p.*(1-p))./sqrt(total);
    
    col = contrast_color(iContrast,:);
    
%    errorbar(offset,p,s,...
%        'LineStyle','none','LineWidth',2,'Color',col,'capsize',0); hold on;
    h(iContrast) = plot(offset,p,...
        'LineStyle','none','LineWidth',1,'Marker','o','Color',col,'MarkerFaceColor',col,'MarkerEdgeColor',[1 1 1]); hold on;
    
    legtext{iContrast} = ['c = ' num2str(data.contrasts_vec(iContrast)*100,'%.3g') '%'];
end

% %% Decorate plot
xlim([-N+0.5,N+0.5]);
ylim([0 1]);
set(gca,'TickDir','out');
box off;
set(gcf,'Color','w');
set(gca,'Xtick',-N+1:N);
set(gca,'Ytick',0:0.5:1);
% set(gca,'XTickLabel',xtext);
xlabel('Trial with respect to change point','FontSize',fontsize);
ylabel('P(match change-point block)','FontSize',fontsize);
if ~isempty(titlestr); title(titlestr,'FontSize',fontsize); end
 
if plot_data_flag
     hl = legend(h,legtext{:});
     set(hl,'Location','SouthEast','Box','off','FontSize',fontsize);
end


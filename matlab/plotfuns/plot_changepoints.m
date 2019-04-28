function plot_changepoints(data,params,N,titlestr,fitinfo_flag)
%PLOT_CHANGEPOINTS Plot performance around changepoints.

if nargin < 2; params = []; end
if nargin < 3 || isempty(N); N = 20; end
if nargin < 4; titlestr = []; end
if nargin < 5 || isempty(fitinfo_flag); fitinfo_flag = true; end

fontsize = 14;

if ischar(data); data = read_data_from_csv(data); end

if ~isempty(params)
    [~,output] = nllfun([],params,data);

    % Compute model average if passing multiple samples
    Nparams = numel(output);
    resp_model = zeros(size(output(1).resp_model,1),Nparams);
    for iParam = 1:Nparams
        resp_model(:,iParam) = output(iParam).resp_model;
    end    
end

% List of all change points
cps = find(abs(abs(diff(data.p_true)) - 0.6) < eps)+1;

Ncontrasts = numel(data.contrasts_vec);

% Create table for each changepoint * surrounding trials * contrast * match
tab_match = zeros(size(cps,1),2*N,Ncontrasts,2);
tab_match_model = tab_match;
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
        
        if ~isempty(params)
            if p_block > 0.5
                tab_match_model(iC,tab_index,data.contrasts_idx(idx),match_column) = tab_match_model(tab_index) + mean(resp_model(idx,:),2);
            else
                tab_match_model(iC,tab_index,data.contrasts_idx(idx),match_column) = tab_match_model(tab_index) + 1 - mean(resp_model(idx,:),2);
            end
        end
        
    end    
    
end


contrast_groups = {[1 2 3], [4 5]};
contrast_color{1} = linspace(0,0.6,numel(contrast_groups))'*[0.5 0.5 1];
contrast_color{2} = linspace(0.4,0.8,numel(contrast_groups))'*[1 0.5 0.5];

legtext = {'low contrasts', 'high contrasts'};

plot([0 0],[0,1],'k--','LineWidth',1); hold on;

for iContrast = 1:numel(contrast_groups)
    cc = contrast_groups{iContrast};
    
    for j = 1:2
        match = sum(sum(tab_match(:,:,cc,j),1),3);
        total = sum(sum(tab_counts(:,:,cc,j),1),3);    
        match_model = sum(sum(tab_match_model(:,:,cc,j),1),3);

        %match = sum(sum(tab_match(:,:,iContrast,:),1),4);
        %total = sum(sum(tab_counts(:,:,iContrast,:),1),4);
        offset = -N+1:N;
        idx = total > 0;

        match = match(idx); total = total(idx); offset = offset(idx);

        p = match ./ total;
        s = sqrt(p.*(1-p))./sqrt(total);
        p_model = match_model ./ total;

        col = contrast_color{j}(iContrast,:);

        legtext{iContrast} = [legtext{iContrast} ' ('];
        for i = 1:numel(cc)
            legtext{iContrast} = [legtext{iContrast} num2str(data.contrasts_vec(cc(i))*100,'%.3g') '%'];    
            if i < numel(cc); legtext{iContrast} = [legtext{iContrast} ', ']; end
        end
        legtext{iContrast} = [legtext{iContrast} ')'];

        if ~isempty(params)        
    %        errorbar((1:numel(cc_vec))+offset,model_mean,model_stderr,...
    %            'LineStyle','none','LineWidth',2,'Color',col,'capsize',0); hold on;
            h(iContrast+numel(contrast_groups)) = plot(offset,p_model,...
                'LineStyle','-','LineWidth',2,'Color',col);
            legtext{iContrast+numel(contrast_groups)} = legtext{iContrast};
        end

        errorbar(offset,p,s,...
            'LineStyle','none','LineWidth',2,'Color',col,'capsize',0); hold on;
        h(iContrast) = plot(offset,p,...
            'LineStyle','none','LineWidth',1,'Marker','o','Color',col,'MarkerFaceColor',col,'MarkerEdgeColor',[1 1 1]); hold on;
    
    end    
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
     set(hl,'Location','Best','Box','off','FontSize',fontsize);
end


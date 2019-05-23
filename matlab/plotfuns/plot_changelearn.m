function plot_changelearn()

fontsize = 18;
axesfontsize = 14;

mice_list = get_mice_list();
model_name = 'changelearn_nakarushton';

tmax = 1e4;

ylims = [1 100; 1 100; 0 1; 0 1];
param_names = {'Run-length exponential scale \tau', 'Minimum block length','P(Left) in Right blocks','P(Left) in Left blocks'};
true_vals = [60,20,0.2,0.8];

t = unique(round(exp(linspace(log(1),log(tmax),1e3))));

mu_all = nan(numel(mice_list),numel(t),4);

for i = 1:numel(mice_list)
    mouse_name = mice_list{i};    
    data = read_data_from_csv(mouse_name);
    params = load_model_fit(mouse_name,model_name);
    
    if isempty(params); continue; end
    
    p_grid = params.output.p_grid;
    NumTrials = size(p_grid,1);

    % Compute posterior over change-point parameters
    loglike = zeros(size(p_grid));
    loglike(data.C(1:NumTrials) == 1,:) = log(p_grid(data.C(1:NumTrials) == 1,:));
    loglike(data.C(1:NumTrials) ~= 1,:) = log(1-p_grid(data.C(1:NumTrials) ~= 1,:));

    % Flat prior for the moment
    logpost_grid = cumsum(loglike,1);    
    post_grid = exp(bsxfun(@minus,logpost_grid,max(logpost_grid,[],2)));
    post_grid = bsxfun(@rdivide,post_grid,sum(post_grid,2));
    
    % Only keep subset of trials    
    post_grid = post_grid(t(t<=NumTrials),:);
    
    tt = t(t<=NumTrials);

    p3(1,:,:) = params.output.param_grid;
    mu = squeeze(sum(post_grid.*p3,2));
        
    for iParam = 1:4
        subplot(2,2,iParam);
        
        
        if iParam < 3; mu(:,iParam) = exp(mu(:,iParam)); end
        
        mu_all(i,1:numel(tt),iParam) = mu(:,iParam);        
        
        % s2 = max(0,squeeze(sum(post_grid.*p3.^2,2)) - mu.^2);

        if i == 1
            h(1) = plot([1,tmax],true_vals(iParam)*[1 1],'-b','LineWidth',3);
            hold on;
        end
        
        idx = iParam; 
        h(2) = plot(tt,mu(:,idx),'-','Color',0.5*[1 1 1],'LineWidth',1);
        hold on;
        %plot(tt,mu(:,idx) + sqrt(s2(:,idx)),':k');
        %plot(tt,mu(:,idx) - sqrt(s2(:,idx)),':k');
        box off;
        set(gca,'Xscale','log','TickDir','out');
        set(gcf,'Color','w');
        set(gca,'XTick',[1,10,100,1e3,1e4],...
            'XTickLabel',{'1','10','100','1000','10000'},...
            'Fontsize',axesfontsize);
        xlim([1,tmax]);
        ylim(ylims(iParam,:));
        title(param_names{iParam},'Fontsize',fontsize);
        xlabel('Trial','Fontsize',fontsize);
        ylabel('Posterior mean','Fontsize',fontsize);
        
        if iParam < 3
            set(gca,'YTick',0:20:100);
        else
            set(gca,'YTick',[0,0.2,0.5,0.8,1]);            
        end
        
    end
    
    drawnow
end

for iParam = 1:4
    subplot(2,2,iParam);
    h(3) = plot(t,nanmean(mu_all(:,:,iParam),1),'-k','LineWidth',3);
    if iParam == 4
        hl = legend(h,{'True value','Single datasets','Grand average'});
        set(hl,'Box','off','Fontsize',axesfontsize,'Location','best');
    end
end
function output = plot_performance(data,params,titlestr,fitinfo_flag)
%PLOT_PERFORMANCE Plot histogram of performance.

if nargin < 2; params = []; end
if nargin < 3; titlestr = []; end
if nargin < 4 || isempty(fitinfo_flag); fitinfo_flag = true; end

fontsize = 14;

if ~isempty(params)
    [~,output] = nllfun([],params,data);

    % Compute model average if passing multiple samples
    Nparams = numel(output);
    resp_model = zeros(size(output(1).resp_model,1),Nparams);
    for iParam = 1:Nparams
        resp_model(:,iParam) = output(iParam).resp_model;
    end
    
    model_pcorrect = resp_model.*(data.C == 1) + (1 - resp_model).*(data.C == 2);
end

plot_data_flag = any(isfinite(data.resp_obs));

cc_vec = unique(data.contrasts);
cc_vec = [cc_vec(cc_vec>0)',Inf];

for i_contrast = 1:numel(cc_vec)
    if isfinite(cc_vec(i_contrast))
        idx = data.contrasts == cc_vec(i_contrast);  
    else
        idx = data.contrasts > 0;
    end
    pp = mean(data.resp_correct(idx));
    data_mean(i_contrast) = pp;
    data_stderr(i_contrast) = sqrt(pp*(1-pp)/sum(idx));
    
    if ~isempty(params)
        
        model_mean(i_contrast) = mean(mean(model_pcorrect(idx,:),2));
        model_stderr(i_contrast) = stderr(mean(model_pcorrect(idx,:),1));
    end
    
    if isfinite(cc_vec(i_contrast))
        xtext{i_contrast} = [num2str(cc_vec(i_contrast)*100,'%.3g') '%'];
    else
        xtext{i_contrast} = 'all > 0%';        
    end
end

if plot_data_flag
    col = [0 0 0];
    errorbar(1:numel(cc_vec),data_mean,data_stderr,...
        'LineStyle','none','LineWidth',2,'Color',col,'capsize',0); hold on;
    h(1) = plot(1:numel(cc_vec),data_mean,...
        'LineStyle','none','LineWidth',1,'Marker','o','Color',col,'MarkerFaceColor',col,'MarkerEdgeColor',[1 1 1]);
    legtext{1} = ['data'];
end
 
if ~isempty(params)
    col = [0 0 1];
    offset = 0.1;
    errorbar((1:numel(cc_vec))+offset,model_mean,model_stderr,...
        'LineStyle','none','LineWidth',2,'Color',col,'capsize',0); hold on;
    h(2) = plot((1:numel(cc_vec))+offset,model_mean,...
        'LineStyle','none','LineWidth',1,'Marker','o','Color',col,'MarkerFaceColor',col,'MarkerEdgeColor',[1 1 1]);
    legtext{2} = params.model_desc;
end

output.data_mean = data_mean;
output.model_mean = model_mean;

% %% Loop over block types
% for i_prob = 1:numel(p_true_vec)
%     idx_trial = data.p_true == p_true_vec(i_prob);    
%     resp = double(data.resp_obs(idx_trial) == 1);
%     model_resp = 1 - resp_model(idx_trial);
%         
%     cc_vec = unique(cc);
%     
%     mean_resp = zeros(1,numel(cc_vec));
%     std_resp = zeros(1,numel(cc_vec));
%     mean_model = zeros(1,numel(cc_vec));
%     
%     for i_contrast = 1:numel(cc_vec)
%         idx_contrast = cc == cc_vec(i_contrast);        
%         rr = resp(idx_contrast);
%         mean_resp(i_contrast) = mean(rr);
%         std_resp(i_contrast) = sqrt(mean_resp(i_contrast)*(1-mean_resp(i_contrast))/numel(rr));
%         
%         mean_model(i_contrast) = mean(model_resp(idx_contrast));
%     end
%     
%     if numel(p_true_vec) > 1
%         jitter = ((i_prob-1)/(numel(p_true_vec)-1)-0.5)*0.01;
%     else
%         jitter = 0;
%     end
%     
%     if p_true_vec(i_prob) == 0.5
%         col = [0 0 0];
%     elseif p_true_vec(i_prob) > 0.5
%         col = [0.7 0 0.7];
%     else
%         col = [0.5 0.5 0];        
%     end
%     
%     if ~plot_data_flag; col = 1 - 0.15*(1-col); end
%     
%     if strcmp(params.model_name,'psychofun')
%         % Plot continuous psychometric function
%         Nc = 200;
%         params1.psycho_mu = params.psycho_mu(i_prob);
%         params1.psycho_sigma = params.psycho_sigma(i_prob);
%         params1.psycho_gammalo = params.psycho_gammalo(i_prob);
%         params1.psycho_gammahi = params.psycho_gammahi(i_prob);        
%         cc_psy = linspace(-1,1,Nc);
%         data1.signed_contrasts = cc_psy';
%         data1.p_true = p_true_vec(i_prob)*ones(Nc,1);
%         data1.resp_obs = -1*ones(Nc,1);        
%         [~,output1] = psycho_nll(params1,data1);
%         mean_psy = 1 - output1.resp_model;
%         plot(cc_psy,mean_psy,'LineStyle','-','Color',col,'LineWidth',2);        
%     else    
%         plot(cc_vec,mean_model,'LineStyle','-','Color',col,'LineWidth',2);
%     end
% 
%     hold on;
%     if plot_data_flag
%         errorbar(cc_vec+jitter,mean_resp,std_resp,...
%             'LineStyle','none','LineWidth',2,'Color',col,'capsize',0);
%         h(i_prob) = plot(cc_vec+jitter,mean_resp,...
%             'LineStyle','none','LineWidth',1,'Marker','o','Color',col,'MarkerFaceColor',col,'MarkerEdgeColor',[1 1 1]);
%         legtext{i_prob} = ['probabilityLeft ' num2str(p_true_vec(i_prob))];
%     end
% end

% %% Decorate plot
xlim([0.5 numel(cc_vec)+0.5]);
ylim([0.5 1]);
set(gca,'TickDir','out');
box off;
set(gcf,'Color','w');
set(gca,'Xtick',1:numel(cc_vec));
set(gca,'Ytick',0.5:0.1:1);
set(gca,'XTickLabel',xtext);
xlabel('Contrast level','FontSize',fontsize);
ylabel('P(correct)','FontSize',fontsize);
if ~isempty(titlestr); title(titlestr,'FontSize',fontsize); end
 
if plot_data_flag
     hl = legend(h,legtext{:});
     set(hl,'Location','SouthEast','Box','off','FontSize',fontsize);
end
% 
% if fitinfo_flag
%     aic = 2*nLL + 2*numel(params.theta);
%     text(0.1,0.8,['AIC = ' num2str(aic,'%.2f')],'FontSize',fontsize,'Units','Normalized');
%     text(0.1,0.875,['LL = ' num2str(-nLL,'%.2f')],'FontSize',fontsize,'Units','Normalized');
%     text(0.1,0.95,['# params = ' num2str(numel(params.theta))],'FontSize',fontsize,'Units','Normalized');
% 
%     % Plot some parameters
%     if strncmp(params.model_name,'changepoint',numel('changepoint'))
%         text(0.1,0.5,['tau = ' num2str(params.runlength_tau,'%.2f')],'FontSize',fontsize,'Units','Normalized');
%         text(0.1,0.425,['min-rl = ' num2str(params.runlength_min,'%.2f')],'FontSize',fontsize,'Units','Normalized');
%         text(0.1,0.35,['p-low = ' num2str(min(params.p_vec),'%.2f')],'FontSize',fontsize,'Units','Normalized');
%     end
% end

end
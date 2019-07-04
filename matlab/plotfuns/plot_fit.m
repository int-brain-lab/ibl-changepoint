function output = plot_fit(data,params,titlestr,fitinfo_flag)
%PLOT_FIT Plot psychometric curves of responding RIGHT vs. contrast level.

if nargin < 3; titlestr = []; end
if nargin < 4 || isempty(fitinfo_flag); fitinfo_flag = true; end

fontsize = 14;

[nLL,output] = nllfun([],params,data);
nLL = nLL(1);
resp_model = output(1).resp_model;

% Plot psychometric curve
p_true_vec = unique(data.p_true);

plot_data_flag = any(isfinite(data.resp_obs));

%% Loop over block types
for i_prob = 1:numel(p_true_vec)
    idx_trial = data.p_true == p_true_vec(i_prob);    
    resp = double(data.resp_obs(idx_trial) == 1);
    cc = data.signed_contrasts(idx_trial);    
    model_resp = 1 - resp_model(idx_trial);
        
    cc_vec = unique(cc);
    
    mean_resp = zeros(1,numel(cc_vec));
    std_resp = zeros(1,numel(cc_vec));
    mean_model = zeros(1,numel(cc_vec));
    
    for i_contrast = 1:numel(cc_vec)
        idx_contrast = cc == cc_vec(i_contrast);        
        rr = resp(idx_contrast);
        mean_resp(i_contrast) = mean(rr);
        std_resp(i_contrast) = sqrt(mean_resp(i_contrast)*(1-mean_resp(i_contrast))/numel(rr));
        
        mean_model(i_contrast) = mean(model_resp(idx_contrast));
    end
    
    if numel(p_true_vec) > 1
        jitter = ((i_prob-1)/(numel(p_true_vec)-1)-0.5)*0.01;
    else
        jitter = 0;
    end
    
    if p_true_vec(i_prob) == 0.5
        col = [0 0 0];
    elseif p_true_vec(i_prob) > 0.5
        col = [0.7 0 0.7];
    else
        col = [0.5 0.5 0];        
    end
    
    %if ~plot_data_flag; col = 1 - 0.15*(1-col); end
    
    if strcmp(params.model_name,'psychofun')
        % Plot continuous psychometric function
        
        Nsamples = 100;
        theta = get_posterior_samples(params,Nsamples);
        if isempty(theta); theta = params.theta; end
        
        [mean_psy,std_psy,cc_psy] = psychofun(theta,params,i_prob,p_true_vec);
        
        if any(std_psy > 0)
            fill([cc_psy,fliplr(cc_psy)],[mean_psy+std_psy,fliplr(mean_psy-std_psy)],col,'FaceAlpha',0.2,'EdgeColor','none'); hold on;
        end
        
        plot(cc_psy,mean_psy,'LineStyle','-','Color',col,'LineWidth',2);
        
    else    
        plot(cc_vec,mean_model,'LineStyle','-','Color',col,'LineWidth',2);
    end

    hold on;
    if plot_data_flag
        errorbar(cc_vec+jitter,mean_resp,std_resp,...
            'LineStyle','none','LineWidth',2,'Color',col,'capsize',0);
        h(i_prob) = plot(cc_vec+jitter,mean_resp,...
            'LineStyle','none','LineWidth',1,'Marker','o','Color',col,'MarkerFaceColor',col,'MarkerEdgeColor',[1 1 1]);
        legtext{i_prob} = ['probabilityLeft ' num2str(p_true_vec(i_prob))];
    end
end

%% Decorate plot
xlim([-1.025 1.025]);
ylim([-0.1 1.1]);
set(gca,'TickDir','out');
box off;
set(gcf,'Color','w');
set(gca,'Xtick',-1:0.5:1);
set(gca,'Ytick',[0:0.2:1]);
xlabel('Signed contrast','FontSize',fontsize);
ylabel('Fraction choice Right','FontSize',fontsize);
if ~isempty(titlestr); title(titlestr,'FontSize',fontsize); end

if plot_data_flag
    hl = legend(h,legtext{:});
    set(hl,'Location','SouthEast','Box','off','FontSize',fontsize);
end

if fitinfo_flag
    aic = 2*nLL + 2*numel(params.theta);
    text(0.1,0.8,['AIC = ' num2str(aic,'%.2f')],'FontSize',fontsize,'Units','Normalized');
    text(0.1,0.875,['LL = ' num2str(-nLL,'%.2f')],'FontSize',fontsize,'Units','Normalized');
    text(0.1,0.95,['# params = ' num2str(numel(params.theta))],'FontSize',fontsize,'Units','Normalized');

    % Plot some parameters
    if strncmp(params.model_name,'changepoint',numel('changepoint'))
        text(0.1,0.5,['tau = ' num2str(params.runlength_tau,'%.2f')],'FontSize',fontsize,'Units','Normalized');
        text(0.1,0.425,['min-rl = ' num2str(params.runlength_min,'%.2f')],'FontSize',fontsize,'Units','Normalized');
        text(0.1,0.35,['p-low = ' num2str(min(params.p_vec),'%.2f')],'FontSize',fontsize,'Units','Normalized');
    end
end

end

%--------------------------------------------------------------------------
function [mean_psy,std_psy,cc_psy] = psychofun(theta,params,i_prob,p_true_vec)

Nc = 200;
cc_psy = linspace(-1,1,Nc);
data1.signed_contrasts = cc_psy';
data1.p_true = p_true_vec(i_prob)*ones(Nc,1);
data1.resp_obs = -1*ones(Nc,1);        

for i = 1:size(theta,1)

    params = setup_params(theta(i,:),params);

    % Compute continuous psychometric function
    params1.psycho_mu = params.psycho_mu(i_prob);
    params1.psycho_sigma = params.psycho_sigma(i_prob);
    params1.psycho_gammalo = params.psycho_gammalo(i_prob);
    params1.psycho_gammahi = params.psycho_gammahi(i_prob);        
    [~,output1] = psycho_nll(params1,data1);
    psy(i,:) = 1 - output1.resp_model;
    % plot(cc_psy,mean_psy,'LineStyle','-','Color',col,'LineWidth',2);        
end

mean_psy = mean(psy,1);
std_psy = std(psy,[],1);

end


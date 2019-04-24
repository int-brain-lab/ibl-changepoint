function [bias_shift,bias_prob] = plot_predicted_bias_shift(metric,mice_list,frac_flag,avg_mouse_flag)
%PLOT_PREDICTED_BIAS_SHIFT Plot results of analysis of optimal bias shift.

if nargin < 1 || isempty(metric)
    metric = 'bias_shift';
end

if nargin < 2 || isempty(mice_list)        
    % mice_list = {'CSHL_005','CSHL_007','IBL-T1','IBL-T4','ibl_witten_04','ibl_witten_05'};   % Example mice
    mice_list = {'CSHL_003','CSHL_005','CSHL_007','CSHL_008','CSHL_010','IBL-T1','IBL-T4','ZM_1084','ZM_1085','ZM_1086','ZM_1091','ZM_1092','ZM_1093','ZM_1097','ZM_1098','ibl_witten_04','ibl_witten_05','ibl_witten_06'};
end

if nargin < 3 || isempty(frac_flag)
    frac_flag = false;
end

if nargin < 4 || isempty(avg_mouse_flag)
    avg_mouse_flag = true;
end

fontsize = 14;


Nmice = numel(mice_list);
Nsamples = 1e3;

bias_shift = [];
bias_prob = [];

for iMouse = 1:Nmice
    temp = load([mice_list{iMouse} '_bias_shift_endtrain.mat'],'bias_shift','bias_prob');
    
    models = fields(temp.bias_shift)';
    for iModel = 1:numel(models)
        model_shift = temp.bias_shift.(models{iModel});
        model_pmatch = temp.bias_prob.(models{iModel}).MatchProbability;
        
        % Generate posterior over bias shift and matching probability for the data
        if strcmp(models{iModel},'data') && numel(model_shift) == 1
            temp2 = load([mice_list{iMouse} '_fits.mat'],'modelfits');
            idx = find(cellfun(@(p) strcmp(p.model_name,'psychofun'),temp2.modelfits.params),1);
            theta_rnd = get_posterior_samples(temp2.modelfits.params{idx},Nsamples);
            if ~isempty(theta_rnd)
                sampled_shift = theta_rnd(:,3) - theta_rnd(:,1);            
                model_shift = [model_shift,sampled_shift'];
                
                % Psychometric curves at zero contrast for biased blocks
                sampled_pmatch = zeros(1,Nsamples);
                
                % Psychometric curve at zero contrast for a given block
                psycho0 = @(theta,offset) psychofun(0,[theta(1+offset),exp(theta(4+offset)),theta(7+offset),theta(10+offset)]);
                
                for iSample = 1:Nsamples
                    pRblock = psycho0(theta_rnd(iSample,:),0);
                    pLblock = psycho0(theta_rnd(iSample,:),2);                    
                    sampled_pmatch(iSample) = 0.5*(pRblock + 1 - pLblock);
                end
                model_pmatch = [model_pmatch,sampled_pmatch];
            end
        end
        
        if strcmp(models{iModel},'data')
            data_shift = model_shift;
            data_bias_prob = model_pmatch;
        end
        
        % Compute bunch of summary statistics
        bias_shift = metric_stats(bias_shift,model_shift,data_shift,models{iModel},iMouse,Nsamples);
        bias_prob = metric_stats(bias_prob,model_pmatch,data_bias_prob,models{iModel},iMouse,Nsamples);
                
    end
end

% Make plot

colors = repmat(linspace(0,0.5,Nmice)',[1,3]);

color_palette = {'003f5c','2f4b7c','665191','a05195','d45087','f95d6a','ff7c43','ffa600'};

for iColor = 1:numel(color_palette)
    colors(iColor,1) = hex2dec(color_palette{iColor}(1:2))/255;
    colors(iColor,2) = hex2dec(color_palette{iColor}(3:4))/255;
    colors(iColor,3) = hex2dec(color_palette{iColor}(5:6))/255;
end

colors = [colors; 0.5*ones(20,3)];

if frac_flag
    model_offset = 0;
else
    model_offset = 0;
end

if avg_mouse_flag
    for iModel = 1:numel(models)
        bias_shift = metric_add_avg_mouse(bias_shift,models{iModel},Nmice);
        bias_prob = metric_add_avg_mouse(bias_prob,models{iModel},Nmice);
    end
    Nmice = Nmice + 1;
    mice_list{end+1} = 'average';
    colors(Nmice,:) = [0 0 0];
end

yy_max = 0;

for iModel = 1+model_offset:numel(models)
    switch metric
        case 'bias_shift'
            tab = bias_shift.(models{iModel});
            metricname = 'bias shift';
            ybottom = 0;
        case 'pmatch'
            tab = bias_prob.(models{iModel});
            metricname = 'P(match) on 0% trials';
            ybottom = 0.3;
    end
    
    dx = 0.5/(Nmice-1);
    offset = iModel - model_offset - dx*(Nmice+0.5)/2;
    
    for iMouse = 1:Nmice
        if frac_flag
            y_mle = tab.frac_mle(iMouse);
            y_post = tab.frac_mean(iMouse)+tab.frac_std(iMouse)*[-1,1];
        else
            y_mle = tab.mle(iMouse);            
            y_post = tab.mean(iMouse)+tab.std(iMouse)*[-1,1];
        end
        
        x = offset + iMouse*dx;
        plot([x x],y_post,'Color',colors(iMouse,:));
        hold on;
        if avg_mouse_flag && iMouse == Nmice
            marker = 's';
        else
            marker = 'o';
        end
        
        h(iMouse) = plot(x,y_mle,['k' marker],'MarkerFaceColor',colors(iMouse,:));
        
        yy_max = max([yy_max,y_mle,y_post]);
    end
    
    name = models{iModel};
    name(name == '_') = '-';    
    modeltick{iModel-model_offset} = name;
end

plot([0,numel(models)+1-model_offset],[1 1],'k--','LineWidth',1);
if strcmp(metric,'pmatch')
    plot([0,numel(models)+1-model_offset],0.5*[1 1],'k--','LineWidth',1);
end

xlim([0,numel(models)+1-model_offset]);
if frac_flag
    ylim([0 1.2]);
    set(gca,'YTick',[0 0.5 1]);
    ystring = ['Fraction of ' metricname];
else
    ytop = ceil(yy_max*10)/10;
    ylim([ybottom ytop]);
    if ytop < 0.5; set(gca,'YTick',ybottom:0.05:1); else; set(gca,'YTick',ybottom:0.1:1); end
    ystring = ['Optimal ' metricname];
end

box off;
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'XTick',1:numel(models)-model_offset);
set(gca,'XTickLabel',modeltick);
set(gca,'FontSize',fontsize-2);

xlabel(['Model used to predict optimal ' metricname],'FontSize',fontsize);
ylabel(ystring,'FontSize',fontsize);

for iMouse = 1:Nmice
    mice_names{iMouse} = mice_list{iMouse};
    mice_names{iMouse}(mice_names{iMouse} == '_') = '-';
end
hl = legend(h,mice_names{:});
set(hl,'Box','off','FontSize',fontsize,'Location','SouthEast');


end

%--------------------------------------------------------------------------
function metric_struct = metric_stats(metric_struct,metric_samples,data_samples,model,iMouse,Nsamples)
%METRIC_STATS Compute MLE and other summary statistics of a given metric.

metric_struct.(model).mle(iMouse) = metric_samples(1);
metric_struct.(model).frac_mle(iMouse) = ...
    metric_struct.data.mle(iMouse)/metric_struct.(model).mle(iMouse);

if numel(metric_samples) > 1
    metric_struct.(model).mean(iMouse) = mean(metric_samples(2:end));
    metric_struct.(model).std(iMouse) = std(metric_samples(2:end));

    rdata = randi(numel(data_samples)-1,[Nsamples,1])+1;
    rmodel = randi(numel(metric_samples)-1,[Nsamples,1])+1;           
    frac = data_samples(rdata)./metric_samples(rmodel);

    metric_struct.(model).frac_mean(iMouse) = mean(frac);
    metric_struct.(model).frac_std(iMouse) = std(frac);                        
else
    metric_struct.(model).mean(iMouse) = NaN;
    metric_struct.(model).std(iMouse) = NaN;            
    metric_struct.(model).frac_mean(iMouse) = NaN;
    metric_struct.(model).frac_std(iMouse) = NaN;                        
end
end

%--------------------------------------------------------------------------
function metric_struct = metric_add_avg_mouse(metric_struct,model,Nmice)
%METRIC_ADD_AVG_MOUSE Add average mouse to metric struct.

metric_struct.(model).mle(Nmice+1) = mean(metric_struct.(model).mle(1:Nmice));
metric_struct.(model).frac_mle(Nmice+1) = mean(metric_struct.(model).frac_mle(1:Nmice));
metric_struct.(model).mean(Nmice+1) = metric_struct.(model).mle(Nmice+1);
metric_struct.(model).frac_mean(Nmice+1) = metric_struct.(model).frac_mle(Nmice+1);
metric_struct.(model).std(Nmice+1) = std(metric_struct.(model).mle(1:Nmice))/sqrt(Nmice);
metric_struct.(model).frac_std(Nmice+1) = std(metric_struct.(model).frac_mle(1:Nmice))/sqrt(Nmice);

end
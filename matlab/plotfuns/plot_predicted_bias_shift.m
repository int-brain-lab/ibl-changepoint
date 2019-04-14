function bias_shift = plot_predicted_bias_shift(example_mice,frac_flag,avg_mouse_flag)
%PLOT_PREDICTED_BIAS_SHIFT Plot results of analysis of optimal bias shift.

if nargin < 2 || isempty(frac_flag)
    frac_flag = false;
end

if nargin < 3 || isempty(avg_mouse_flag)
    avg_mouse_flag = true;
end

fontsize = 14;

if nargin < 1 || isempty(example_mice)
    example_mice = {'CSHL_005','CSHL_007','IBL-T1','IBL-T4','ibl_witten_04','ibl_witten_05'};
end

Nmice = numel(example_mice);
Nsamples = 1e3;

for iMouse = 1:Nmice
    temp = load([example_mice{iMouse} '_bias_shift.mat'],'bias_shift');
    
    models = fields(temp.bias_shift)';
    for iModel = 1:numel(models)
        model_shift = temp.bias_shift.(models{iModel});
        
        % Generate posterior over bias shift for the data
        if strcmp(models{iModel},'data') && numel(model_shift) == 1
            temp2 = load([example_mice{iMouse} '_fits.mat'],'modelfits');
            idx = find(cellfun(@(p) strcmp(p.model_name,'psychofun'),temp2.modelfits.params),1);
            theta_rnd = get_posterior_samples(temp2.modelfits.params{idx},Nsamples);
            if ~isempty(theta_rnd)
                sampled_shift = theta_rnd(:,3) - theta_rnd(:,1);            
                model_shift = [model_shift,sampled_shift'];
            end
        end
        
        if strcmp(models{iModel},'data')
            data_shift = model_shift;
        end
        
        bias_shift.(models{iModel}).mle(iMouse) = model_shift(1);
        bias_shift.(models{iModel}).frac_mle(iMouse) = ...
            bias_shift.data.mle(iMouse)/bias_shift.(models{iModel}).mle(iMouse);
                
        if numel(model_shift) > 1
            bias_shift.(models{iModel}).mean(iMouse) = mean(model_shift(2:end));
            bias_shift.(models{iModel}).std(iMouse) = std(model_shift(2:end));
            
            rdata = randi(numel(data_shift)-1,[Nsamples,1])+1;
            rmodel = randi(numel(model_shift)-1,[Nsamples,1])+1;           
            frac = data_shift(rdata)./model_shift(rmodel);
            
            bias_shift.(models{iModel}).frac_mean(iMouse) = mean(frac);
            bias_shift.(models{iModel}).frac_std(iMouse) = std(frac);                        
        else
            bias_shift.(models{iModel}).mean(iMouse) = NaN;
            bias_shift.(models{iModel}).std(iMouse) = NaN;            
            bias_shift.(models{iModel}).frac_mean(iMouse) = NaN;
            bias_shift.(models{iModel}).frac_std(iMouse) = NaN;                        
        end        
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

if frac_flag
    model_offset = 0;
else
    model_offset = 0;
end

if avg_mouse_flag
    for iModel = 1:numel(models)
        bias_shift.(models{iModel}).mle(Nmice+1) = mean(bias_shift.(models{iModel}).mle(1:Nmice));
        bias_shift.(models{iModel}).frac_mle(Nmice+1) = mean(bias_shift.(models{iModel}).frac_mle(1:Nmice));
        bias_shift.(models{iModel}).mean(Nmice+1) = bias_shift.(models{iModel}).mle(Nmice+1);
        bias_shift.(models{iModel}).frac_mean(Nmice+1) = bias_shift.(models{iModel}).frac_mle(Nmice+1);
        bias_shift.(models{iModel}).std(Nmice+1) = std(bias_shift.(models{iModel}).mle(1:Nmice))/sqrt(Nmice);
        bias_shift.(models{iModel}).frac_std(Nmice+1) = std(bias_shift.(models{iModel}).frac_mle(1:Nmice))/sqrt(Nmice);
    end
    Nmice = Nmice + 1;
    example_mice{end+1} = 'average';
    colors(Nmice,:) = [0 0 0];
end

yy_max = 0;

for iModel = 1+model_offset:numel(models)
    tab = bias_shift.(models{iModel});
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

xlim([0,numel(models)+1-model_offset]);
if frac_flag
    ylim([0 1.2]);
    set(gca,'YTick',[0 0.5 1]);
    ystring = 'Fraction of optimal bias shift';
else
    ylim([0 ceil(yy_max*10)/10]);
    set(gca,'YTick',0:0.05:1);
    ystring = 'Optimal bias shift (contrast units)';
end

box off;
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'XTick',1:numel(models)-model_offset);
set(gca,'XTickLabel',modeltick);
set(gca,'FontSize',fontsize-2);

xlabel('Model used to predict optimal bias shift','FontSize',fontsize);
ylabel(ystring,'FontSize',fontsize);

for iMouse = 1:Nmice
    mice_names{iMouse} = example_mice{iMouse};
    mice_names{iMouse}(mice_names{iMouse} == '_') = '-';
end
hl = legend(h,mice_names{:});
set(hl,'Box','off','FontSize',fontsize,'Location','SouthEast');


end
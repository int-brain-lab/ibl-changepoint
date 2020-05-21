function plot_optimal_performance(data,model)
%PLOT_OPTIMAL_PERFORMANCE Plot actual performance vs. model-based performance

if nargin < 2 || isempty(model); model = 'changepoint_contrastnoise_runlength_probs'; end

for iData = 1:numel(data)
    data1 = data{iData};        
    if ischar(data1)
        dataset = read_data_from_csv(data1);
    else
        dataset = data1;
    end
    
    params = load_model_fit(data1,model);
%    [~,output] = nllfun([],params,dataset);
%    resp_model = output.resp_model;   
%    pcorrect(1,iData) = mean(resp_model .* (dataset.tab(:,5) < 0) + (1 - resp_model) .* (dataset.tab(:,5) >= 0));    

    
    pcorrect(1,iData) = calculate_pcorrect(model,params.theta,dataset);
    
    % mean(resp_model .* (dataset.tab(:,5) < 0) + (1 - resp_model) .* (dataset.tab(:,5) >= 0));    
    theta = params.theta(1:4);
    pcorrect(2,iData) = calculate_pcorrect('changepoint_contrastnoise',theta,dataset);
        
end

fontsize = 18;
axesfontsize = 14;

subplot(1,2,1);

scatter(pcorrect(1,:),pcorrect(2,:),'k'); 
hold on; 
plot([0.5 1],[0.5 1],'k--')
axis square;
box off;
    
set(gcf,'Color','w');
set(gca,'TickDir','out','FontSize',axesfontsize);
set(gca,'XTick',0.5:0.1:1,'YTick',0.5:0.1:1);

xlabel('Change-point observer (flexible)','FontSize',fontsize);
ylabel('Change-point observer (ideal)','FontSize',fontsize);
title('Average reward','FontSize',fontsize);

subplot(1,2,2);
histogram(pcorrect(2,:) - pcorrect(1,:),20);
hold on;

xlabel('\Delta Average reward','FontSize',fontsize);
ylabel('Count','FontSize',fontsize);
box off;
set(gca,'TickDir','out','FontSize',axesfontsize);
set(gca,'XTick',-0.1:0.02:0.1,'YTick',0:5:100);
    
yy = get(gca,'Ylim');
plot([0 0],yy,'k--')


end


function pc = calculate_pcorrect(model,theta,data)

params = params_new(model,data);
params = setup_params(theta,params);
[~,output] = nllfun([],params,data);
resp_model = output.resp_model;
pc = mean(resp_model .* (data.tab(:,5) < 0) + (1 - resp_model) .* (data.tab(:,5) >= 0));

end
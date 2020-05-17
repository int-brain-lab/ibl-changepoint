% PLOT GROUP PARAMETERS

batch_compute_parameters;

fontsize = 18;
axesfontsize = 14;

true_values = [60 20 0.2 0.8];

nbins = 50;
xmax = [20,10,1,1];

xx{1} = linspace(1,xmax(1),min(xmax(1),nbins));
xx{2} = linspace(1,xmax(2),min(xmax(2),nbins));
xx{3} = linspace(0,1,nbins);
xx{4} = linspace(0,1,nbins);

params_list = {'runlength_tau','runlength_min','prob_low','prob_high'};
logflags = [1 1 0 0];
Np = numel(params_list);

tt = zeros(size(theta,1),Np);

for iParam = 1:Np
    idx = find(strcmp(params.names,params_list{iParam}),1);
    tt(:,iParam) = theta(:,idx);
    if logflags(iParam); tt(:,iParam) = exp(tt(:,iParam)); end    
end

for iParam = 1:Np
    subplot(2,2,iParam);
    hist(tt(:,iParam),xx{iParam});
    h = findobj(gca,'Type','patch');
    h(1).FaceColor = [0 0.5 0.5];
    h(1).EdgeColor = 'w';
    box off;
    
    xtext = params_list{iParam};
    switch xtext
        case 'lapse_rate'; xtext = 'lapse rate';
        case 'runlength_tau'; xtext = 'block length time constant \tau'; xlim([0 xmax(1)]);
        case 'runlength_min'; xtext = 'minimum block length'; xlim([0 xmax(2)]);
        case 'prob_low'; xtext = 'P(Left) in Right block'; xlim([0 1]);
        case 'prob_high'; xtext = 'P(Left) in Left block'; xlim([0 1]);
        otherwise
            xtext(xtext == '_') = '-';
    end
    
    xlabel(xtext,'Fontsize',fontsize);
    ylabel('Counts');
    
    set(gca,'TickDir','out');
    set(gcf,'Color','w');
    set(gca,'Fontsize',axesfontsize);
    
    yy = get(gca,'Ylim');
    hold on;
    plot(true_values(iParam)*[1,1],yy,'--k','LineWidth',1);
    
end
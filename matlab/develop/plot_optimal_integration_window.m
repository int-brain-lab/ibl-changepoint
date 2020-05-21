% PLOT "OPTIMAL" INTEGRATION WINDOW

model_name = 'exponential_contrastnoise';
batch_compute_parameters;

fontsize = 18;
axesfontsize = 14;

true_values = [5.9248];

nbins = 50;
xmax = [10];

xx{1} = linspace(1,xmax(1),min(xmax(1),nbins));

params_list = {'runlength_tau'};
Np = numel(params_list);
logflags = 1;

tt = zeros(size(theta,1),Np);

for iParam = 1:Np
    idx = find(strcmp(params.names,params_list{iParam}),1);
    tt(:,iParam) = theta(:,idx);
    if logflags(iParam); tt(:,iParam) = exp(tt(:,iParam)); end
    idx_param(iParam) = idx;
end

figure(1);
for iParam = 1:Np
    hist(tt(:,iParam),xx{iParam});
    h = findobj(gca,'Type','patch');
    h(1).FaceColor = [0 0.5 0.5];
    h(1).EdgeColor = 'w';
    box off;
    
    xtext = params_list{iParam};
    switch xtext
        case 'runlength_tau'; xtext = 'integration time constant \tau'; xlim([0 xmax(1)]);
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

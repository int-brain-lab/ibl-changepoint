% PLOT GROUP PARAMETERS

batch_compute_parameters;

fontsize = 18;
axesfontsize = 14;

true_values = [60 20 0.2 0.8];

nbins = 50;
xx{1} = linspace(1,200,nbins);
xx{2} = linspace(1,30,nbins);
xx{3} = linspace(0,1,nbins);
xx{4} = linspace(0,1,nbins);

tt = theta(:,8:11);
tt(:,[1 2]) = exp(tt(:,[1 2]));
for i = 1:4
    subplot(2,2,i);
    hist(tt(:,i),xx{i});
    h = findobj(gca,'Type','patch');
    h(1).FaceColor = [0 0.5 0.5];
    h(1).EdgeColor = 'w';
    box off;
    
    xtext = params.names{i+7};
    switch xtext
        case 'lapse_rate'; xtext = 'lapse rate';
        case 'runlength_min'; xtext = 'minimum block length'; xlim([0 30]);
        case 'runlength_tau'; xtext = 'block length time constant \tau'; xlim([0 200]);
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
    plot(true_values(i)*[1,1],yy,'--k','LineWidth',1);
    
end
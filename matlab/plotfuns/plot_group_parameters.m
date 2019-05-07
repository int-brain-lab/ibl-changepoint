% PLOT GROUP PARAMETERS

batch_compute_parameters;

fontsize = 18;
axesfontsize = 14;

true_values = [NaN 60 20 0.2 0.8];

tt = theta(:,7:11);
tt(:,[2 3]) = exp(tt(:,[2 3]));
for i = 1:5
    subplot(2,3,i);
    hist(tt(:,i),50);
    h = findobj(gca,'Type','patch');
    h(1).FaceColor = [0 0.5 0.5];
    h(1).EdgeColor = 'w';
    box off;
    
    xtext = params.names{i+6};
    switch xtext
        case 'lapse_rate'; xtext = 'lapse rate';
        case 'runlength_min'; xtext = 'minimum block length';
        case 'runlength_tau'; xtext = 'block length time constant \tau';
        case 'prob_low'; xtext = 'P(Left) in Right block'; xlim([0 0.5]);
        case 'prob_high'; xtext = 'P(Left) in Left block'; xlim([0.5 1]);
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
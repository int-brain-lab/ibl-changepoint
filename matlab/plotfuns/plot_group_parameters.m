% PLOT GROUP PARAMETERS

batch_compute_parameters;

fontsize = 18;
axesfontsize = 14;

true_values = [NaN 60 20 0.2 0.8];

tt = theta(:,7:11);
tt(:,[2 3]) = exp(tt(:,[2 3]));
for i = 1:5
    subplot(2,3,i);
    hist(tt(:,i),100);
    h = findobj(gca,'Type','patch');
    h(1).FaceColor = [0 0.5 0.5];
    h(1).EdgeColor = 'w';
    box off;
    
    xtext = params.names{i+6};
    xtext(xtext == '_') = '-';
    xlabel(xtext,'Fontsize',fontsize);
    ylabel('Counts');
    
    set(gca,'TickDir','out');
    set(gcf,'Color','w');
    set(gca,'Fontsize',axesfontsize);
    
    yy = get(gca,'Ylim');
    hold on;
    plot(true_values(i)*[1,1],yy,'--k','LineWidth',1);
    
end
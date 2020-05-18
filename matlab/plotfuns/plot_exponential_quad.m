function plot_exponential_quad(params)

fontsize = 18;
axesfontsize = 14;

xx = linspace(0,59,100);

for iParam = 1:numel(params)
    
    tau = params{iParam}.runlength_tau;
    tau_quad = params{iParam}.tau_quad;
    tau_quadmu = params{iParam}.tau_quadmu;
        
    tau_y(iParam,:) = max(1,tau + tau_quad.*(log(xx+1) - log(tau_quadmu)).^2);
end

for ii = 1:size(tau_y,1)
    plot(xx,tau_y(ii,:),'-','Color',0.5*[1 1 1],'LineWidth',0.5); hold on;
end
plot(xx,mean(tau_y,1),'k-','LineWidth',3);

% yy = nanmean(tau_y,1);
% yy_serr = nanstd(tau_y,[],1);
% yyerr_down = yy - 1.96*yy_serr;
% yyerr_up = yy + 1.96*yy_serr;
% xxerr = [xx, fliplr(xx)];
% yyerr = [yyerr_down, fliplr(yyerr_up)];
% 
% fill(xxerr, yyerr,'k','FaceAlpha',0.5,'LineStyle','none'); hold on;
% h(1) = plot(xx,yy,'k-','LineWidth',2);

set(gcf,'Color','w');
set(gca,'TickDir','out','FontSize',axesfontsize);
box off;
xlabel('Trials from block start','FontSize',fontsize);
ylabel('\tau estimated from data','FontSize',fontsize);

xlim([0,60]);
ylim([1,30]);
set(gca,'XTick',0:20:100,'YTick',5:5:100);




end
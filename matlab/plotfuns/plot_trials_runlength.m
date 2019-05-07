function plot_trials_runlength(data,params,titlestr,xlims)

if nargin < 4 || isempty(xlims)
    first_session = data.tab(:,2) == data.tab(1,2);
    xlims = [0,sum(first_session)];
end

nx = 9;

subplot(4,nx,[1:nx-1,nx+(1:nx-1),2*nx+(1:nx-1)]);
plotchoice_flag = false;
plotbycontrast_flag = true;
plot_trials(data,params,titlestr,plotchoice_flag,plotbycontrast_flag);
xlabel('');
set(gca,'XTickLabel',[]);
xlim(xlims);
set(gca,'YTick',[0 0.2 0.5 0.8 1]);

subplot(4,nx,3*nx+(1:nx-1));
plot_runlength_posterior(data,params,'Posterior over run lengths');
xlim(xlims);
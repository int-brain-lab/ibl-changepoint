function plot_trials_runlength(data,params,titlestr,xlims,rl_flag)

if nargin < 4 || isempty(xlims)
    first_session = data.tab(:,2) == data.tab(1,2);
    xlims = [0,sum(first_session)];
end

if nargin < 5 || isempty(rl_flag); rl_flag = true; end

nx = 9;
fontsizes = [24 18];

if rl_flag; subplot(4,nx,[1:nx-1,nx+(1:nx-1),2*nx+(1:nx-1)]); end
plotchoice_flag = false;
plotbycontrast_flag = true;
plot_trials(data,params,titlestr,plotchoice_flag,plotbycontrast_flag,fontsizes);
xlim(xlims);
set(gca,'YTick',[0 0.2 0.5 0.8 1]);

if rl_flag
    xlabel('');
    set(gca,'XTickLabel',[]);
    
    subplot(4,nx,3*nx+(1:nx-1));
    plot_runlength_posterior(data,params,'Posterior probability over run lengths',fontsizes);
    xlim(xlims);
    colorbar;
end
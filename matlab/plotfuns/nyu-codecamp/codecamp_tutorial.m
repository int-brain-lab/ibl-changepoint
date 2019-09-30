% NYC Code Camp tutorial - September 03, 2019

% Load example dataset (only training)
dataset = 'NYU-01_endtrain';
data = read_data_from_csv(dataset);
params = load_model_fit(dataset,'psychofun');

%% Make basic plot
axesfontsize = 22; fontsize = 28;

plot_fit(data,params,[],false);
legend off;
ylim([0 1]);
set(gca,'XTick',[-1 -0.5 -0.25 -0.125 0 0.125 0.25 0.5 1]);
set(gca,'YTick',[0,0.5,1]);
set(gca,'FontSize',axesfontsize);
xlabel('Signed contrast','FontSize',fontsize);
ylabel('Fraction choice ''Right''','FontSize',fontsize);
title('Training sessions performance (example mouse)','FontSize',fontsize);


%% Plot posterior
params.logspace = logical([0 1 0 0]);
plot_posterior(params);

dataset1 = 'NYU-01_sess1_unbiased';
data1 = read_data_from_csv(dataset1);
params1 = fit_model('psychofun',data1,[1,1],0,1);

params1.logspace = logical([0 1 0 0]);
plot_posterior(params1);

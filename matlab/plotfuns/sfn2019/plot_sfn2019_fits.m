%PLOT_SFN2019_FITS Plot data and fits for one representative mouse

close all;
data_list = get_mice_list('strict_sep2019');
model_list = {'changepoint_contrastnoise', 'changepoint_contrastnoise_runlength_probs', ...
    'exponential_contrastnoise', 'omniscient_contrastnoise_fixedfreeprior'};
tab = collect_model_comparison(data_list, model_list);
mouse_name = 'IBL-T1';
data = read_data_from_csv(mouse_name);

fontsizes = 32;
for iModel = 1:4
    iPlot = iModel; if iPlot > 2; iPlot = iPlot + 1; end
    subplot(2,3,iPlot);
    params = load_model_fit(mouse_name,tab.models{iModel});
    plot_changepoints(data,params,20,[],'l',false,fontsizes)
    xlabel(''); ylabel('');
    legend off;
end

subplot(2,3,3);
plot_changepoints(data,params,20,[],'l',false,fontsizes)
xlabel(''); ylabel(''); axis off;
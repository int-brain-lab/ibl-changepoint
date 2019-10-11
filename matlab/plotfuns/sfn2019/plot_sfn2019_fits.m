%PLOT_SFN2019_FITS Plot data and fits for one representative mouse

data_list = get_mice_list('strict_sep2019');
model_list = {'changepoint_contrastnoise', 'changepoint_contrastnoise_runlength_probs', ...
    'exponential_contrastnoise', 'omniscient_contrastnoise_fixedfreeprior'};
tab = collect_model_comparison(data_list, model_list);

mouse_name = 'IBL-T1';
for iModel = 1:4
    params = load_model_fit(mouse_name,tab.models{4});
    data = read_data_from_csv(mouse_name);
    plot_changepoints(data,params,20,[],'l')
end

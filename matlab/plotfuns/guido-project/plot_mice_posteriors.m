% Plot results from model fits

mice_names = {'ZM_1895','ZM_2102','ZM_2108'};
mice_list = get_mice_list('guido');
model_name = 'exponential_contrastnoise';
color_list = [0 0.8 0; 0 0.8 0; 0.3 1 0.3; 0 0 0.8; 0 0 0.8; 0.3 0.3 1; 0.8 0 0; 0.8 0 0; 1 0.3 0.3];
symbol_list = {'-',':','-','-',':','-','-',':','-'};
fig_pos = [1,41,1920,963];

for iMouse = 1:numel(mice_names)
    figure(iMouse);
    session_idx = strmatch(mice_names{iMouse},mice_list);
    [X{iMouse},pnames{iMouse},params{iMouse},data_names{iMouse}] = plot_posterior_comparison(mice_list(session_idx),model_name,color_list,symbol_list);
    set(gcf,'Position',fig_pos);
    drawnow;
    saveas(gcf,[mice_names{iMouse} '_posterior.png'],'png');
end

date_str = lower(date);
date_str(date_str == '-') = '';
filename = ['guido_analysis_' date_str '.mat'];
save(filename,'data_names','X','params','pnames','mice_list');
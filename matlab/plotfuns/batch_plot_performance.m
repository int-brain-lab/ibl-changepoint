% PLOT all mice performance

if ~exist('mice_list','var') || isempty(mice_list)
    mice_list = get_mice_list();
end

training_data = 'endtrain';

for iMouse = 1:numel(mice_list)
    subplot(4,5,iMouse);
    mouse_name = mice_list{iMouse};
    fprintf('%s ',mouse_name);
    
    load([mouse_name '_' training_data '_fits.mat']);
    data_all = read_data_from_csv(mouse_name);
    % suffix = '';
    % suffix = '_biasedlapse';
    suffix = '_altnoise';
    
    params = params_new(['changepoint' suffix],data_all);
    idx = find(cellfun(@(p) strcmp(p.model_name,['omniscient' suffix]),modelfits.params),1);
    theta_rnd = modelfits.params{idx}.theta;
    params1 = setup_params(theta_rnd(1,:),params);
    
    titlestr = mouse_name;
    titlestr(titlestr=='_') = '-';
    plot_performance(data_all,params1,titlestr);
    legend off;
end
fprintf('\n');

mypath = which('savefigure.m');
savefigure([fileparts(mypath) filesep() 'predicted_performance_with_' training_data]);

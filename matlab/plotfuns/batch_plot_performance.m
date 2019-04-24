% PLOT all mice performance

if ~exist('mice_list','var') || isempty(mice_list)
    % mice_list = {'CSHL_005','CSHL_007','IBL-T1','IBL-T4','ibl_witten_04','ibl_witten_05'}; % Example mice
    mice_list = {'CSHL_003','CSHL_005','CSHL_007','CSHL_008','CSHL_010','IBL-T1','IBL-T4','ZM_1084','ZM_1085','ZM_1086','ZM_1091','ZM_1092','ZM_1093','ZM_1097','ZM_1098','ibl_witten_04','ibl_witten_05','ibl_witten_06'};
end

training_data = 'endtrain';

for iMouse = 1:numel(mice_list)
    subplot(4,5,iMouse);
    mouse_name = mice_list{iMouse};
    fprintf('%s ',mouse_name);
    
    load([mouse_name '_' training_data '_fits.mat']);
    data_all = read_data_from_csv(mouse_name);
    % suffix = '_biasedlapse';
    suffix = '';
    
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
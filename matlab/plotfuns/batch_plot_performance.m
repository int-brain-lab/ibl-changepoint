% PLOT all mice performance

if ~exist('mice_list','var') || isempty(mice_list)
    mice_list = {'CSHL_005','CSHL_007','IBL-T1','IBL-T4','ibl_witten_04','ibl_witten_05'}; % Example mice
    % mice_list = {'ibl_witten_06', 'IBL-T4', 'ibl_witten_05', 'ibl_witten_04', 'IBL-T1', 'ZM_1092', 'CSHL_007', 'ZM_1093', 'CSHL_010', 'CSHL_008', 'ZM_1085', 'CSHL_003', 'ZM_1084', 'ZM_1098', 'ZM_1091', 'CSHL_005', 'ZM_1097', 'ZM_1086'};
end

for iMouse = 1:numel(mice_list)
    subplot(2,3,iMouse);
    mouse_name = mice_list{iMouse};    
    load([mouse_name '_unbiased_fits.mat']);
    data_all = read_data_from_csv(mouse_name);
    % suffix = '_biasedlapse';
    suffix = '';
    
    params = params_new(['changepoint' suffix],data_all);
    idx = find(cellfun(@(p) strcmp(p.model_name,['omniscient' suffix]),modelfits.params),1);
    theta_rnd = modelfits.params{idx}.theta;
    params1 = setup_params(theta_rnd(1,:),params);
    plot_performance(data_all,params1)
end
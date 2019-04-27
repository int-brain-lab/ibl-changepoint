% PLOT all mice performance

if ~exist('mice_list','var') || isempty(mice_list)
    mice_list = get_mice_list();
end

test_model_base = 'changepoint';
% test_model_base = 'omniscient';
% test_model_base = 'omniscient_fixedprior';
% training_data = 'unbiased';
training_data = 'endtrain';

output = [];

for iMouse = 1:numel(mice_list)
    subplot(4,5,iMouse);
    mouse_name = mice_list{iMouse};
    fprintf('%s ',mouse_name);
    
    data_all = read_data_from_csv(mouse_name);
    % suffix = '';
    % suffix = '_biasedlapse';
    suffix = '_doublenoise';

    % Get training model parameters
    params_train = load_model_fit([mouse_name '_' training_data],['omniscient' suffix]);
    theta_train = params_train.theta;
    
    test_model = [test_model_base suffix];    
    params = params_new(test_model,data_all);
    params1 = setup_params(theta_train(1,:),params);
    
    titlestr = mouse_name;
    titlestr(titlestr=='_') = '-';
    output{iMouse} = plot_performance(data_all,params1,titlestr);
    legend off;
end
fprintf('\n');

mypath = which('savefigure.m');
savefigure([fileparts(mypath) filesep() 'predicted_performance_with_' training_data '_and_' test_model]);

function params = params_transfer(train_params,test_model,test_data,train_theta)
%PARAMS_TRANSFER Transfer parameters between two models.

if nargin < 4; train_theta = []; end

% Update training model if a parameter vector is explicitly passed
if ~isempty(train_theta)
    train_params = setup_params(train_theta,train_params);
end

% Read dataset if passed as string
if ischar(test_data)
    test_data = read_data_from_csv(test_data);
end

% Create new model
params = params_new(test_model,test_data);
params = setup_params(train_params.theta,params);

% Ensure that lapse is transferred correctly
if contains(test_model,'lapse')
    params.lapse_rate = train_params.lapse_rate;
    params.lapse_bias = train_params.lapse_bias;    
end




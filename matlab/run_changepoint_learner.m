%RUN_CHANGEPOINT_LEARNER Script to execute and save a changepoint learner.

if ~exist('mouse_name','var')
    error('Need to specify a mouse name.')
end

tmax = 1e4; % Maximum number of trials considered

data = read_data_from_csv(mouse_name); % Load mouse data

% Load perceptual decision-making parameters from flexible model fit
params1 = load_model_fit(mouse_name,'changepoint_nakarushton_runlength_probs');
theta = params1.theta(1:6);

params = params_new('changelearn_nakarushton',data);
params = setup_params(theta,params);
params.Ngrid = 10;

% Compute model
[nLL,params.output] = changelearn_bayesian_nll(params,data,tmax);

% Save negative log likelihood
params.mle_fits.x0 = theta;
params.mle_fits.x = theta;
params.mle_fits.nll = nLL;
params.mle_fits.x_best = theta;
params.mle_fits.nll_best = theta;
params.mle_fits.ndata = size(data.tab,1);

% Save fit
save_model_fit(mouse_name,params);
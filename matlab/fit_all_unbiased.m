if ~exist('example_mice','var') || isempty(example_mice)
    example_mice = {'CSHL_005','CSHL_007','IBL-T1','IBL-T4','ibl_witten_04','ibl_witten_05'};
end

clear train_models test_models;

Nsamples = 20;  % Approximate posterior samples for model predictions

train_models{1} = 'psychofun';
train_models{2} = 'omniscient';
train_models{3} = 'omniscient_lapse';
train_models{4} = 'omniscient_biasedlapse';
% train_models{5} = 'omniscient_softmax';

test_models{1} = [];
test_models{2} = 'changepoint';
test_models{3} = 'changepoint_lapse';
test_models{4} = 'changepoint_biasedlapse';
% test_models{5} = 'changepoint_softmax';

% Psychometric curve at zero contrast for a given block
psycho0 = @(psy,block) psychofun(0,[psy.psycho_mu(block),psy.psycho_sigma(block),psy.psycho_gammalo(block),psy.psycho_gammahi(block)]);

for iMouse = 1:numel(example_mice)
    
    % Fit psychometric curves for all blocks
    modelfits_psy = batch_model_fit('psychofun',example_mice{iMouse},3,1,0);
    idx = find(cellfun(@(p) strcmp(p.model_name,'psychofun'),modelfits_psy.params),1);
    psy_data = modelfits_psy.params{idx};
    bias_shift.data(iMouse,1) = psy_data.psycho_mu(3) - psy_data.psycho_mu(1);
    
    % Psychometric curves at zero contrast for biased blocks
    pRblock = psycho0(psy_data,1);
    pLblock = psycho0(psy_data,3);
    bias_prob.data.RightBlock(iMouse,1) = pRblock;
    bias_prob.data.LeftBlock(iMouse,1) = pLblock;
    bias_prob.data.MatchProbability(iMouse,1) = 0.5*(pRblock + 1 - pLblock);
    
    % Fit all models on unbiased blocks only
    modelfits = batch_model_fit(train_models,[example_mice{iMouse} '_unbiased'],5,1,0);
    
    % Simulate change-point models on all blocks w/ parameters from unbiased blocks
    data_all = read_data_from_csv(example_mice{iMouse});    % Get mouse data
    
    psy_model_mle = []; gendata_mle = [];
    
    for iModel = 1:numel(test_models)
        if isempty(test_models{iModel}); continue; end
        
        params = params_new(test_models{iModel},data_all);
        idx = find(cellfun(@(p) strcmp(p.model_name,train_models{iModel}),modelfits.params),1);

        % First sample is maximum-likelihood estimate (MLE)
        theta_rnd = modelfits.params{idx}.theta;
        
        % Generate additional samples from the approximate posterior
        theta_rnd = [theta_rnd; ...
            get_posterior_samples(modelfits.params{idx},Nsamples)];
                
        theta_rnd
        
        for iSample = 1:size(theta_rnd,1)
            
            iSample
            
            % Setup parameter struct for changepoint model
            params1 = setup_params(theta_rnd(iSample,:),params);

            % Generate data
            nreps = 10;
            gendata = model_gendata(params1,data_all,nreps);

            % Fit psychometric curve to model-generated data
            psy_model = fit_model('psychofun',gendata,1,0,1,psy_data);
            bias_shift.(test_models{iModel})(iMouse,iSample) = ...
                psy_model.psycho_mu(3) - psy_model.psycho_mu(1);
            
            % Psychometric curves at zero contrast for biased blocks
            pRblock = psycho0(psy_model,1);
            pLblock = psycho0(psy_model,3);
            bias_prob.(test_models{iModel}).RightBlock(iMouse,iSample) = pRblock;
            bias_prob.(test_models{iModel}).LeftBlock(iMouse,iSample) = pLblock;
            bias_prob.(test_models{iModel}).MatchProbability(iMouse,iSample) = 0.5*(pRblock + 1 - pLblock);
                        
            % Store maximum-likelihood psychometric function and data
            if iSample == 1
                gendata_mle{iModel} = gendata;
                psy_model_mle{iModel} = psy_model;
            end
        end
    end
    
    save([example_mice{iMouse} '_bias_shift.mat'],'bias_shift','bias_prob','psy_model_mle','gendata_mle','train_models','test_models');
end

% Make plots
close all;
subplot(2,1,1);
bias_shift = plot_predicted_bias_shift([],0);
subplot(2,1,2);
plot_predicted_bias_shift([],1);
mypath = which('savefigure.m');
savefigure([fileparts(mypath) filesep() 'bias_shift']);
close all;
plot_predicted_psychometric_shift(example_mice);
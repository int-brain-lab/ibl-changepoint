if ~exist('mice_list','var') || isempty(mice_list)
    % Example mice
    % mice_list = {'CSHL_005','CSHL_007','IBL-T1','IBL-T4','ibl_witten_04','ibl_witten_05'};   
    mice_list = {'CSHL_003','CSHL_005','CSHL_007','CSHL_008','CSHL_010','IBL-T1','IBL-T4','ZM_1084','ZM_1085','ZM_1086','ZM_1091','ZM_1092','ZM_1093','ZM_1097','ZM_1098','ibl_witten_04','ibl_witten_05','ibl_witten_06'};
end

clear train_models test_models;

% Train on unbiased sessions of the biased blocks
% train_set = 'unbiased';

% Train on final three training sessions (before biased protocol)
train_set = 'endtrain';

Nsamples = 20;  % Approximate posterior samples for model predictions

train_models{1} = 'psychofun';
train_models{2} = 'omniscient';
train_models{3} = 'omniscient_lapse';
train_models{4} = 'omniscient_biasedlapse';
train_models{5} = 'omniscient_altnoise';
% train_models{6} = 'omniscient_softmax';

test_models{1} = [];
test_models{2} = 'changepoint';
test_models{3} = 'changepoint_lapse';
test_models{4} = 'changepoint_biasedlapse';
train_models{5} = 'omniscient_altnoise';
% test_models{6} = 'changepoint_softmax';

% Compute posterior distributions over parameters?
%compute_posteriors_flag = true;
compute_posteriors_flag = false;

% Psychometric curve at zero contrast for a given block
psycho0 = @(psy,block) psychofun(0,[psy.psycho_mu(block),psy.psycho_sigma(block),psy.psycho_gammalo(block),psy.psycho_gammahi(block)]);

for iMouse = 1:numel(mice_list)
    
    %% First, fit psychometric curve of all biased sessions
    
    % Fit psychometric curves for all blocks
    modelfits_psy = batch_model_fit('psychofun',mice_list{iMouse},3,compute_posteriors_flag,0);
    idx = find(cellfun(@(p) strcmp(p.model_name,'psychofun'),modelfits_psy.params),1);
    psy_data = modelfits_psy.params{idx};
    bias_shift.data(iMouse,1) = psy_data.psycho_mu(3) - psy_data.psycho_mu(1);
    
    % Psychometric curves at zero contrast for biased blocks
    pRblock = psycho0(psy_data,1);
    pLblock = psycho0(psy_data,3);
    bias_prob.data.RightBlock(iMouse,1) = pRblock;
    bias_prob.data.LeftBlock(iMouse,1) = pLblock;
    bias_prob.data.MatchProbability(iMouse,1) = 0.5*(pRblock + 1 - pLblock);
    
    %% Second, fit all training models on training data only
    
    train_filename = [mice_list{iMouse} '_' train_set];
    
    % Fit all models on training data
    modelfits = batch_model_fit(train_models,train_filename,5,compute_posteriors_flag,0);
    
    %% Third, simulate ideal change-point observer given training models
    
    % Simulate change-point models on all blocks w/ parameters from unbiased blocks
    data_all = read_data_from_csv(mice_list{iMouse});    % Get mouse data
    
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
    
    save([mice_list{iMouse} '_bias_shift_' train_set '.mat'],'bias_shift','bias_prob','psy_model_mle','gendata_mle','train_models','test_models');
end

% Make plots
close all;
metrics = {'bias_shift','pmatch'};

subplot(2,1,1);
bias_shift = plot_predicted_bias_shift(metrics{1},[],0);
subplot(2,1,2);
plot_predicted_bias_shift(metrics{1},[],1);
mypath = which('savefigure.m');
savefigure([fileparts(mypath) filesep() metrics{1}]);
close all;
[~,pmatch] = plot_predicted_bias_shift(metrics{2},[],0);
mypath = which('savefigure.m');
savefigure([fileparts(mypath) filesep() metrics{2}]);
close all;
plot_predicted_psychometric_shift(mice_list);
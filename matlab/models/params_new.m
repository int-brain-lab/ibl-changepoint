function params = params_new(model_name,data)
%PARAMS_NEW Return new parameter struct for a given model.

if isempty(model_name); model_name = 'changepoint_extra'; end

params.model_name = model_name;
params.model_startingfit = [];  % Model used to initialize fits
extra_features = [];

%% Properties of base models

idx = find(model_name=='_',1);
if isempty(idx); idx = numel(model_name)+1; end
base_model = model_name(1:idx-1);

switch base_model
    case 'psychofun'
        params.model_nLLfun = @psycho_nll;
        params.model_desc = 'Psychometric functions';
    case 'omniscient'
        params.model_nLLfun = @omniscient_bayesian_nll;
        params.model_desc = 'Omniscient Bayesian observer';
    case 'changepoint'        
        params.model_nLLfun = @changepoint_bayesian_nll;
        params.model_desc = 'Change-point Bayesian observer';
    otherwise
        error('Unknown base model.');
end        

%% Global parameters

if strcmp(base_model,'psychofun')
    % Psychometric function model
    
    Nprobs = numel(unique(data.p_true));
    
    % Psychometric function parameters
    params.psycho_mu = 0*ones(1,Nprobs);
    params.psycho_sigma = 0.1*ones(1,Nprobs);
    params.psycho_gammalo = 0.05*ones(1,Nprobs);
    params.psycho_gammahi = 0.05*ones(1,Nprobs);
    
    for iProb = 1:Nprobs
        params.names{iProb} = 'psycho_mu';
        params.names{Nprobs+iProb} = 'psycho_sigma';
        params.names{Nprobs*2+iProb} = 'psycho_gammalo';
        params.names{Nprobs*3+iProb} = 'psycho_gammahi';
    end
        
else
    % Bayesian observer models

    % Noise parameters (~discrimination threshold at given contrast)
    if contains(model_name,'altnoise')
        params.precision = [log(1/10^2), log(1/10^2)];
        params.precision_power = [1 1];
        for iParam = 1:numel(params.precision)
            params.names{iParam} = 'precision';
        end
        for iParam = 1:numel(params.precision)
            params.names{end+1} = 'precision_power';
        end
        extra_features{end+1} = 'asymmetric precision';
        
    else
        params.sigma_contrasts(1,:) = [0.05 0.15 0.45]; % data.contrasts_vec(data.contrasts_vec~=0);
        params.sigma = linspace(20,5,numel(params.sigma_contrasts));
        if contains(model_name,'doublenoise')
            params.sigma = repmat(params.sigma,[1 2]);
            extra_features{end+1} = 'left/right noise';
        end
        for iParam = 1:numel(params.sigma)
            params.names{iParam} = 'sigma';
        end
    end
    
    % Load noise parameters from another fit
    if contains(model_name,'loadnoise')
        
        if contains(model_name,'altnoise')
            noise_model_name = 'omniscient_altnoise';
        elseif contains(model_name,'doublenoise')
            noise_model_name = 'omniscient_doublenoise';
        else
            noise_model_name = 'omniscient';            
        end        
        
        % Load existing fits
        mypath = fileparts(mfilename('fullpath'));
        fits_path = [mypath filesep 'fits'];
        addpath(fits_path);
        
        idx_str = strfind(model_name,'loadnoise');
        str = model_name(idx_str+numel('loadnoise'):end);
        idx = find([str '_'] == '_',1);        
        noisefile = str(1:idx-1);

        params_noise = load_model_fit([data.name '_' noisefile],noise_model_name);
        
        if contains(model_name,'altnoise')
            params.precision = params_noise.precision;
            params.precision_power = params_noise.precision_power;
        elseif contains(model_name,'doublenoise')
            params.sigma = params_noise.sigma;
            params.sigma_poly_left = params_noise.sigma_poly_left;
            params.sigma_poly_right = params_noise.sigma_poly_right;
        else
            params.sigma = params_noise.sigma;
            params.sigma_poly = params_noise.sigma_poly;
        end
        
        % Remove SIGMA from the parameters
        params.names(cellfun(@(x)strcmp(x,'sigma'),params.names)) = [];
    end
    

    % Lapse rate
    params.lapse_rate = 0;
    
    % Lapse bias (probability of responding "Left" on a lapse)
    params.lapse_bias = 0.5;
    
    % Attentional shift (multiplicative/divisive factor to SIGMA for left/right)
    params.attention_factor = 1;

    % Softmax parameters for probabilistic mapping from posterior to choice
    params.softmax_eta = Inf;     % Deterministic choice, take argmax
    params.softmax_bias = 0;      % No bias by default
    
    if contains(model_name,'_biasedlapse')
        params.names{end+1} = 'lapse_rate';
        params.names{end+1} = 'lapse_bias';
        extra_features{end+1} = 'biased lapse';

    elseif contains(model_name,'_fixedbiasedlapse')
        extra_features{end+1} = 'fixed biased lapse';
        
        left_hi = 1-data.resp_correct(data.contrasts == 1 & data.S < 0);
        right_hi = data.resp_correct(data.contrasts == 1 & data.S > 0);        
        lapse_left = 1 - mean(left_hi);
        lapse_right = 1 - mean(right_hi);
        Nleft = numel(left_hi);
        Nright = numel(right_hi);
        
        params.lapse_rate = 1+mean(left_hi)-mean(right_hi);
        params.lapse_bias = mean(left_hi)/params.lapse_rate;
        
        params
        
    elseif contains(model_name,'_lapse')
        params.names{end+1} = 'lapse_rate';
        extra_features{end+1} = 'lapse';        
    end
    
    if contains(model_name,'_softmax')
        params.names{end+1} = 'softmax_eta'; 
        params.names{end+1} = 'softmax_bias';
        extra_features{end+1} = 'softmax';
    end
    
    if contains(model_name,'_attention')
        params.names{end+1} = 'attention_factor';
        extra_features{end+1} = 'attention';
    end
    
    if contains(model_name,'fixedprior')
        params.fixed_prior = 0.5;
        extra_features{end+1} = 'fixed prior';
    end
    
    %% Parameters for change-point observer only

    if strncmp(model_name,'changepoint',numel('changepoint'))
        % Min and max run lengths
        params.runlength_min = 20;
        params.runlength_max = 100;

        % Function handle to change-point prior
        params.runlength_tau = 60;
        params.runlength_prior = ['@(t) exp(-t/' num2str(params.runlength_tau,'%.8f') ')'];

        % Probability levels in the session
        p_true_vec = unique(data.p_true);
        
        % Set transition states
        switch numel(p_true_vec)
            case 3  % Start with biased block, then alternating
                params.p_true_vec = [p_true_vec(1) p_true_vec(3) p_true_vec(2) 0.5];
                params.p0 = [0 0 0 1];
                params.Tmat = [0 1 0 0; 1 0 0 0; 0.5 0.5 0 0; 0 0 1 0];
            case 2  % Only alternating blocks
                params.p_true_vec = [p_true_vec(1) p_true_vec(2)];
                params.p0 = [0.5 0.5];
                params.Tmat = [0 1; 1 0];
            case 1  % Only one p -- must be the "unbiased" case
                params.p_true_vec = [0.2 0.8 0.5 0.5];
                params.p0 = [0 0 0 1];
                params.Tmat = [0 1 0 0; 1 0 0 0; 0.5 0.5 0 0; 0 0 1 0];
        end
        
        % Probability states (representing "probability of left stimulus")
        params.p_vec = params.p_true_vec;
        Nprobs = numel(params.p_vec);
        
        % Probabilities for low and high-probability blocks
        params.prob_low = min(params.p_vec);
        params.prob_high = max(params.p_vec);
        
%         % Probability states (representing "probability of left stimulus")
%         defaults.p_vec = [0.2,0.5,0.8];
%         Nprobs = numel(defaults.p_vec);
% 
%         % Starting belief over state
%         defaults.p0 = ones(1,Nprobs)/Nprobs;
% 
%         % Default transition matrix (equal-probability change)
%         defaults.Tmat = (ones(Nprobs,Nprobs) - eye(Nprobs)) / (Nprobs-1);

        % Beta hyperprior on observations after changepoint (bias towards L/R)
        params.beta_hyp = [0,0];      % No bias by default
        
        if contains(model_name,'_runlength')
            params.names{end+1} = 'runlength_tau';
            params.names{end+1} = 'runlength_min';
            extra_features{end+1} = 'runlengths';
        end
        if contains(model_name,'_freesym')
            params.names{end+1} = 'prob_low';
            extra_features{end+1} = 'symmetric probs';
        end
        if contains(model_name,'_probs')
            params.names{end+1} = 'prob_low';
            params.names{end+1} = 'prob_high';
            extra_features{end+1} = 'probs';
        end
    end
    
    % Assign extra features to model description
    if ~isempty(extra_features)
        extra_desc = ' (';
        for iFeat = 1:numel(extra_features)
            extra_desc = [extra_desc, extra_features{iFeat}];
            if iFeat < numel(extra_features)
                extra_desc = [extra_desc,', '];
            end
        end
        extra_desc = [extra_desc,')'];        
        params.model_desc = [params.model_desc,extra_desc];        
    end    
    
    % Assign starting fits
    
    if strcmp(base_model,'changepoint')
        params.model_startingfit = 'omniscient';
        if numel(model_name) > numel(base_model)
            params.model_startingfit = ...
                [params.model_startingfit, model_name(numel(base_model)+1:end)];
        end
    end

end
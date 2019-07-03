%FIT_PSYCHO_SESSIONS Fit psychometric functions separately on all sessions.

if ~exist('mice_list','var') || isempty(mice_list)
    mice_list = get_mice_list('criteria_jul2019');
end

Nsamples = 20; % Approximate posterior samples for model predictions
Nopts_psy = 5; % # optimization restarts for fitting psychometric curves

% Compute posterior distributions over parameters?
compute_posteriors_flag = true;
%compute_posteriors_flag = false;

empirical_list = mice_list; % Datasets for empirical Bayes

% Psychometric curve at zero contrast for a given block
psycho0 = @(psy,block) psychofun(0,[psy.psycho_mu(block),psy.psycho_sigma(block),psy.psycho_gammalo(block),psy.psycho_gammahi(block)]);

for iMouse = 1:numel(mice_list)
    
    mouse_name = mice_list{iMouse};
    
    %% Fit psychometric curve of all biased sessions
        
    % First, fit psychometric curves for all sessions
    psyfit_all = batch_model_fit('psychofun',mouse_name,Nopts_psy,compute_posteriors_flag,0,empirical_list);
    
    % Then, fit separately each session
    data = read_data_from_csv(mouse_name);
    sessions = unique(data.tab(:,2));
    
    % Take only session that have all block types
    for iSession = 1:numel(sessions)
        idx = data.tab(:,2) == sessions(iSession);
        p = unique(data.tab(idx,3));
        if numel(p) ~= 3;  sessions(iSession) = NaN; end
    end    
    sessions = sessions(isfinite(sessions));
    
    for iSession = 1:numel(sessions)
        psyfit_sessions{iSession} = batch_model_fit('psychofun', ...
            [mouse_name '_sess' num2str(sessions(iSession))],Nopts_psy,compute_posteriors_flag,0,empirical_list);
    end     
end
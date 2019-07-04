%FIT_PSYCHO_SESSIONS Fit psychometric functions separately on all sessions.

if ~exist('mice_list','var') || isempty(mice_list)
    mice_list = get_mice_list('strict_jul2019');
end

Nsamples = 20;              % Approximate posterior samples for model predictions
Nopts_psy = 5;              % # restarts for fitting psychometric curves
refit_flag = [false,true];  % Refit already existing fits?

% Compute posterior distributions over parameters?
compute_posteriors_flag = true;

empirical_list = get_mice_list('strict_jul2019'); % Datasets for empirical Bayes

empirical_endtrain = [];
for iMouse = 1:numel(empirical_list)
    empirical_endtrain{iMouse} = [empirical_list{iMouse} '_endtrain'];
end

for iMouse = 1:numel(mice_list)
    
    mouse_name = mice_list{iMouse};
    
    %% Fit psychometric curve of all biased sessions    
    data_name = mouse_name;

    % First, fit psychometric curves for all sessions
    psyfit_all = batch_model_fit('psychofun',data_name,Nopts_psy,compute_posteriors_flag,refit_flag,empirical_list);

    % Then, fit separately each session
    data = read_data_from_csv(data_name);
    sessions = unique(data.tab(:,2));

    % Take only session that have all block types
    for iSession = 1:numel(sessions)
        idx = data.tab(:,2) == sessions(iSession);
        p = unique(data.tab(idx,3));
        if numel(p) ~= 3;  sessions(iSession) = NaN; end
    end    
    sessions = sessions(isfinite(sessions));

    for iSession = 1:numel(sessions)
        fprintf('Fitting session %d / %d...\n',sessions(iSession),numel(sessions));
        psyfit_sessions{iSession} = batch_model_fit('psychofun', ...
            [data_name '_sess' num2str(sessions(iSession))],Nopts_psy,compute_posteriors_flag,refit_flag,empirical_list);
    end
    
    %% Fit psychometric curves of last three training sessions
    
    data_name = [mouse_name '_endtrain'];
    
    % First, fit psychometric curves for all sessions
    psyfit_all_endtrain = batch_model_fit('psychofun',data_name,Nopts_psy,compute_posteriors_flag,refit_flag,empirical_endtrain);
    
    % Then, fit separately each session
    data = read_data_from_csv(data_name);
    sessions = unique(data.tab(:,2));
        
    for iSession = 1:numel(sessions)
        fprintf('Fitting training session %d / %d...\n',sessions(iSession),numel(sessions));
        psyfit_sessions_endtrain{iSession} = batch_model_fit('psychofun', ...
            [data_name '_sess' num2str(sessions(iSession))],Nopts_psy,compute_posteriors_flag,refit_flag,empirical_endtrain);
    end


end
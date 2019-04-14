function [nLL,output] = psycho_nLL(params,data,sigma)
%PSYCHO_NLL Negative log likelihood of psychometric functions.

MIN_P = 1e-4;   % Minimum lapse/error probability

p_true_vec = unique(data.p_true);

PChatL = zeros(size(data.resp_obs));

% Fit each probability level separately
for i_prob = 1:numel(p_true_vec)
    psycho_mu = params.psycho_mu(i_prob);
    psycho_sigma = params.psycho_sigma(i_prob);
    
    idx_trial = data.p_true == p_true_vec(i_prob);
    cc = data.signed_contrasts(idx_trial);
    
    % Psychometric curve for choosing RIGHT
    p = 0.5*(1 + erf((cc-psycho_mu)/(sqrt(2)*psycho_sigma)));
    
    % Add left (low) and right (high) lapses
    p = params.psycho_gammalo(i_prob) + ...
        max(0,(1 - params.psycho_gammalo(i_prob) - params.psycho_gammahi(i_prob)))*p;
        
    % Probability or choice LEFT is 1 - psychometric curve
    PChatL(idx_trial) = 1 - p;
end

% Log probability of observer responses (ignore no-response trials)
log_PChat = log(PChatL).*(data.resp_obs == -1) + log(1-PChatL).*(data.resp_obs == 1);

% Correct for NaNs or Infs
log_PChat(~isfinite(log_PChat)) = log(MIN_P);

% Sum negative log likelihood
nLL = -nansum(log_PChat);

if nargout > 1
    output.resp_model = PChatL;    
end

end
function [nLL,output] = psycho_nll(params,data,mu,sigma)
%PSYCHO_NLL Negative log likelihood of psychometric functions.

MIN_P = 1e-4;   % Minimum lapse/error probability

p_true_vec = unique(data.p_true);

PChatL = zeros(size(data.resp_obs));

% If requested, fit each probability level separately
Nprobs = numel(params.psycho_mu);
for i_prob = 1:Nprobs
    psycho_mu = params.psycho_mu(i_prob);
    psycho_sigma = params.psycho_sigma(i_prob);
    psycho_gammalo = params.psycho_gammalo(i_prob);
    psycho_gammahi = params.psycho_gammahi(i_prob);
    
    if Nprobs > 1
        % Extract trials corresponding to distinct psychometric curves
        idx_trial = data.p_true == p_true_vec(i_prob);
        cc = data.signed_contrasts(idx_trial);
        
        % Psychometric curve for choosing RIGHT
        p = psychofun(cc,[psycho_mu,psycho_sigma,psycho_gammalo,psycho_gammahi]);

        % Probability or choice LEFT is 1 - psychometric curve
        PChatL(idx_trial) = 1 - p;        
    else        
        % Psychometric curve for choosing RIGHT
        p = psychofun(data.signed_contrasts,[psycho_mu,psycho_sigma,psycho_gammalo,psycho_gammahi]);
        PChatL = 1 - p;        
    end
    
end

% Log probability of observer responses (ignore no-response trials)
log_PChat = log(PChatL).*(data.resp_obs == -1) + log(1-PChatL).*(data.resp_obs == 1);

% Correct for NaNs or Infs
log_PChat(~isfinite(log_PChat)) = log(MIN_P);

% Negative log likelihood per trial
nLL = -log_PChat(:);

if nargout > 1
    output.resp_model = PChatL;    
end

end
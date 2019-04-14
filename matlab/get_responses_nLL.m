function [nLL,PChatL] = get_responses_nLL(dhat,softmax_eta,softmax_bias,lapse_rate,lapse_bias,W,resp_obs)
%GET_NLL Compute negative log likelihood of responses

MIN_P = 1e-4;   % Minimum lapse/error probability

% Compute probability given decision variable DHAT
PCx_soft = 1./(1+exp(-softmax_eta*(dhat + softmax_bias)));

% Marginalize over noisy measurements and add lapse
PChatL = qtrapz(bsxfun(@times,PCx_soft,W),2);
PChatL = lapse_rate*lapse_bias + (1-lapse_rate)*PChatL;

% Log probability of observer responses (ignore no-response trials)
log_PChat = log(PChatL).*(resp_obs == -1) + log(1-PChatL).*(resp_obs == 1);

% Correct for NaNs or Infs
log_PChat(~isfinite(log_PChat)) = log(MIN_P);

% Sum negative log likelihood
nLL = -nansum(log_PChat);

end
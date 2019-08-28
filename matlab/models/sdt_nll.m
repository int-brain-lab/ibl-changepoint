function [nLL,PChatL] = sdt_nll(params,data,priorL)
%SDT_NLL Compute negative log likelihood with signal-detection theory model.

if all(isnan(data.resp_obs))
    % Skip computation of the log-likelihood
    nLL = NaN;
    PChatL = NaN(size(data.resp_obs));
    return;
end

% Compute log prior odds per trial
logprior_odds = log(priorL./(1-priorL));

% Interpolate response probability from precomputed table
if isfield(params,'contrast_sigma')
    scidx = data.contrasts_idx + (data.S > 0)*numel(data.contrasts_vec);
else
    scidx = data.signed_contrasts_idx;
end
PChatL = zeros(size(data.resp_obs));

method = 'pchip';
for iCnd = 1:size(params.PChatL_grid,1)
    PChatL(scidx == iCnd) = ...
        interp1(params.lp_odds_grid,params.PChatL_grid(iCnd,:), ...
        logprior_odds(scidx == iCnd),method);
end

% Ensure probability is within bounds
PChatL = min(1,max(PChatL,0));

%% Compute negative log likelihood of responses

MIN_P = 1e-4;   % Minimum lapse/error probability

lapse_rate = max(MIN_P,params.lapse_rate);     % Minimum lapse to avoid numerical trouble
lapse_bias = params.lapse_bias;

% Add lapse
PChatL = lapse_rate*lapse_bias + (1-lapse_rate)*PChatL;

% Log probability of observer responses (ignore no-response trials)
log_PChat = log(PChatL).*(data.resp_obs == -1) + log(1-PChatL).*(data.resp_obs == 1);

% Correct for NaNs or Infs
log_PChat(~isfinite(log_PChat)) = log(MIN_P);

% Negative log likelihood
nLL = -log_PChat;

end

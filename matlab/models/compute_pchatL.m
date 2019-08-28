function PChatL = compute_pchatL(Nc,Nx,Np,loglikediff_t,logprior_odds,softmax_eta,softmax_bias,W)
%COMPUTE_PCHATL Compute probability of responding Left.
%
% ================ INPUT VARIABLES ====================
% NC: number of signed contrast levels. [scalar] (integer)
% NX: number of grid points. [scalar] (integer)
% NP: number of probability points. [scalar] (integer)
% LOGLIKEDIFF_T: log likelihood of Left minus Right. [Nx,Nc] (double)
% LOGPRIOR_ODDS: log prior odds. [1,Np] (double)
% SOFTMAX_ETA: softmax inverse temperature. [scalar] (double)
% SOFTMAX_BIAS: softmax bias. [scalar] (double)
% W: grid weight vector. [1,Nx] (double)
% 
% ================ OUTPUT VARIABLES ==================
% PCHATL: probability of responding left. [Nc,Np] (double)

logprior_odds3(1,1,:) = logprior_odds;

% Decision variable (3-D table) for contrast level, measurement, and log prior odds
dhat = bsxfun(@plus, loglikediff_t', logprior_odds3);

% Compute probability given decision variable DHAT
PCx_soft = 1./(1+exp(-softmax_eta*(dhat + softmax_bias)));
PCx_soft(dhat == -softmax_bias) = 0.5;

% Marginalize over noisy measurements
PChatL(:,:) = qtrapz(bsxfun(@times,PCx_soft,W),2);

end

%COMPUTE_PCHATL_TEST
%   Test script for MEX-file compute_pchatL.
%
%   Template MATLAB code generated on 26-Aug-2019 with MEXXER v0.2 
%   (https://github.com/lacerbi/mexxer).

TolErr = sqrt(eps);	% Maximum error tolerance per array element

% Define array sizes (choose reasonable values)
Nc = 100;
Np = 150;
Nx = 200;

% Randomly initialize input variables
% (or write here alternative initializations)
loglikediff_t = 10*rand([Nx,Nc]);	%LOGLIKEDIFF_T: log likelihood of Left minus Right.
logprior_odds = 10*rand([1,Np]);	%LOGPRIOR_ODDS: log prior odds.
if rand() < 0.5
    softmax_eta = 10*rand();	%SOFTMAX_ETA: softmax inverse temperature.
else
    softmax_eta = Inf;
end
softmax_bias = 10*rand([1,1]);	%SOFTMAX_BIAS: softmax bias.
W = 10*rand([1,Nx]);	%W: grid weight vector.

fprintf('=======================\n');
fprintf('Testing compute_pchatL:\n');
fprintf('=======================\n');

% Call MATLAB and MEX functions
tic; [PChatL] = compute_pchatL(Nc,Nx,Np,loglikediff_t,logprior_odds,softmax_eta,softmax_bias,W); t = toc;
tic; [PChatL_mex] = compute_pchatL_mex(Nc,Nx,Np,loglikediff_t,logprior_odds,softmax_eta,softmax_bias,W); t_mex = toc;

% Correctness check
PChatL_err = sum(abs(PChatL(:)-PChatL_mex(:)));
fprintf('Total error (PChatL): %g\n', PChatL_err);
if PChatL_err > TolErr*numel(PChatL);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in PChatL.');
end

% Runtime analysis
fprintf('Time for MATLAB code: %.3f s\n', t);
fprintf('Time for MEX file: %.3f s\n', t_mex);
fprintf('Speed gain: %.2f\n', t/t_mex);

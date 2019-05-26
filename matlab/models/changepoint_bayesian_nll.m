function [nLL,output] = changepoint_bayesian_nll(params,data,pflag)
%CHANGEPOINT_BAYESIAN_NLL Bayesian online changepoint detection observer.
% (Documentation to be written.)
%
% Author:   Luigi Acerbi
% Email:    luigi.acerbi@gmail.com
% Date:     Apr/4/2019

% If PFLAG is TRUE only compute P_ESTIMATE
if nargin < 3 || isempty(pflag); pflag = false; end

sessions = unique(data.tab(:,2));

% Pre-compute response probability as a function of signed contrast level 
% and log prior odds for speed
if ~pflag && (~isfield(params,'PChatL_grid') || isempty(params.PChatL_grid))
    np = 501;
    pgrid = linspace(params.prob_low-sqrt(eps),params.prob_high+sqrt(eps),np);

    [params.PChatL_grid,params.lp_odds_grid] = ...
        precompute_sdt(params,data,pgrid);
end

%% Split multiple sessions
if numel(sessions) > 1
    nLL = zeros(1,numel(sessions));
    output = [];
    for iSession = 1:numel(sessions)
        idx_session = data.tab(:,2) == sessions(iSession);
        data1_tab = data.tab(idx_session,:);
        data1 = format_data(data1_tab,data.filename,[data.fullname '_session' num2str(sessions(iSession))]);
        if nargout > 1
            [nLL(iSession),output1] = changepoint_bayesian_nll(params,data1,pflag);
            if isempty(output)
                output = output1;
            else
                output.p_estimate = [output.p_estimate; output1.p_estimate];
                output.runlength_post = [output.runlength_post; output1.runlength_post];
                output.fullpost = [output.fullpost; output1.fullpost];
                output.resp_model = [output.resp_model; output1.resp_model];
            end
        else
            nLL(iSession) = changepoint_bayesian_nll(params,data1,pflag);
        end
    end
    nLL = sum(nLL);
    return;
end

%% Assign observer model parameters

p_vec = params.p_vec;
Nprobs = numel(p_vec);          % # states

runlength_min = params.runlength_min;
runlength_max = params.runlength_max;
runlength_prior = str2func(params.runlength_prior);

%% Initialize and run change-point algorithm

% Construct hazard function from change-point prior (prior over run lengths)
rlprior = runlength_prior(floor(runlength_min):ceil(runlength_max));

% Trick to allow non-integer run length min and max
frac_min = 1 - (runlength_min - floor(runlength_min));
rlprior(1) = rlprior(1)*frac_min;
frac_max = 1 - (ceil(runlength_max) - runlength_max);
rlprior(end) = rlprior(end)*frac_max;
H = rlprior(:)./flipud(cumsum(flipud(rlprior(:))));

% Posterior over run lengths (from 0 to RUNLENGTH_MAX-1)
post = zeros(ceil(runlength_max),Nprobs);
post(1,:) = params.p0;  % Change in the first trial

% Transition matrix
Tmat(1,:,:) = params.Tmat;
p_vec3(1,1,:) = p_vec;

% Store full posterior at each trial?
save_fullpost = isfield(params,'save_fullpost') && params.save_fullpost ...
    && nargout > 1;

% Iterative Bayesian update over trials
[post,P,runlength_post,fullpost] = ...
    bayesianOCPDloop(data.C,post,Tmat,H,p_vec3,save_fullpost,pflag);    

%% Compute log likelihood and response probability

% Compute predictive posterior over Left
priorL = sum(bsxfun(@times,P,p_vec(:)'),2);
priorL = priorL(1:end-1); % Cut prediction for trial + 1

% Compute negative log likelihood and probability of responding L
if ~pflag
    [nLL,PChatL] = sdt_nll(params,data,priorL);
else
    nLL = NaN;  PChatL = [];
end

if nargout > 1
    output.p_estimate = priorL;
    output.runlength_post = runlength_post;
    output.resp_model = PChatL;
    output.fullpost = fullpost;
end

end

%--------------------------------------------------------------------------
function [post,P,runlength_post,fullpost] = bayesianOCPDloop(C,post,Tmat,H,p_vec,save_fullpost,pflag)
%BAYESIANCPDUPDATE Bayesian online changepoint detection update

NumTrials = size(C,1);
Nprobs = numel(p_vec);

% Table of binomial count/probability
Psi = ones(size(post,1),1,Nprobs);
Psi(1,1,:) = 1/Nprobs;

P = zeros(NumTrials+1,Nprobs);      % Posterior over state
P(1,:) = post(1,:);

runlength_post = [];    fullpost = [];
if ~pflag; runlength_post = zeros(NumTrials,size(post,1)); end

% Store full posterior at each iteration
if save_fullpost; fullpost = zeros(NumTrials,size(post,1),Nprobs); end

for t = 1:NumTrials        
    %t
    % [post,Psi,pi_post] = bayesianOCPDupdate(data.C(t),post,Psi,Tmat,H,p_vec3);


    %if mod(t,100) == 0
    %    t
    %end
    L = size(H,1);  % Width of prior over run lengths
    % Slice with nonzero hazard function (probability of change)
    idxrange = (size(post,1)-L+1:size(post,1));

    % Posterior over pi_i before observing C_t
    pi_post = bsxfun(@times, Psi, Tmat);
    predCatL(:,:) = sum(bsxfun(@times, pi_post, p_vec),3);
    predCatR(:,:) = sum(bsxfun(@times, pi_post, 1 - p_vec),3);
        
    %----------------------------------------------------------------------
    % 4a. Evaluate predictive probability
    if C(t) == 1
        predCat = predCatL ./ (predCatL + predCatR);
    else
        predCat = predCatR ./ (predCatL + predCatR);         
    end
    
    % 4b. Multiply posterior by predictive probability
    post = bsxfun(@times, post, predCat);
        
    % 4c. Evaluate posterior probability over state (only for relevant range)
    if C(t) == 1
        pi_postC = bsxfun(@times, pi_post(idxrange,:,:), p_vec);
    else
        pi_postC = bsxfun(@times, pi_post(idxrange,:,:), 1-p_vec);
    end
    pi_postC = bsxfun(@rdivide, pi_postC, sum(pi_postC,3));

    %----------------------------------------------------------------------    
    % 5a. Calculate unnormalized changepoint probabilities
    slice = post(idxrange,:);   % Slice out trials where change is possible
    currenttrial = sum( ...
        bsxfun(@times, sum(bsxfun(@times, slice, pi_postC),2), H), ...
        1);
    
    % 5b. Calculate unnormalized growth probabilities    
    post(idxrange,:) = bsxfun(@times, slice, 1-H); 
    
    % Shift posterior by one step
    post = circshift(post,1,1);
    post(1,:) = currenttrial;
    post(isnan(post)) = 0;
    
    % 5c. Calculate normalization
    Z = sum(post(:));
    
    % 5d. Normalize posterior
    post = post ./ Z;
    
    %----------------------------------------------------------------------
    % 6a. Update sufficient statistics
    Psi = circshift(Psi,1,1);
    if C(t) == 1
        Psi = bsxfun(@times, Psi, p_vec);
    else
        Psi = bsxfun(@times, Psi, 1 - p_vec);
    end
    
    Psi(1,1,:) = 1;
    Psi = bsxfun(@rdivide, Psi, sum(Psi,3));
    
    % 6b. Store predictive posterior over pi_t
    pi_post = sum(sum(bsxfun(@times,bsxfun(@times, Psi, Tmat), post),2),1);
    pi_post = pi_post / sum(pi_post);

    % 7. More bookkeeping
    
    tt = sum(post,2);

    % The predictive posterior is about the next trial
    P(t+1,:) = pi_post;
    if ~pflag; runlength_post(t,:) = tt/sum(tt); end

    if save_fullpost; fullpost(t,:,:) = post; end
end
    
    
end

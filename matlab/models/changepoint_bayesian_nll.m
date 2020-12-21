function [nLL,output] = changepoint_bayesian_nll(params,data,pflag,debug_flag)
%CHANGEPOINT_BAYESIAN_NLL Bayesian online changepoint detection observer.
% (Documentation to be written.)
%
% This is the corrected version (thanks to Charles Findling for finding a
% bug with the previous implementation, and proposing a solution).
%
% Author:   Luigi Acerbi
% Email:    luigi.acerbi@gmail.com
% Date:     Dec/20/2020

% If PFLAG is TRUE only compute P_ESTIMATE
if nargin < 3 || isempty(pflag); pflag = false; end

% Debug the algorithm using particle filtering
if nargin < 4 || isempty(debug_flag); debug_flag = false; end

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
    nLL = [];
    output = [];
    for iSession = 1:numel(sessions)
        idx_session = data.tab(:,2) == sessions(iSession);
        data1_tab = data.tab(idx_session,:);
        data1 = format_data(data1_tab,data.filename,[data.fullname '_session' num2str(sessions(iSession))]);
        
        if nargout > 1
            [nLL_temp,output1] = changepoint_bayesian_nll(params,data1,pflag);
            if isempty(output)
                output = output1;
            else
                output.p_estimate = [output.p_estimate; output1.p_estimate];
                output.runlength_post = [output.runlength_post; output1.runlength_post];
                output.fullpost = [output.fullpost; output1.fullpost];
                output.resp_model = [output.resp_model; output1.resp_model];
            end
        else
            nLL_temp = changepoint_bayesian_nll(params,data1,pflag);
        end
        nLL = [nLL; nLL_temp];
    end
    return;
end

%% Assign observer model parameters

p_vec = params.p_vec;
Nprobs = numel(p_vec);          % # states

% If 50/50 blocks are used, the mouse has knowledge of that
p0 = params.p0;
if ischar(p0) && strcmpi(p0,'ideal')
    if data.p_true(1) == 0.5
        p0 = [0 0 1];
    else
        p0 = [0.5 0.5 0];
    end
end

runlength_min = params.runlength_min;
runlength_max = params.runlength_max;
runlength_prior = str2func(params.runlength_prior);

%% Initialize and run change-point algorithm

% Construct hazard function from change-point prior (prior over run lengths)
rlprior = runlength_prior(1:ceil(runlength_max));
rlprior(1:floor(runlength_min)-1) = 0;

% Trick to allow non-integer run length min and max
frac_min = 1 - (runlength_min - floor(runlength_min));
rlprior(floor(frac_min)) = rlprior(floor(frac_min))*frac_min;
frac_max = 1 - (ceil(runlength_max) - runlength_max);
rlprior(end) = rlprior(end)*frac_max;

% Compute hazard function
H = rlprior(:)./flipud(cumsum(flipud(rlprior(:))));

% Posterior over block lengths (from 1 to RUNLENGTH_MAX)
post = zeros(ceil(runlength_max),Nprobs);
post(1,:) = p0;  % Change in the first trial

% Transition matrix
Tmat = params.Tmat;

% Store full posterior at each trial?
save_fullpost = isfield(params,'save_fullpost') && params.save_fullpost ...
    && nargout > 1;

% Iterative Bayesian update over trials
if debug_flag
    [~,priorL,runlength_post,fullpost] = ...
        particlefilter(data.C,post,Tmat,H,p_vec,save_fullpost,pflag);
else
    [~,priorL,runlength_post,fullpost] = ...
        bayesianOCPDloop(data.C,post,Tmat,H,p_vec,save_fullpost,pflag);    
end

%% Compute log likelihood and response probability

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
function [alpha_t,priorL,runlength_post,fullpost] = bayesianOCPDloop(C,alpha_t,Tmat,H,p_vec,save_fullpost,pflag)
%BAYESIANCPDUPDATE Bayesian online changepoint detection update

NumTrials = size(C,1);      % # trials
Nprobs = numel(p_vec);      % # states
L = size(H,1);              % Max block length

priorL = zeros(NumTrials,1);      % Predictive prior over current Left

runlength_post = [];    fullpost = [];
if ~pflag; runlength_post = zeros(NumTrials,L); end

% Store full posterior at each iteration?
if save_fullpost; fullpost = zeros(NumTrials,L,Nprobs); end

% Create block length transition matrix
p_block_growth = circshift(diag(1-H(:)),1,2);
p_block_growth(:,1) = H(:);
P_l_lp1(:,1,:) = p_block_growth;

% Create block type transition matrix
P_ltp1bt_btp1 = zeros(L,Nprobs,Nprobs);
P_ltp1bt_btp1(1,:,:) = Tmat;
eye3(1,:,:) = eye(Nprobs);
P_ltp1bt_btp1(2:end,:,:) = repmat(eye3,[L-1,1,1]);

for t = 1:NumTrials        
    if t >= 2
        temp(:,:) = sum(bsxfun(@times,h_t,P_l_lp1),1);
        alpha_t = squeeze(sum(bsxfun(@times,temp',P_ltp1bt_btp1),2));        
    end
        
    if C(t) == 1
        h_t = bsxfun(@times,alpha_t,p_vec);
    else
        h_t = bsxfun(@times,alpha_t,1-p_vec);        
    end
    
    zeta_t = alpha_t ./ sum(alpha_t(:));
    
    % Predictive posterior
    priorL(t,:) = sum(sum(bsxfun(@times,zeta_t,p_vec),1),2);
    
    if ~pflag; runlength_post(t,:) = sum(zeta_t,2)'; end    
    if save_fullpost; fullpost(t,:,:) = zeta_t; end
end

end
    
%--------------------------------------------------------------------------
function [alpha_t,priorL,runlength_post,fullpost] = particlefilter(C,alpha_t,Tmat,H,p_vec,save_fullpost,pflag)

NumTrials = size(C,1);      % # trials
    
priorL = zeros(NumTrials,1);      % Predictive prior over current Left

runlength_post = []; fullpost = [];
Ns = 5e3;

% Sample initial state
idx = catrnd(alpha_t(1,:),Ns);

% Initial particles and weights
X = [ones(Ns,1), idx(:)];

for t = 1:NumTrials
    
    change_idx = rand(Ns,1) <= H(X(:,1));
    X(change_idx,1) = 1;
    X(~change_idx,1) = X(~change_idx,1) + 1;
      
    idx = find(change_idx);
    for i = 1:numel(idx)
        ii = catrnd(Tmat(X(idx(i),2),:),1);
        X(idx(i),2) = ii;        
    end
    
    if C(t) == 1
        W = p_vec(X(:,2));
    else
        W = 1 - p_vec(X(:,2));
    end
    
    priorL(t) = sum(p_vec(X(:,2)))/Ns;
    
    % Resample    
    idx = catrnd(W,Ns);
    
    X = X(idx,:);
    
end


end

%--------------------------------------------------------------------------
function x = catrnd(p,n)
%CATRND Sample from categorical distribution.

maxel = 1e6;
Nel = n*numel(p);
stride = ceil(maxel/numel(p));

cdf(1,:) = cumsum(p);
u = rand(n,1)*cdf(end);

% Split for memory reasons
if Nel <= maxel
    x = sum(bsxfun(@lt, cdf, u),2) + 1;
else
    x = zeros(n,1);
    idx_min = 1;
    while idx_min <= n
        idx_max = min(idx_min+stride-1,n);
        idx = idx_min:idx_max;
        x(idx) = sum(bsxfun(@lt, cdf, u(idx)),2) + 1;
        idx_min = idx_max+1;
    end
end

end

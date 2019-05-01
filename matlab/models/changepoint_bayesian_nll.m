function [nLL,output] = changepoint_bayesian_nll(params,data,mu,sigma)
%CHANGEPOINT_BAYESIAN_NLL Bayesian online changepoint detection observer.
% (Documentation to be written.)
%
% Author:   Luigi Acerbi
% Email:    luigi.acerbi@gmail.com
% Date:     Apr/4/2019

% Potentially different SIGMA for left and right stimuli, copy if only one
% is passed
if size(sigma,2) == 1; sigma = repmat(sigma,[1,2]); end

sessions = unique(data.tab(:,2));

% Split multiple sessions
if numel(sessions) > 1
    nLL = zeros(1,numel(sessions));
    output = [];
    for iSession = 1:numel(sessions)
        idx_session = data.tab(:,2) == sessions(iSession);
        data1_tab = data.tab(idx_session,:);
        data1 = format_data(data1_tab,data.filename,[data.fullname '_session' num2str(sessions(iSession))]);
        if nargout > 1
            [nLL(iSession),output1] = changepoint_bayesian_nll(params,data1,mu(idx_session,:),sigma(idx_session,:));
            if isempty(output)
                output = output1;
            else
                output.p_estimate = [output.p_estimate; output1.p_estimate];
                output.rmse = [output.rmse; output1.rmse];
                output.post = output1.post;
                output.resp_model = [output.resp_model; output1.resp_model];
            end
        else
            nLL(iSession) = changepoint_bayesian_nll(params,data1,mu(idx_session,:),sigma(idx_session,:));            
        end
    end
    nLL = sum(nLL);
    if nargout > 1
        output.rmse = sqrt(mean(output.rmse.^2));
    end
    return;
end

NumTrials = size(data.C,1);

%% Assign observer model parameters

p_vec = params.p_vec;
Nprobs = numel(p_vec);          % # states
beta_hyp = params.beta_hyp;
if isscalar(beta_hyp); beta_hyp = beta_hyp*[1,1]; end

runlength_min = params.runlength_min;
runlength_max = params.runlength_max;
runlength_prior = str2func(params.runlength_prior);


%% Initialize inference

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

% Table of binomial count/probability
Psi = zeros(size(post,1),1,Nprobs);
Psi(1,1,:) = 1/Nprobs;

% Transition matrix
Tmat(1,:,:) = params.Tmat;
p_vec3(1,1,:) = p_vec;

%% Begin loop over trials
P = zeros(NumTrials+1,Nprobs);      % Posterior over state
P(1,:) = ones(1,Nprobs)/Nprobs;
last = zeros(NumTrials,size(post,1));

for t = 1:NumTrials
    %t
    [post,Psi,pi_post] = bayesianOCPDupdate(data.C(t),post,Psi,Tmat,H,p_vec3,beta_hyp);
    tt = nansum(post,2);

    % The predictive posterior is about the next trial
    P(t+1,:) = pi_post;
    last(t,:) = tt/sum(tt);    
end

%% Compute log likelihood and response probability

% Compute predictive posterior over Left
priorL = sum(bsxfun(@times,P,p_vec(:)'),2);
priorL = priorL(1:end-1); % Cut prediction for trial + 1

% Compute negative log likelihood and probability of responding L
[nLL,PChatL] = sdt_nll(data.S,mu,sigma,priorL,data.resp_obs,params);

if nargout > 1
    % RMSE between predictive posterior probability and true category probability
    meanP = sum(bsxfun(@times,p_vec,P(1:NumTrials,:)),2);
    output.p_estimate = meanP;
    output.rmse = sqrt(mean((meanP - data.p_true).^2));
    output.post = post;
    output.resp_model = PChatL(1:NumTrials);
end

end

%--------------------------------------------------------------------------
function [post,Psi,pi_post] = bayesianOCPDupdate(C,post,Psi,Tmat,H,p_vec,beta_hyp)
%BAYESIANCPDUPDATE Bayesian online changepoint detection update

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
    if C == 1
        predCat = predCatL ./ (predCatL + predCatR);
    else
        predCat = predCatR ./ (predCatL + predCatR);         
    end
    
    % 4b. Multiply posterior by predictive probability
    post = bsxfun(@times, post, predCat);
        
    % 4c. Evaluate posterior probability over state (only for relevant range)
    if C == 1
        pi_postC = bsxfun(@times, pi_post(idxrange,:,:), p_vec);
    else
        pi_postC = bsxfun(@times, pi_post(idxrange,:,:), 1-p_vec);
    end
    pi_postC = bsxfun(@rdivide, pi_postC, sum(pi_postC,3));

    %----------------------------------------------------------------------    
    % 5a. Calculate unnormalized changepoint probabilities
    slice = post(idxrange,:);   % Slice out trials where change is possible
    currenttrial = nansum( ...
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
    if C == 1
        Psi = bsxfun(@times, Psi, p_vec);
    else
        Psi = bsxfun(@times, Psi, 1 - p_vec);
    end
    
    % Hyperprior
    Psi(1,1,:) = exp(beta_hyp(1).*log(p_vec) + beta_hyp(2).*log(1-p_vec));
    Psi = bsxfun(@rdivide, Psi, sum(Psi,3));
    
    % 6b. Store predictive posterior over pi_t
    pi_post = nansum(sum(bsxfun(@times,bsxfun(@times, Psi, Tmat), post),2),1);
    pi_post = pi_post / sum(pi_post);

end

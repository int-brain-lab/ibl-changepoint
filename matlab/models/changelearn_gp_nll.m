function [nLL,output] = changelearn_gp_nll(params,data,tmax)
%CHANGELEARN_GP_NLL Bayesian online changepoint detection learner observer.
% (Documentation to be written.)
%
% Author:   Luigi Acerbi
% Email:    luigi.acerbi@gmail.com
% Date:     May/22/2019

if nargin < 3 || isempty(tmax); tmax = Inf; end


%% Trim dataset if requested

if isfinite(tmax)
    fprintf('Trimming data down to %d trials.\n',tmax); 
    data.tab = data.tab(1:min(end,tmax),:);
    data = format_data(data.tab,data.filename,data.fullname);    
end

%% Initialize change-point learner

% Get parameter bounds
bounds = setup_params([],params);
idx = 5:8;
bounds.LB = bounds.LB(idx);
bounds.UB = bounds.UB(idx);
bounds.PLB = bounds.PLB(idx);
bounds.PUB = bounds.PUB(idx);
bounds.x0 = bounds.x0(idx);

Ntheta = size(param_grid,1);
p_grid = NaN(size(data.tab,1),Ntheta);

%% Loop over change-point parameters
for iIter = 1:Ntheta
    if mod(iIter,10) == 0; fprintf('Iteration %d / %d\n',iIter,Ntheta); end
    
    % Create observer with specific change-point parameters
    theta = [params.theta,param_grid(iIter,:)];
    params1 = setup_params(theta,params.changeparams);
    
    % Run observer and get predictions for each trial
    [~,output1] = changepoint_bayesian_nll(params1,data,1);
    p_grid(:,iIter) = output1.p_estimate;    
end

%% Compute posterior over change-point parameters

loglike = zeros(size(p_grid));
loglike(data.C == 1,:) = log(p_grid(data.C == 1,:));
loglike(data.C ~= 1,:) = log(1-p_grid(data.C ~= 1,:));

% Flat prior for the moment
loglike_grid = cumsum(loglike,1);
post_grid = exp(bsxfun(@minus,loglike_grid,max(loglike_grid,[],2)));
post_grid = bsxfun(@rdivide,post_grid,sum(post_grid,2));

% Compute predictive posterior over Left
priorL = sum(bsxfun(@times,post_grid,p_grid),2);

%% Compute log likelihood and response probability

% Pre-compute response probability as a function of signed contrast level 
% and log prior odds for speed

np = 501;
pgrid = linspace(prob_low(1)-sqrt(eps),prob_high(end)+sqrt(eps),np);

[params.PChatL_grid,params.lp_odds_grid] = ...
    precompute_sdt(params,data,pgrid);

% Compute negative log likelihood and probability of responding L
[nLL,PChatL] = sdt_nll(params,data,priorL);

if nargout > 1
    % Time steps
    tt = unique(round(exp(linspace(log(1),log(1e5),1509))));
    tt = tt(tt <= size(priorL,1));
    
    output.p_estimate = priorL;
    output.resp_model = PChatL;    
    output.param_grid = param_grid;
    output.trials = tt;
    output.loglike_grid_trials = loglike_grid(tt,:);
end

end
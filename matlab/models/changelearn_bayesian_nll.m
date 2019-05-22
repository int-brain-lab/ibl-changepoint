function [nLL,output] = changelearn_bayesian_nll(params,data)
%CHANGELEARN_NLL Bayesian online changepoint detection learner observer.
% (Documentation to be written.)
%
% Author:   Luigi Acerbi
% Email:    luigi.acerbi@gmail.com
% Date:     May/22/2019

%% Initialize change-point learner

% Create parameter grid
if isfield(params,'Ngrid') && ~isempty(params.Ngrid)
    Ngrid = params.Ngrid;
else
    Ngrid = 10;
end
runlength_tau = log(round(exp(linspace(log(1),log(100),Ngrid))));
runlength_min = log(round(exp(linspace(log(1),log(90),Ngrid))));
prob_low = linspace(0.1,0.9,Ngrid-1);
prob_high = linspace(0.1,0.9,Ngrid-1);

param_grid = combvec(runlength_tau,runlength_min,prob_low,prob_high)';
param_grid(param_grid(:,4)<param_grid(:,3),:) = [];

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

% Flat prior for the moment
logpost_grid = cumsum(log(p_grid),1);
post_grid = exp(bsxfun(@minus,logpost_grid,max(logpost_grid,[],2)));
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
    output.p_estimate = priorL;
    output.resp_model = PChatL;    
    output.param_grid = param_grid;
    output.p_grid = p_grid;
end

end
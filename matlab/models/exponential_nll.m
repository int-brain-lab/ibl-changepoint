function [nLL,output] = exponential_nll(params,data)
%EXPONENTIAL_NLL Exponential-averaging observer.
% (Documentation to be written.)
%
% Author:   Luigi Acerbi
% Email:    luigi.acerbi@gmail.com
% Date:     May/6/2019

sessions = unique(data.tab(:,2));

PMIN = 0.01;
PMAX = 0.99;

% Pre-compute response probability as a function of signed contrast level 
% and log prior odds for speed
if ~isfield(params,'PChatL_grid') || isempty(params.PChatL_grid)
    np = 501;
    pgrid = linspace(PMIN-sqrt(eps),PMAX+sqrt(eps),np);

    [params.PChatL_grid,params.lp_odds_grid] = ...
        precompute_sdt(params,data,pgrid);
end

% Split multiple sessions
if numel(sessions) > 1
    nLL = zeros(1,numel(sessions));
    output = [];
    for iSession = 1:numel(sessions)
        idx_session = data.tab(:,2) == sessions(iSession);
        data1_tab = data.tab(idx_session,:);
        data1 = format_data(data1_tab,data.filename,[data.fullname '_session' num2str(sessions(iSession))]);
        if nargout > 1
            [nLL(iSession),output1] = exponential_nll(params,data1);
            if isempty(output)
                output = output1;
            else
                output.p_estimate = [output.p_estimate; output1.p_estimate];
                output.resp_model = [output.resp_model; output1.resp_model];
            end
        else
            nLL(iSession) = exponential_nll(params,data1);            
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

beta_hyp = params.beta_hyp;
if isscalar(beta_hyp); beta_hyp = beta_hyp*[1,1]; end
tau = params.runlength_tau;

%% Initialize inference

windowSize = 100;
Lcounts = data.C == 1;
Rcounts = data.C ~= 1;

ff = exp(-(0:windowSize)/tau);   % Exponential filter

Lexp_avg = filter(ff,1,Lcounts);
Rexp_avg = filter(ff,1,Rcounts);

priorcountsL = beta_hyp(1);
priorcountsR = beta_hyp(2);

% Posterior mean
priorL = (Lexp_avg + priorcountsL) ./ ...
    (Lexp_avg + priorcountsL + Rexp_avg + priorcountsR);

priorL = min(max(priorL,PMIN),PMAX);

%% Compute log likelihood and response probability

% Compute negative log likelihood and probability of responding L
[nLL,PChatL] = sdt_nll(params,data,priorL);

if nargout > 1
    output.p_estimate = priorL;
    output.rmse = sqrt(mean((priorL - data.p_true).^2));
    output.resp_model = PChatL(1:NumTrials);
end

end

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
    nLL = [];
    output = [];
    for iSession = 1:numel(sessions)
        idx_session = data.tab(:,2) == sessions(iSession);
        data1_tab = data.tab(idx_session,:);
        data1 = format_data(data1_tab,data.filename,[data.fullname '_session' num2str(sessions(iSession))]);
        if nargout > 1
            [nLL_temp,output1] = exponential_nll(params,data1);
            if isempty(output)
                output = output1;
            else
                output.p_estimate = [output.p_estimate; output1.p_estimate];
                output.resp_model = [output.resp_model; output1.resp_model];
            end
        else
            nLL_temp = exponential_nll(params,data1);            
        end
        nLL = [nLL; nLL_temp];
    end
    if nargout > 1
        output.rmse = sqrt(mean(output.rmse.^2));
    end
    return;
end

%% Observer model

NumTrials = size(data.C,1);
tau = params.runlength_tau;

windowSize = 100;

if contains(params.model_name,'linweight')
    ww = params.contrastweights;
    p = polyfit([0,1],[ww,1],1);
        
    Lcounts = zeros(size(data.C));
    idx = data.C == 1;
    Lcounts(idx) = data.contrasts(idx)*p(1) + p(2);
    
    Rcounts = zeros(size(data.C));
    idx = data.C ~= 1;
    Rcounts(idx) = data.contrasts(idx)*p(1) + p(2);    
else
    Lcounts = data.C == 1;
    Rcounts = data.C ~= 1;
end

if isfield(params,'tau_quad')
    [data.trial_block,data.trials2end_block] = ...
        get_trial_number_within_blocks(data.tab);
    
    % Compute trial-dependent TAU
    tau_trial = tau + params.tau_quad.*(log(data.trial_block+1) - log(params.tau_quadmu));
    tau_trial = max(tau_trial,1);
    
    tau_vec = linspace(min(tau_trial), max(tau_trial),100);
    
    for iTau = 1:numel(tau_vec)
        ff = exp(-(0:windowSize)/tau_vec(iTau));   % Exponential filter        
        Lexp_avg_mat(:,iTau) = filter(ff,1,Lcounts);
        Rexp_avg_mat(:,iTau) = filter(ff,1,Rcounts);         
    end
    
    [~,idx_tau] = min(abs(bsxfun(@minus,tau_trial,tau_vec)),[],2);
    
    idx_vv = sub2ind([numel(tau_trial),numel(tau_vec)],(1:numel(tau_trial))',idx_tau);
    Lexp_avg = Lexp_avg_mat(idx_vv);
    Rexp_avg = Rexp_avg_mat(idx_vv);
else
    ff = exp(-(0:windowSize)/tau);   % Exponential filter
    Lexp_avg = filter(ff,1,Lcounts);
    Rexp_avg = filter(ff,1,Rcounts);
end


if isfield(params,'lnp_hyp')
    MIN_P = 1e-6;
    np = 201;
    nhyp = numel(params.lnp_hyp);
    lnp_hyp = params.lnp_hyp;
    p = [lnp_hyp(1:nhyp/2),0,lnp_hyp(nhyp/2+1:end)];
    %p = [log(MIN_P) lnp_hyp 0 fliplr(lnp_hyp) log(MIN_P)];
    px = min(max(MIN_P,linspace(0,1,numel(p))),1-MIN_P);    
    xx = linspace(MIN_P,1-MIN_P,np);
    lnpvec = interp1(px,p,xx,'pchip');
    
    likeL = bsxfun(@times,Lexp_avg,log(xx));
    likeR = bsxfun(@times,Rexp_avg,log(1-xx));
    
    lnpost = bsxfun(@plus,likeL + likeR,lnpvec);
    
    % Normalize posterior
    nZ = max(lnpost,[],2);
    postL = exp(bsxfun(@minus,lnpost,nZ));
    postL = bsxfun(@rdivide,postL,sum(postL,2));    
    
    % Mean of the posterior
    priorL = sum(bsxfun(@times,postL,xx),2);
    
elseif isfield(params,'beta_w')
    % Mixture of Beta distributions prior
    
    beta_w = params.beta_w;
    beta_hyp = params.beta_hyp;
    if isscalar(beta_hyp); beta_hyp = beta_hyp*[1,1]; end
       
    MIN_P = 1e-6;
    np = 201;
    xx = linspace(MIN_P,1-MIN_P,np);
    pvec = beta_w*betapdf(xx,beta_hyp(1),beta_hyp(2)) + ...
        (1-beta_w)*betapdf(xx,beta_hyp(3),beta_hyp(4));
    lnpvec = log(pvec);
    
    likeL = bsxfun(@times,Lexp_avg,log(xx));
    likeR = bsxfun(@times,Rexp_avg,log(1-xx));    
    lnpost = bsxfun(@plus,likeL + likeR,lnpvec);
    
    % Normalize posterior
    nZ = max(lnpost,[],2);
    postL = exp(bsxfun(@minus,lnpost,nZ));
    postL = bsxfun(@rdivide,postL,sum(postL,2));    
    
    % Mean of the posterior
    priorL = sum(bsxfun(@times,postL,xx),2);        
    
else
    beta_hyp = params.beta_hyp;
    if isscalar(beta_hyp); beta_hyp = beta_hyp*[1,1]; end

    priorcountsL = beta_hyp(1);
    priorcountsR = beta_hyp(2);
    
    % Posterior mean
    priorL = (Lexp_avg + priorcountsL) ./ ...
        (Lexp_avg + priorcountsL + Rexp_avg + priorcountsR);
    
end


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

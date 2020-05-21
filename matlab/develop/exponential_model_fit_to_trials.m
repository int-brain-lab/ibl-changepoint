function [tau,beta,tau0,beta0] = exponential_model_fit_to_trials(data)
%EXPONENTIAL_MODEL_FIT_TO_TRIALS Fit exponential model to specific trials


if ~iscell(data); data = {data}; end


for iData = 1:numel(data)
    try
        if ischar(data{iData})
            dataset = read_data_from_csv(data{iData});
        else
            dataset = data{iData};
        end
        [tau_hat,beta_hat] = exponential_model_fit(dataset);
    catch
        tau_hat = NaN(1,size(tau,2)); beta_hat = NaN(1,size(beta,2));
    end
    tau(iData,:) = tau_hat(1:end-1);
    beta(iData,:) = beta_hat(1:end-1);
    
    tau0(iData) = tau_hat(end);
    beta0(iData) = beta_hat(end);
end

% Make plots
subplot(1,2,1);
curveplot(tau,'Estimated \tau');
xlim([0,60]);
ylim([1,30]);
set(gca,'XTick',0:20:100,'YTick',5:5:100);

subplot(1,2,2);
curveplot(beta,'Estimated \beta');
xlim([0,60]);
ylim([0,2]);
set(gca,'XTick',0:20:100,'YTick',0:0.5:2);


end

function h = curveplot(Y,ytext)

fontsize = 18;
axesfontsize = 14;

N = size(Y,2);

xx = (-20:N-21);
yy = nanmean(Y,1);
yy_serr = nanstd(Y,[],1);

yyerr_down = yy - 1.96*yy_serr;
yyerr_up = yy + 1.96*yy_serr;
xxerr = [xx, fliplr(xx)];
yyerr = [yyerr_down, fliplr(yyerr_up)];
    
fill(xxerr, yyerr,'k','FaceAlpha',0.5,'LineStyle','none'); hold on;
h(1) = plot(xx,yy,'k-','LineWidth',2);

set(gcf,'Color','w');
set(gca,'TickDir','out','FontSize',axesfontsize);
box off;
xlabel('Trials from block start','FontSize',fontsize);
ylabel(ytext,'FontSize',fontsize);


end



function [tau_hat,beta_hat] = exponential_model_fit(data)

lambda = 1e-3;

params = params_new('changepoint_runlength',data);
[~,output] = changepoint_bayesian_nll(params,data,1);
p_model = output.p_estimate;



tau_grid = exp(linspace(log(1),log(60),70));
Lcounts = data.C == 1;
Rcounts = data.C ~= 1;

for i = 1:numel(tau_grid)
    windowSize = 100;
    
    tau = tau_grid(i);
    ff = exp(-(0:windowSize)/tau);   % Exponential filter
    
    Lexp_avg(:,i) = filter(ff,1,Lcounts);
    Rexp_avg(:,i) = filter(ff,1,Rcounts);
end

betahyp_grid = linspace(sqrt(0),sqrt(2),13).^2;

for i = 1:numel(betahyp_grid)    
    priorcountsL = betahyp_grid(i);
    priorcountsR = betahyp_grid(i);    
    priorL(:,:,i) = (Lexp_avg + priorcountsL) ./ ...
        (Lexp_avg + priorcountsL + Rexp_avg + priorcountsR);    
end

X = data.tab;
X(:,4) = X(:,4).*sign(X(:,5));      % Convert contrast to signed contrast
sessions = unique(X(:,2))';

X = [X, zeros(size(X,1),2)];        % Add two empty columns

for iSession = sessions
    idx_session = (X(:,2) == iSession);    
    [trial_block,trials2end_block] = ...
        get_trial_number_within_blocks(X(idx_session,:));
    X(idx_session,end-1) = trial_block;
    X(idx_session,end) = trials2end_block;
end

biased_block = X(:,3) == 0.2 | X(:,3) == 0.8;

q = min(1-lambda,max(lambda,priorL));

for iiTrial = 0:60    
    idx = (X(:,end-1) == iiTrial) & biased_block;        
    KL = squeeze(sum(- p_model(idx).*log(q(idx,:,:)) ...
        - (1 - p_model(idx)).*(log(1-q(idx,:,:))),1));
        
    [~,idx_best] = min(KL(:));
    
    idx1 = mod(idx_best-1,numel(tau_grid))+1;
    idx2 = floor((idx_best-1)./numel(tau_grid))+1;
    
    tau_hat(iiTrial+21) = tau_grid(idx1);
    beta_hat(iiTrial+21) = betahyp_grid(idx2);    
end

for iiTrial = 1:20
    idx = (X(:,end) == iiTrial) & biased_block;        
    KL = squeeze(sum(- p_model(idx).*log(q(idx,:,:)) ...
        - (1 - p_model(idx)).*(log(1-q(idx,:,:))),1));
        
    [~,idx_best] = min(KL(:));
    
    idx1 = mod(idx_best-1,numel(tau_grid))+1;
    idx2 = floor((idx_best-1)./numel(tau_grid))+1;
    
    tau_hat(21 - iiTrial) = tau_grid(idx1);
    beta_hat(21 - iiTrial) = betahyp_grid(idx2);    
end


idx = biased_block;        
KL = squeeze(sum(- p_model(idx).*log(q(idx,:,:)) ...
    - (1 - p_model(idx)).*(log(1-q(idx,:,:))),1));

[~,idx_best] = min(KL(:));

idx1 = mod(idx_best-1,numel(tau_grid))+1;
idx2 = floor((idx_best-1)./numel(tau_grid))+1;

tau_hat(end+1) = tau_grid(idx1);
beta_hat(end+1) = betahyp_grid(idx2);    



end
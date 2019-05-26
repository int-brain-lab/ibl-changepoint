function explearn_test()
%EXPLEARN_TEST Test speed of learning of truncated exponential.

fontsize = 18;
axesfontsize = 14;

Ng = 201;
Nreps = 24;
Nblocks = 220;
tmax = 1e4;

tau_grid = exp(linspace(log(1),log(100),Ng));
loglike = zeros(Nreps,Ng);

tau0 = 60;
min_rl = 20;
max_rl = 100;


x = Inf(Nreps*Nblocks,1);
idx = true(size(x));

% Get samples from truncated exponential
while any(idx)
    x(idx) = min_rl + exprnd(tau0,[sum(idx),1]);
    idx = x > max_rl;
end
x = reshape(x,[Nreps,Nblocks]);
x = x + 3*randn(size(x));

% Compute normalization constant of truncated exponential
lnZ = zeros(1,Ng);
xx = linspace(0,max_rl-min_rl,1e4);
dxx = xx(2)-xx(1);
for i = 1:Ng
    temp = exp(-xx./tau_grid(i));
    lnZ(i) = log(qtrapz(temp)*dxx);
end

% Compute posterior over tau for each repetition
mu = zeros(Nreps,Nblocks);

post = exp(loglike - max(loglike,[],2));
post = post./sum(post,2);    
mu0 = exp(sum(post.*log(tau_grid),2));

for iBlock = 1:Nblocks
    loglike = loglike + bsxfun(@minus,-bsxfun(@rdivide,x(:,iBlock)-min_rl,tau_grid),lnZ);
    
    post = exp(loglike - max(loglike,[],2));
    post = post./sum(post,2);    
    mu(:,iBlock) = exp(sum(post.*log(tau_grid),2));
end

xall = 0:10:1e4;
mu_all = zeros(Nreps,numel(xall));

subplot(1,2,1);
for iRep = 1:Nreps    
    xrep = [1,cumsum(x(iRep,:))];
    murep = [mu0(iRep),mu(iRep,:)];
    plot(xrep,murep,'-','Color',0.5*[1 1 1],'LineWidth',1);
    
    mu_all(iRep,:) = interp1(xrep,murep,xall);
    hold on;
end
plot(xall,mean(mu_all,1),'-k','LineWidth',3);
xlim([0 tmax]);

box off;
set(gca,'Xscale','log','TickDir','out');
set(gcf,'Color','w');
set(gca,'XTick',[1,10,100,1e3,1e4],...
    'XTickLabel',{'1','10','100','1000','10000'},...
    'Fontsize',axesfontsize);
xlim([10,tmax]);
ylim([0,100]);
title('Inferring exponential time constant \tau','Fontsize',fontsize);
xlabel('Trial','Fontsize',fontsize);
ylabel('Posterior mean \tau','Fontsize',fontsize);

subplot(1,2,2);
for iRep = 1:Nreps
%    plot(tau_grid,loglike(iRep,:)-max(loglike(iRep,:)),'-','Color',0.5*[1 1 1],'LineWidth',1); hold on; ylim([-100 0])    
    plot(tau_grid,post(iRep,:),'-','Color',0.5*[1 1 1],'LineWidth',1); hold on;
end

box off;
set(gca,'TickDir','out');
set(gca,'Fontsize',axesfontsize);
xlim([1 100]);
%xlim([10,tmax]);
%ylim([0,100]);
title('Posterior after ~10000 trials','Fontsize',fontsize);
xlabel('\tau','Fontsize',fontsize);
ylabel('Posterior probability','Fontsize',fontsize);





end
% Scatterplot comparison

close all;
figure(1);

model_name = 'exponential_contrastnoise';

batch_compute_parameters;
Nparams = size(theta_mle,2);
pnames = params.names;
logflag = false(1,Nparams);
switch model_name
    case 'exponential_contrastnoise'; logflag([1:2,5]) = true;
    case 'changepoint_contrastnoise_runlength_probs'; logflag([1:2,5:6]) = true;    
end

for iParam = 1:Nparams
    subplot(2,4,iParam);
    
    xx = theta(:,iParam);
    yy = theta_mle(:,iParam);
    if logflag(iParam)
        xx = exp(xx);
        yy = exp(yy);
    end
    
    scatter(xx,yy); hold on;
    box off;
    axis square;
    
    name = pnames{iParam};
    name(name == '_') = '-';
    
    xlabel([name ' (posterior mean)']);
    ylabel([name ' (MLE)']);
    
    plot(xlim,xlim,'k:','LineWidth',1);
    
    if iParam == 1
        title('Comparison of maximum-likelihood estimation and posterior mean fits');
    end        
end

set(gcf,'Color','w');
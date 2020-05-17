function plot_exponential_quad(params)

xx = linspace(0,59,100);

for iParam = 1:numel(params)
    
    tau = params{iParam}.runlength_tau;
    tau_quad = params{iParam}.tau_quad;
    tau_quadmu = params{iParam}.tau_quadmu;
        
    tau_y(iParam,:) = max(1,tau + tau_quad.*(log(xx+1) - log(tau_quadmu))).^2;
end

plot(xx,tau_y,'k-');








end
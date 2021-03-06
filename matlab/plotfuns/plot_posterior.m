function plot_posterior(params)
%PLOT_POSTERIOR Plot (approximate) parameter posterior.

if ~isfield(params,'vbmc_fits') || isempty(params.vbmc_fits)
    error('No approximate posterior in PARAMS struct!');
end

xrnd = get_posterior_samples(params,1e6);

bounds = setup_params([],params);
bb = [bounds.PLB; bounds.PUB];

if isfield(params,'logspace')
    logspace = logical(params.logspace);
    xrnd(:,logspace) = exp(xrnd(:,logspace));
    bb(:,logspace) = exp(bb(:,logspace));
end


pnames = params.names;
pnames = cellfun(@(str) strrep(str,'_','-'),pnames,'UniformOutput',false);

cornerplot(xrnd,pnames,[],bb);
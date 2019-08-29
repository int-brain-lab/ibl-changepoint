function [x,fval,x0] = optimize_model(params,data,x0)
%OPTIMIZE_MODEL Fit model by optimization (maximum-likelihood or MAP).

if nargin < 3; x0 = []; end             % Starting point

bounds = setup_params([],params);       % Get parameter bounds
nvars = numel(bounds.PLB);  % Number of dimensions

% Fits all models but psychometric curves with BADS    
bads_flag = ~contains(params.model_name,{'psychofun'});

if isempty(x0); x0 = NaN(1,nvars); end
N0 = size(x0,1);
xrnd = bsxfun(@plus,bsxfun(@times,rand(N0,nvars),bounds.PUB-bounds.PLB),bounds.PLB);
for ii = 1:N0
    idx = ~isfinite(x0(ii,:));
    x0(ii,idx) = xrnd(ii,idx);
end

% Evaluate nLL on a bunch of quasirandom points, start from best
Ninit = 5*nvars;    
P = sobolset(nvars);
P = scramble(P,'MatousekAffineOwen');
xx = bsxfun(@plus,bsxfun(@times,net(P,Ninit),bounds.PUB-bounds.PLB),bounds.PLB);    
x0_list = [x0; xx];

nll0 = zeros(size(x0_list,1),1);
for i0 = 1:size(x0_list,1)
    nll0(i0) = sum(nllfun(x0_list(i0,:),params,data));
end
[~,idx0] = min(nll0);
x0 = x0_list(idx0,:);

% Run optimization from best initial point
fun = @(x_) sum(nllfun(x_,params,data));
if bads_flag
%         fminbayesopts = fminbayes('defaults');
%         fminbayes(fun,x0,bounds.LB,bounds.UB,bounds.PLB,bounds.PUB,fminbayesopts);

    MaxFunEvals = 2e3;
    badopts = bads('defaults');
    badopts.MaxFunEvals = MaxFunEvals;        
    [x,fval] = bads(fun,x0,...
        bounds.LB,bounds.UB,bounds.PLB,bounds.PUB,[],badopts);
else
    fminopts.Display = 'iter';
    [x,fval] = fmincon(fun,x0,...
        [],[],[],[],bounds.LB,bounds.UB,[],fminopts);            
end












end
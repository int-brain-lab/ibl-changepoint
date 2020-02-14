function [X,pnames,params,data_names] = plot_posterior_comparison(data_names,model,color,style,mark)
%PLOT_POSTERIOR_COMPARISON

% Example: [X,pnames,params,data_names] = plot_posterior_comparison(get_mice_list('guido'),'exponential_contrastnoise',[0 0.8 0; 0 0.8 0; 0.3 1 0.3; 0 0 0.8; 0 0 0.8; 0.3 0.3 1; 0.8 0 0; 0.8 0 0; 1 0.3 0.3],{'-',':','-','-',':','-','-',':','-'})

if nargin < 3; color = []; end
if nargin < 4; style = []; end
if nargin < 5; mark = []; end

Nsamples = 2e4;     % Posterior samples used to make plots
Nkde = 2^12;        % Points used for kernel density estimation

% Load all model fits
for i = 1:numel(data_names)
    params{i} = load_model_fit(data_names{i},model);
    X{i} = get_posterior_samples(params{i},Nsamples);
    X{i} = setup_params(X{i},params{i});
    % Get empirical parameter bounds from samples
    if i == 1
        lb = min(X{1});             ub = max(X{1});
    else
        lb = min(lb,min(X{i}));     ub = max(ub,max(X{i}));
    end
end

% Set parameter bounds for plotting
%bounds = setup_params([],params{1});
%bb = setup_params([bounds.PLB;bounds.PUB],params{1});
%lb = min(lb,bb(1,:));
%ub = max(ub,bb(2,:));

pnames = params{1}.names;
pnames = cellfun(@(str) strrep(str,'_','-'),pnames,'UniformOutput',false);

Ntheta = numel(pnames); % Number of parameters

Nrows = floor(sqrt(Ntheta));
Ncols = ceil(sqrt(Ntheta));
if Nrows*Ncols < Ntheta; Ncols = Ncols+1; end

for d = 1:Ntheta
    subplot(Nrows,Ncols,d);
    for i = 1:numel(X)
        [~,p,xmesh] = kde1d(X{i}(:,d),Nkde);
        if ~isempty(mark); mm = mark{i}; else; mm = 'none'; end
        if ~isempty(style); ss = style{i}; else; ss = '-'; end
        if ~isempty(color); col = color(i,:); else; col = [0 0 0]; end
        h(i) = plot(xmesh,p,'Color',col,'LineWidth',2,'LineStyle',ss,'Marker',mm); hold on;
    end
    xlim([lb(d),ub(d)]);
    set(gca,'TickDir','out','YTick',0);
    box off;
    xlabel(pnames{d});
    ylabel('Posterior pdf');
    if d == 1
        ltext = cellfun(@(str) strrep(str,'_','-'),data_names,'UniformOutput',false);        
        hl = legend(h,ltext{:});
        set(hl,'Box','off');
    end
end
set(gcf,'Color','w');
function modelfits = batch_model_fit(model_names,data_name,Nopts,methods_flags,refit_flags,empirical_list)
%BATCH_MODEL_FIT Fit a batch of models and save results.

% TODO: Change saving method to separate files for each data/model

if nargin < 3 || isempty(Nopts); Nopts = [10,5]; end
if numel(Nopts) > 1; Nvbmc = Nopts(2); else; Nvbmc = ceil(Nopts(1)/2); end
if nargin < 4; methods_flags = []; end
if nargin < 5 || isempty(refit_flags); refit_flags = false(1,4); end
if nargin < 6; empirical_list = []; end

close all;

% Models to be fitted
if isempty(model_names)
    model_names = get_model_list('default');
elseif ischar(model_names)
    model_names = {model_names};
end

% Plotting grid
switch numel(model_names)
    case 1; plotrows = 1; plotcols = 1;
    case 2; plotrows = 1; plotcols = 2;
    case 3; plotrows = 1; plotcols = 3;
    case 4; plotrows = 2; plotcols = 2;
    case 5; plotrows = 2; plotcols = 3;
    case 6; plotrows = 2; plotcols = 3;
    otherwise
        plotrows = 2;
        plotcols = ceil(numel(model_names)/2);
end

% Load data file
data = read_data_from_csv(data_name);

modelfits.data = data;
modelfits.params = [];

for iModel = 1:numel(model_names)
    
    fprintf('Fitting dataset: %s. Model: %s.\n\n',data.fullname,model_names{iModel});
    
    if isempty(model_names{iModel}); continue; end
    
    % Fit models
    params = fit_model(model_names{iModel},data, ...
        [Nopts(1),Nvbmc],methods_flags,refit_flags,[],1,empirical_list);       
        
    % Store fits
    modelfits.params{end+1} = params;

    % Plot
    subplot(plotrows,plotcols,iModel);
    plot_fit(data,params,params.model_desc);
end

mypath = which('savefigure.m');
savefigure([fileparts(mypath) filesep() data_name]);

end
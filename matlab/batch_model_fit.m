function modelfits = batch_model_fit(model_names,data_filename,Nopts,vbmc_flag,refit_flag)
%BATCH_MODEL_FIT Fit a batch of models and save results.

% TODO: Change saving method to separate files for each data/model

if nargin < 3 || isempty(Nopts); Nopts = 10; end
if nargin < 4 || isempty(vbmc_flag); vbmc_flag = false; end
if nargin < 5 || isempty(refit_flag); refit_flag = false; end

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
data = read_data_from_csv(data_filename);

modelfits.data = data;
modelfits.params = [];

for iModel = 1:numel(model_names)
    
    if isempty(model_names{iModel}); continue; end
    
    % Fit and save model
    params = fit_model(model_names{iModel},data,Nopts,vbmc_flag,refit_flag);        
    save_model_fit(data_filename,params);
        
    % Store fits
    modelfits.params{end+1} = params;

    % Plot
    subplot(plotrows,plotcols,iModel);
    plot_fit(data,params,params.model_desc);
end

mypath = which('savefigure.m');
savefigure([fileparts(mypath) filesep() data_filename]);

end
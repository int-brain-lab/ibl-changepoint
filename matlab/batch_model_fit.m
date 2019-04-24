function modelfits = batch_model_fit(model_names,filename,Nopts,vbmc_flag,refit_flag)
%BATCH_MODEL_FIT Fit a batch of models and save results.

if nargin < 3 || isempty(Nopts); Nopts = 10; end
if nargin < 4 || isempty(vbmc_flag); vbmc_flag = false; end
if nargin < 5 || isempty(refit_flag); refit_flag = false; end

mypath = fileparts(mfilename('fullpath'));
fits_path = [mypath filesep 'fits'];
addpath(fits_path);
addpath([mypath filesep 'utils']);

close all;

% Models to be fitted
if isempty(model_names)
    model_names{1} = 'psychofun';
    model_names{2} = 'omniscient';
    model_names{3} = 'omniscient_lapse';
    model_names{4} = 'omniscient_biasedlapse';
    % model_names{5} = 'changepoint_biasedlapse_runlength';
    % model_names{6} = 'changepoint_biasedlapse_runlength_freesym';
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
data = read_data_from_csv(filename);

modelfits.data = data;
modelfits.params = [];

% Read existing MODELFITS if present
matfilename = [filename '_fits.mat'];
if exist(matfilename,'file')
    fprintf('Found existing file ''%s'', loading previous fits.\n', matfilename);
    for iTry = 1:10
        try load(matfilename); break;
        catch; fprintf('Try %d: I/O error. Waiting a few second before retrying...\n', iTry); pause(5 + 5*rand());
        end
    end
end

for iModel = 1:numel(model_names)
    
    if isempty(model_names{iModel}); continue; end
    
    params = fit_model(model_names{iModel},data,Nopts,vbmc_flag,refit_flag);
    
    % Add fit to MODELFITS cell array
    
    % Check if model was already fitted to this session
    idx_params = [];
    for iFit = 1:numel(modelfits.params)
        pp = modelfits.params{iFit}; 
        if strcmp(pp.model_name,model_names{iModel})
            idx_params = iFit;
        end
    end
    
    % Save fits
    if ~isempty(idx_params)
        modelfits.params{idx_params} = params;
    else
        modelfits.params{end+1} = params;
    end    
    
    save([fits_path filesep() filename '_fits.mat'],'modelfits');
    
    % Plot
    subplot(plotrows,plotcols,iModel);
    plot_fit(data,params,params.model_desc);
end

mypath = which('savefigure.m');
savefigure([fileparts(mypath) filesep() filename]);

end
function batch_model_plots(model_names,filename,model_names_train,filename_train)
%BATCH_MODEL_PLOTS Plot a batch of models results, possibly with transfer.

if nargin < 3 || isempty(model_names_train); model_names_train = model_names; end
if nargin < 4 || isempty(filename_train); filename_train = filename; end

mypath = fileparts(which('changepoint_bayesian_nll'));
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
    model_names{5} = 'omniscient_altnoise';
elseif ischar(model_names)
    model_names = {model_names};
end

if isempty(model_names_train); model_names_train = model_names; end

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

% Read training MODELFITS
matfilename_train = [filename_train '_fits.mat'];
train_fits = load(matfilename_train);

% Read test MODELFITS (only for psychometric curve plotting)
matfilename_test = [filename '_fits.mat'];
try
    test_fits = load(matfilename_test);
catch
    test_fits = [];
end

for iModel = 1:numel(model_names)
    if contains(model_names{iModel},'psychofun')
        % Look for psychometric model of the TEST data
        idx_params = [];
        for iFit = 1:numel(test_fits.modelfits.params)
            pp = test_fits.modelfits.params{iFit}; 
            if strcmp(pp.model_name,model_names_train{iModel})
                idx_params = iFit;
            end
        end
        params = test_fits.modelfits.params{idx_params};
        
    else
        % Look for training model in TRAINING fits
        idx_params = [];
        for iFit = 1:numel(train_fits.modelfits.params)
            pp = train_fits.modelfits.params{iFit}; 
            if strcmp(pp.model_name,model_names_train{iModel})
                idx_params = iFit;
            end
        end

        if ~strcmp(filename,filename_train)    
            theta = train_fits.modelfits.params{idx_params}.theta;        
            params = params_new(model_names{iModel},data);
            params = setup_params(theta,params);
        else
            params = train_fits.modelfits.params{idx_params};
        end
    end
    
    % Plot
    subplot(plotrows,plotcols,iModel);
    plot_fit(data,params,params.model_desc);    
end

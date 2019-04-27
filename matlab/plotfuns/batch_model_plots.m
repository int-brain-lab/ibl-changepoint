function batch_model_plots(model_names,filename,model_names_train,filename_train)
%BATCH_MODEL_PLOTS Plot a batch of models results, possibly with transfer.

if nargin < 3 || isempty(model_names_train); model_names_train = model_names; end
if nargin < 4 || isempty(filename_train); filename_train = filename; end

close all;

% Models to be fitted
if isempty(model_names)
    model_names = get_model_list('default');
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

for iModel = 1:numel(model_names)
    if contains(model_names{iModel},'psychofun')
        % Look for psychometric model of the TEST data
        params = load_model_fit(filename,model_names{iModel});        
    else
        params_train = load_model_fit(filename_train,model_names_train{iModel});        
        if ~strcmp(filename,filename_train)
            % If training and test models are different, transfer parameters
            params = params_transfer(params_train,model_names{iModel},filename);            
        else
            params = params_train;
        end
    end
    % Plot
    subplot(plotrows,plotcols,iModel);
    plot_fit(data,params,params.model_desc);    
end

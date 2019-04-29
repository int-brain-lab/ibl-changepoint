function data = read_data_from_csv(fullname)
%READ_DATA_FROM_CSV Read from file and return DATA struct.

mypath = fileparts(mfilename('fullpath'));
addpath([mypath filesep '..' filesep 'data']);

% Add these modifiers to file name (preceded by underscore) to load a 
% specific subset of the data, e.g. 'CSHL_003_unbiased'
modifiers = {'unbiased','half1','half2'};

data_modifiers = {''};

% Look for data modifier tokens in the provided filename string
filename = fullname;
for iMod = 1:numel(modifiers)
    idx = strfind(filename,['_' modifiers{iMod}]);
    if ~isempty(idx)
        data_modifiers{end+1} = modifiers{iMod};
        filename = [filename(1:idx-1),filename(idx+numel(modifiers{iMod})+1:end)];
    end
end

filename_csv = [filename '.csv'];

% Read session data from CSV file (skip first row)
data_tab = csvread(filename_csv,1,0);

if any(cellfun(@(x)strcmp(x,'unbiased'),data_modifiers))
    % Remove biased trials
    fprintf('Removing biased trials...\n');
    data_tab(data_tab(:,3) ~= 0.5,:) = [];
end

if contains(fullname,'endtrain')
    fprintf('Removing biased trials...\n');
    data_tab(data_tab(:,3) ~= 0.5,:) = [];
end

% Keep only half among all sessions (to check for learning)
if any(cellfun(@(x)strcmp(x,'half1'),data_modifiers)) || ...
    any(cellfun(@(x)strcmp(x,'half2'),data_modifiers))    
    sessions = unique(data_tab(:,2))';
    first_half = sessions(1:floor(numel(sessions)/2));
    second_half = sessions(floor(numel(sessions)/2)+1:end);
    if any(cellfun(@(x)strcmp(x,'half1'),data_modifiers))
        fprintf('Keeping first half among all sessions...\n');
        idx = any(bsxfun(@eq,data_tab(:,2),first_half),2);
    else
        fprintf('Keeping second half among all sessions...\n');
        idx = any(bsxfun(@eq,data_tab(:,2),second_half),2);    
    end
    data_tab = data_tab(idx,:);
end    
    
data = format_data(data_tab,filename,fullname);

end
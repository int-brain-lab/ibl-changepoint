function data = read_data_from_csv(filename)
%READ_DATA_FROM_CSV Read from file and return DATA struct.

mypath = fileparts(mfilename('fullpath'));
addpath([mypath filesep 'data']);

% Check if loading only unbiased trials
idx = strfind(filename,'_unbiased');
if ~isempty(idx)
    unbiased_flag = true;
    filename_csv = [filename(1:idx-1) '.csv'];
else
    unbiased_flag = false;
    filename_csv = [filename '.csv'];
end

% Read session data from CSV file (skip first row)
data_tab = csvread(filename_csv,1,0);
if unbiased_flag    % Remove biased trials
    data_tab(data_tab(:,3) ~= 0.5,:) = [];
end

[~,name] = fileparts(filename);

data = format_data(data_tab,name);

end
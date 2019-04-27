function params = load_model_fit(data_filename,model_name)
%LOAD_MODEL_FIT Load fit for a given dataset and model.

mypath = fileparts(mfilename('fullpath'));
fits_path = [mypath filesep 'fits' filesep];
filename = [data_filename '_' model_name];

% Try loading file, if it does not exist or loading fails, return empty
try
    temp = load([fits_path filename '.mat'],'params');
    params = temp.params;
catch
    params = [];
end

end
function params = save_model_fit(data_filename,params)
%SAVE_MODEL_FIT Save fit for a given dataset and model.

mypath = fileparts(mfilename('fullpath'));
fits_path = [mypath filesep 'fits' filesep];

model_name = params.model_name;
filename = [data_filename '_' model_name];

root_name = get_mouse_name(data_filename);

% Create dataset folder if it does not exist
if ~(exist([fits_path root_name]) == 7)
    mkdir(fits_path,root_name);
end
save([fits_path root_name filesep filename '.mat'],'params');

end

%--------------------------------------------------------------------------
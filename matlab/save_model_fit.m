function params = save_model_fit(data_filename,params)
%SAVE_MODEL_FIT Save fit for a given dataset and model.

mypath = fileparts(mfilename('fullpath'));
fits_path = [mypath filesep 'fits' filesep];

model_name = params.model_name;
filename = [data_filename '_' model_name];

save([fits_path filename '.mat'],'params');

end

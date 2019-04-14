%IBL_CHANGE_ADD2PATH Temporarily repo full path to Matlab path.
%
%   Call this Matlab script at the beginning of a Matlab session, before 
%   running any file in the repository.
%
%   It is *not* recommended to permanently add the 'matlab' folder and sub-
%   folders of this repository to the Matlab path. Since Matlab does not 
%   have any way to deal with function duplicates, it might cause clashes 
%   with other projects.

% Luigi Acerbi 2019

repo_path = fileparts(mfilename('fullpath'));
matlab_path = [repo_path filesep 'matlab'];
addpath(matlab_path);

folders_list = {'data','figures','fits','plotfuns','scripts','utils'};

for iFolder = 1:numel(folders_list)
    folder_path = [matlab_path filesep folders_list{iFolder}];
    addpath(folder_path);
end

clear matlab_path folder_path iFolder folders_list matlab_path;
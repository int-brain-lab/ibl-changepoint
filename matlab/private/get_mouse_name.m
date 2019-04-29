function mouse_name = get_mouse_name(datafile)
%GET_MOUSE_NAME Get mouse name out of data file name.

mice_names = get_mice_list('all');

mouse_name = datafile;  % If match not found, return full name

for iMouse = 1:numel(mice_names)
    if strncmp(datafile,mice_names{iMouse},numel(mice_names{iMouse}))
        mouse_name = mice_names{iMouse};
    end
end



end
function mice_list = get_mice_list(id,suffix)
%GET_MICE_LIST Get cell array of standard mice.

if nargin < 1 || isempty(id); id = 'default'; end
if nargin < 2; suffix = []; end

all_mice = {'CSHL_003','CSHL_005','CSHL_007','CSHL_008','CSHL_010','IBL-T1','IBL-T4','ZM_1084','ZM_1085','ZM_1086','ZM_1091','ZM_1092','ZM_1093','ZM_1097','ZM_1098','ibl_witten_04','ibl_witten_05','ibl_witten_06'};

switch id
    case {0,'example'}    % Example mice        
        mice_list = {'CSHL_005','CSHL_007','IBL-T1','IBL-T4','ibl_witten_04','ibl_witten_05'};
        
    case {1,'default','all','apr2019'}        % All mice good for analysis
        mice_list = {'CSHL_003','CSHL_005','CSHL_007','CSHL_008','CSHL_010','IBL-T1','IBL-T4','ZM_1084','ZM_1085','ZM_1086','ZM_1091','ZM_1092','ZM_1093','ZM_1097','ZM_1098','ibl_witten_04','ibl_witten_05','ibl_witten_06'};

    otherwise
        if ischar(id) && any(strcmp(id,all_mice))            
            mice_list{1} = id;      % Single mouse
        else
            error('Unknown mice list identifier.');
        end    
end

if ~isempty(suffix)
   for ii = 1:numel(mice_list)
       mice_list{ii} = [mice_list{ii} '_' suffix];
   end
end
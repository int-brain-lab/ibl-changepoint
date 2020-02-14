function mice_list = get_mice_list(id,suffix)
%GET_MICE_LIST Get cell array of standard mice.

if nargin < 1 || isempty(id); id = 'default'; end
if nargin < 2; suffix = []; end

all_mice = {'CSHL_002','CSHL_003','CSHL_004','CSHL_005','CSHL_007','CSHL_008','CSHL_010','DY_001','IBL-T1','IBL-T2','IBL-T3','IBL-T4','NYU-02','ZM_1084','ZM_1085','ZM_1086','ZM_1091','ZM_1092','ZM_1093','ZM_1097','ZM_1098','ibl_witten_04','ibl_witten_05','ibl_witten_06','CSHL_007','CSHL_010','CSK-les-008','DY_006','IBL-T4','IBL_13','IBL_17','IBL_36','IBL_45','ibl_witten_06','ibl_witten_07','KS003','KS004','NYU-01','ZM_1089','ZM_1092','ZM_1093','ZM_1745','ZM_1746','CSHL_007','CSHL_010','CSK-les-008','DY_006','IBL-T4','IBL_13','IBL_17','IBL_36','IBL_45','ibl_witten_06','ibl_witten_07','KS003','KS004','NYU-01','ZM_1089','ZM_1092','ZM_1093','ZM_1745','ZM_1746'};

switch id
    case {0,'example'}    % Example mice        
        mice_list = {'CSHL_005','CSHL_007','IBL-T1','IBL-T4','ibl_witten_04','ibl_witten_05'};
        
    case {1,'default','all','may2019'}        % All mice good for analysis
        mice_list = all_mice;
        % {'CSHL_002','CSHL_003','CSHL_004','CSHL_005','CSHL_007','CSHL_008','CSHL_010','DY_001','IBL-T1','IBL-T2','IBL-T3','IBL-T4','NYU-02','ZM_1084','ZM_1085','ZM_1086','ZM_1091','ZM_1092','ZM_1093','ZM_1097','ZM_1098','ibl_witten_04','ibl_witten_05','ibl_witten_06'};

    case {2,'criteria_jul2019'}     % All mice that satisfy loose training criteria
        mice_list = {'CSHL_007', 'CSHL_010', 'CSK-les-008', 'DY_006', 'IBL-T4', 'IBL_13', 'IBL_17', 'IBL_36', 'IBL_45', 'ibl_witten_06', 'ibl_witten_07', 'KS003', 'KS004', 'NYU-01', 'ZM_1089', 'ZM_1092', 'ZM_1093', 'ZM_1745', 'ZM_1746', 'CSHL_007', 'CSHL_010', 'CSK-les-008', 'DY_006', 'IBL-T4', 'IBL_13', 'IBL_17', 'IBL_36', 'IBL_45', 'ibl_witten_06', 'ibl_witten_07', 'KS003', 'KS004', 'NYU-01', 'ZM_1089', 'ZM_1092', 'ZM_1093', 'ZM_1745', 'ZM_1746'};
        
    case {3,'strict_jul2019'}     % All mice that satisfy strict training criteria
        mice_list = {'CSHL_007', 'CSHL_010', 'CSK-les-008', 'DY_006', 'IBL-T4', 'IBL_13', 'IBL_17', 'IBL_36', 'IBL_45', 'ibl_witten_06', 'KS003', 'NYU-01', 'ZM_1092', 'ZM_1745', 'ZM_1746'};

    case {4,'strict_sep2019'}     % All mice that satisfy strict training criteria
        mice_list = {'CSHL_002', 'CSHL_003', 'CSHL_005', 'CSHL_008', 'CSHL_010', 'CSHL_014', 'CSHL_015', 'CSH_ZAD_001', 'CSH_ZAD_003', 'CSH_ZAD_004', 'CSH_ZAD_006', 'CSH_ZAD_007', 'CSH_ZAD_010', 'DY_001', 'DY_007', 'IBL-T1', 'IBL-T2', 'IBL-T4', 'IBL_001', 'IBL_002', 'KS002', 'KS003', 'KS004', 'KS005', 'KS014', 'KS015', 'KS016', 'KS017', 'NYU-01', 'NYU-02', 'NYU-06', 'SWC_013', 'ZM_1084', 'ZM_1085', 'ZM_1086', 'ZM_1087', 'ZM_1091', 'ZM_1092', 'ZM_1097', 'ZM_1098', 'ZM_1367', 'ZM_1371', 'ZM_1372', 'ZM_1743', 'ZM_1745', 'ZM_1746', 'ibl_witten_04', 'ibl_witten_05', 'ibl_witten_06', 'ibl_witten_12', 'ibl_witten_14', 'ibl_witten_15', 'ibl_witten_16'};
        
    case {'guido'}
        mice_list = {'ZM_1895_date20191127','ZM_1895_date20191128','ZM_1895_date20191129', ...
            'ZM_2102_date20191126','ZM_2102_date20191127','ZM_2102_date20191128', ...
            'ZM_2108_date20191204','ZM_2108_date20191205','ZM_2108_date20191206', ...
            'ZM_1895_date20191217','ZM_1895_date20191218','ZM_1895_date20191219', ...
            'ZM_2102_date20191218','ZM_2102_date20191219','ZM_2102_date20191220', ...
            'ZM_2108_date20191217','ZM_2108_date20191218','ZM_2108_date20191219' ...
            };
        
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

if nargout == 0
    for i = 1:numel(mice_list)
        fprintf('%s\n',mice_list{i});
    end
end
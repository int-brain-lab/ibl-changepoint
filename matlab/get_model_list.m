function model_names = get_model_list(id)
%GET_MODEL_LIST Get cell array of standard models.

switch id
    case {1,'omniscient','default'}
        model_names{1} = 'psychofun';
        model_names{2} = 'omniscient';
        model_names{3} = 'omniscient_lapse';
        model_names{4} = 'omniscient_biasedlapse';
        model_names{5} = 'omniscient_altnoise';
        
    case {2,'changepoint'}
        model_names{1} = 'psychofun';
        model_names{2} = 'changepoint';
        model_names{3} = 'changepoint_lapse';
        model_names{4} = 'changepoint_biasedlapse';
        model_names{5} = 'changepoint_altnoise';
end





%BATCH_PLOT_PSYCHOFUNS

mice_list = {'CSHL_010', 'CSK-les-008', 'DY_006', 'IBL-T4', 'IBL_13', 'IBL_17', 'IBL_36', 'IBL_45', 'ibl_witten_06', 'KS003', 'NYU-01', 'ZM_1092', 'ZM_1745', 'ZM_1746'};
    % get_mice_list('strict_jul2019');

plot_list = []; plot_sessions = []; sessinfo = []; plot_train_sessions = [];
theta = []; theta_train = [];
    
sessions_info = [[-1.,  0., 25., 31.]; ...
       [ 0., -1.,  1., -1.]; ...
       [-1.,  0.,  6.,  7.]; ...
       [-1., -1.,  0.,  3.]; ...
       [-1., -1.,  0.,  8.]; ...
       [-1., -1.,  0., -1.]; ...
       [-1., -1.,  0., -1.]; ...
       [-1., -1.,  0.,  4.]; ...
       [-1., -1.,  0.,  8.]; ...
       [-1., -1.,  0., 25.]; ...
       [ 0.,  6., 35., -1.]; ...
       [-1., -1.,  0.,  4.]; ...
       [-1., -1.,  0., -1.]; ...
       [-1., -1.,  0.,  3.]] + 1;
sessions_info(sessions_info == 0) = NaN;

training_stages = {'nan','1a','1b','eph'};

for iMouse = 1:numel(mice_list)
    mouse_name = mice_list{iMouse};
    data = read_data_from_csv(mouse_name);
    sessions = unique(data.tab(:,2));
    
    % Take only session that have all block types
    for iSession = 1:numel(sessions)
        idx = data.tab(:,2) == sessions(iSession);
        p = unique(data.tab(idx,3));
        if numel(p) ~= 3;  sessions(iSession) = NaN; end
    end    
    sessions = sessions(isfinite(sessions));   

    n = numel(sessions);
    if n < 3; continue; end
    
    data = read_data_from_csv([mouse_name '_endtrain']);
    train_sessions = unique(data.tab(:,2));    
    
    plot_list{end+1} = mouse_name;
    plot_sessions{end+1} = sessions;
    plot_train_sessions{end+1} = train_sessions;
    sessinfo{end+1} = sessions_info(iMouse,:);
end
        
for iMouse = 1:numel(plot_list)
    mouse_name = plot_list{iMouse};
    n = numel(plot_sessions{iMouse});

    for iSession = 1:n
        nSession = plot_sessions{iMouse}(iSession);
        % training_stage = find(nSession >= [sessinfo{iMouse},Inf],1,'last');
        data_name = [mouse_name '_sess' num2str(nSession)];
        data = read_data_from_csv(data_name);
        params = load_model_fit(data_name,'psychofun');        
        if isempty(params); continue; end
        for iParam = 1:4
            theta{iMouse,iParam}(iSession) = params.theta(4+iParam);
        end
    end

    for iSession = 1:3
        nSession = plot_train_sessions{iMouse}(iSession);
        data_name = [mouse_name '_endtrain_sess' num2str(nSession)];
        data = read_data_from_csv(data_name);
        params = load_model_fit(data_name,'psychofun');
        if isempty(params); continue; end
        for iParam = 1:4
            theta_train{iMouse,iParam}(iSession) = params.theta(iParam);
        end
    end    
end

for iParam = 1:4
    th{iParam} = []; thtrain{iParam} = [];
    for iMouse = 1:size(theta,1)
        th{iParam} = [th{iParam}; theta{iMouse,iParam}(:)];
        thtrain{iParam} = [thtrain{iParam}; theta_train{iMouse,iParam}(:)];
    end
    subplot(1,4,iParam);
    histogram(th{iParam}); hold on;
    histogram(thtrain{iParam}); hold on;
    %[~,d,xx] = kde(th{iParam}); plot(xx,d); hold on;
    %[~,d,xx] = kde(thtrain{iParam}); plot(xx,d);
end
    


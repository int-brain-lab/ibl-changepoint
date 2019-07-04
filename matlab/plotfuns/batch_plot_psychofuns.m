%BATCH_PLOT_PSYCHOFUNS

mice_list = {'CSHL_010', 'CSK-les-008', 'DY_006', 'IBL-T4', 'IBL_13', 'IBL_17', 'IBL_36', 'IBL_45', 'ibl_witten_06', 'KS003', 'NYU-01', 'ZM_1092', 'ZM_1745', 'ZM_1746'};
    % get_mice_list('strict_jul2019');

plot_list = []; plot_sessions = []; sessinfo = [];
    
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
    
    plot_list{end+1} = mouse_name;
    plot_sessions{end+1} = sessions;
    sessinfo{end+1} = sessions_info(iMouse,:);
end
    
Nmax = 6;
Nfigs = ceil(numel(plot_list)/Nmax);
    
for iMouse = 1:numel(plot_list)
    iFigure = floor((iMouse-1)/Nmax)+1;
    iRow = iMouse - (iFigure-1)*Nmax;    
    figure(iFigure);

    mouse_name = plot_list{iMouse}
    n = numel(plot_sessions{iMouse});
    idx = [round(n/3),round(n*2/3),n];

    nrows = Nmax;
    ncols = 6;

    sessinfo{iMouse}

    for iSession = 1:numel(idx)
        subplot(nrows,ncols,(iRow-1)*ncols + 3 + iSession);
        nSession = plot_sessions{iMouse}(idx(iSession))
        training_stage = find(nSession >= [sessinfo{iMouse},Inf],1,'last');
        data_name = [mouse_name '_sess' num2str(nSession)];
        data = read_data_from_csv(data_name);
        data.resp_obs(:) = NaN;
        params = load_model_fit(data_name,'psychofun');
        if isempty(params); continue; end
        plot_fit(data,params,[],false);
        legend off;
        if iRow < nrows || iSession > 1; xlabel(''); end
        if iSession > 1 || iRow < nrows; ylabel(''); end
        set(gca,'Ytick',[0 0.5 1]);
        ylim([0 1]);
        text(0.1,0.8,training_stages{training_stage},'Units','normalized','Fontsize',14);
    end

    for iSession = 1:3
        subplot(nrows,ncols,(iRow-1)*ncols + iSession);
        nSession = iSession;
        data_name = [mouse_name '_endtrain_sess' num2str(nSession)];
        data = read_data_from_csv(data_name);
        data.resp_obs(:) = NaN;
        params = load_model_fit(data_name,'psychofun');
        if iSession == 1
            text_name = mouse_name; text_name(text_name == '_') = '-';
            text(-1,0.5,text_name,'Units','normalized','Fontsize',14);
        end
        if isempty(params); continue; end
        plot_fit(data,params,[],false);
        legend off;
        if iRow < nrows || iSession > 1; xlabel(''); end
        if iSession > 1 || iRow < nrows; ylabel(''); end
        set(gca,'Ytick',[0 0.5 1]);
        ylim([0 1]);
    end    
end

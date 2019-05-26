function avgmu = blocklength_posterior()

mice = get_mice_list();
avgmu = zeros(1,numel(mice));

for iMouse = 1:numel(mice)
    data_name = mice{iMouse};
    fprintf('%s\n',data_name);
    data = read_data_from_csv(data_name);
    model_name = 'changepoint_nakarushton';
    params = load_model_fit(data_name,model_name);

    params.save_fullpost = true;    

    [~,output] = changepoint_bayesian_nll(params,data);

    post = sum(output.fullpost,3);  % Runlength posterior

    % Find change-points, remove changes of sessions
    changepoints = find(diff(data.p_true) ~= 0);
    changesessions = find(diff(data.tab(:,2)));
    changepoints = setdiff(changepoints,changesessions);

    t = 0:size(post,2)-1;

    for iChange = 1:numel(changepoints)
        postmu(iChange) = sum(post(changepoints(iChange),:).*t,2);    
    end
    avgmu(iMouse) = mean(postmu);
end


% Compute mean at average inferred tau from mice
mice_tau = 27;
min_rl = 20;
max_rl = 100;

idx = true(1,1e6);
a = zeros(size(idx));
while any(idx)
    a(idx) = exprnd(mice_tau,[1,sum(idx)])+min_rl;
    idx = a > max_rl;
end
mean(a)

end
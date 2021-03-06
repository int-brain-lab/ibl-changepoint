function data = format_data(data_tab,filename,fullname)
%FORMAT_DATA Create DATA struct from a data matrix and file name NAME.

data.tab = data_tab;
data.filename = filename;
data.fullname = fullname;

data.p_true = data_tab(:,3);
data.contrasts = data_tab(:,4);
data.S = data_tab(:,5);
data.resp_obs = data_tab(:,6);
data.resp_correct = data_tab(:,7);

data.mu = sort(unique(data.S(data.S ~= 0)));
data.C = (data.S > 0) + 1;    % Stimulus category (1 for L, 2 for R)

% Read contrast levels
data.contrasts_vec = sort(unique(data.contrasts));
data.contrasts_idx = data.contrasts;
for ii = 1:numel(data.contrasts_vec)
    idx = data.contrasts == data.contrasts_vec(ii);
    data.contrasts_idx(idx) = ii;
end
data.signed_contrasts = data.contrasts;
data.signed_contrasts(data.C == 1) = -data.signed_contrasts(data.C == 1);

signed_contrasts_vec = sort(unique(data.signed_contrasts));
for ii = 1:numel(signed_contrasts_vec)
    idx = data.signed_contrasts == signed_contrasts_vec(ii);
    data.signed_contrasts_idx(idx) = ii;
end

end

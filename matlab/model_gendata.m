function gendata = model_gendata(params,data,nreps)
%MODEL_GENDATA Generate simulated data from model fit.

if nargin < 3 || isempty(nreps); nreps = 1; end

[~,output] = nllfun([],params,data);

% Compute model average if passing multiple samples
Nparams = numel(output);
resp_model = zeros(size(output(1).resp_model,1),Nparams);
for iParam = 1:Nparams
    resp_model(:,iParam) = output(iParam).resp_model;
end
resp_model = mean(resp_model,2);

data_tab = repmat(data.tab,[nreps,1]);  % Extend dataset
% Correct session number
if nreps > 1
    data_tab(:,2) = cumsum(data_tab(:,1) == 1);
end
resp_model = repmat(resp_model,[nreps,1]);

% Predict Left and Right responses
resp_L = rand(size(resp_model,1),1) < resp_model;
data_tab(resp_L == 1,6) = -1; 
data_tab(resp_L == 0,6) = 1; 

% Correctness feedback
data_tab(:,7) = data_tab(:,6) == sign(data_tab(:,5));

% Reaction times are not modified

% Format simulated data
gendata = format_data(data_tab,data.name);

end
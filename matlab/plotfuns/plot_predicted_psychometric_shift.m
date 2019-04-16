function plot_predicted_psychometric_shift(example_mice)
%PLOT_PREDICTED_PSYCHOMETRIC_SHIFT Plot predicted psychometric functions.

if nargin < 1 || isempty(example_mice)
    example_mice = {'CSHL_005','CSHL_007','IBL-T1','IBL-T4','ibl_witten_04','ibl_witten_05'};
end
if ischar(example_mice); example_mice = {example_mice}; end

if numel(example_mice) > 1
   for iFig = 1:numel(example_mice)
       figure;
       plot_predicted_psychometric_shift(example_mice{iFig});
   end
   return;
end

plotrows = 2;
plotcols = 2;

mouse_name = example_mice{1};
mouse_string = mouse_name;
mouse_string(mouse_string == '_') = '-';

for iModel = 1:4
    subplot(plotrows,plotcols,iModel);
    if iModel == 1  % Real data
        temp = load([mouse_name '_fits.mat']);
        true_data = temp.modelfits.data;
        idx = find(cellfun(@(p) strcmp(p.model_name,'psychofun'),temp.modelfits.params),1);
        true_params = temp.modelfits.params{idx};
        title_string = ['data (' mouse_string ')'];
        plot_fit(true_data,true_params,title_string,0);
        true_data.resp_obs(:) = NaN;
    else            % Model-generated data
        temp = load([mouse_name '_bias_shift.mat']);
        data = temp.gendata_mle{iModel};
        params = temp.psy_model_mle{iModel};
        test_model = temp.test_models{iModel};
        test_model(test_model == '_') = '-';
        title_string = [test_model ' (simulated data + psychometric fit)'];
        % Plot data and psychometric curve
        plot_fit(true_data,true_params,title_string,0);
        plot_fit(data,params,title_string,0);
    end
    
end

filename = [mouse_name '_predictions'];
mypath = which('savefigure.m');
savefigure([fileparts(mypath) filesep() filename]);
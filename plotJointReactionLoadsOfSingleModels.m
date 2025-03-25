clear;
% dataPath = 'C:\Users\Willi\Documents\UniDataLokal\TorsionToolComparison\Simulationsresults';
dataPath = 'C:\Users\willi\ucloud\TreadmillTest\SimulationOutput';
load(fullfile(dataPath, 'dataStruct_ErrorScores_no_trials_removed.mat'));

%%

alpha = 0.1;

section = 'JRL';
fieldsToPlotLeft = {'hip_l_on_femur_l_in_femur_l_fx', 'hip_l_on_femur_l_in_femur_l_fy', 'hip_l_on_femur_l_in_femur_l_fz', 'patellofemoral_l_on_patella_l_in_patella_l_fx', 'patellofemoral_l_on_patella_l_in_patella_l_fy', 'patellofemoral_l_on_patella_l_in_patella_l_fz', 'walker_knee_l_on_tibia_l_in_tibia_l_fx', 'walker_knee_l_on_tibia_l_in_tibia_l_fy', 'walker_knee_l_on_tibia_l_in_tibia_l_fz', 'ankle_l_on_talus_l_in_talus_l_fx', 'ankle_l_on_talus_l_in_talus_l_fy', 'ankle_l_on_talus_l_in_talus_l_fz'};

fieldsToPlotRight = {'hip_r_on_femur_r_in_femur_r_fx', 'hip_r_on_femur_r_in_femur_r_fy', 'hip_r_on_femur_r_in_femur_r_fz', 'patellofemoral_r_on_patella_r_in_patella_r_fx', 'patellofemoral_r_on_patella_r_in_patella_r_fy', 'patellofemoral_r_on_patella_r_in_patella_r_fz', 'walker_knee_r_on_tibia_r_in_tibia_r_fx', 'walker_knee_r_on_tibia_r_in_tibia_r_fy', 'walker_knee_r_on_tibia_r_in_tibia_r_fz', 'ankle_r_on_talus_r_in_talus_r_fx', 'ankle_r_on_talus_r_in_talus_r_fy', 'ankle_r_on_talus_r_in_talus_r_fz'};

plotTitles = fieldsToPlotLeft;
for i = 1 : numel(fieldsToPlotRight)
    if strcmp(fieldsToPlotRight{i}(end-1 : end), '_l')
        fieldsToPlotRight{i} = [fieldsToPlotRight{i}(1 : end-2) '_r'];
        plotTitles{i} = plotTitles{i}(1 : end-2);
    end
end

factorsRight = ones(1, 12);
% specify which fields have to be inverted (e.g. pelvis rotation is
% inverted between left and right side

factorsRight([3 6 9 12]) = -1;

models = fieldnames(data.(section));

for i = 1 : numel(models)
    model = models{i};
    figure('Units', 'normalized', 'Position', [0.1 0.1 0.7 0.7]);
    tiledlayout(4, 3)
    sgtitle(model, 'Interpreter', 'none')

    for f = 1 : numel(fieldsToPlotLeft)
        nexttile(f); hold on;

        % f_getArrayForField returns the required data for the left and right steps
        [tmp_data, ~] = f_getArrayForField(data.(section).(model), fieldsToPlotLeft{f}, 0, [], 'walking');
        tmp_data = tmp_data / (data.(section).(model).model_mass * 9.81);
        stdshade(tmp_data, alpha, [1 0 0]);

        % f_getArrayForField returns the required data for the left and right steps
        [~, tmp_data] = f_getArrayForField(data.(section).(model), fieldsToPlotRight{f}, 0, [], 'walking');
        tmp_data = tmp_data / (data.(section).(model).model_mass * 9.81);
        stdshade(tmp_data * factorsRight(f), alpha, [0 0 1]);

        % f_getArrayForField returns the required data for the left and right steps
        [tmp_data, ~] = f_getArrayForField(data.(section).(model), fieldsToPlotLeft{f}, 0, [], 'running');
        tmp_data = tmp_data / (data.(section).(model).model_mass * 9.81);
        stdshade(tmp_data, alpha, [1 0.5 0]);

        % f_getArrayForField returns the required data for the left and right steps
        [~, tmp_data] = f_getArrayForField(data.(section).(model), fieldsToPlotRight{f}, 0, [], 'running');
        tmp_data = tmp_data / (data.(section).(model).model_mass * 9.81);
        stdshade(tmp_data * factorsRight(f), alpha, [0 0.5 1]);

        title(plotTitles{f}, 'Interpreter', 'none');
        ylabel('Force [BW]');
        xlabel('Gait Cycle [%]');
    end

    leg = legend({'', 'walking left', '', 'walking right', '', 'running left', '', 'running right'}, 'Orientation', 'Horizontal');
    leg.Layout.Tile = 'north';
end
clear;
dataPath = './ExampleOutput';
load(fullfile(dataPath, 'dataStruct_ErrorScores.mat'));
%%

alpha = 0.1;
section = 'IK';
fieldsToPlotLeft = {'pelvis_tilt', 'pelvis_list', 'pelvis_rotation', 'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l', 'knee_angle_l', 'ankle_angle_l', 'subtalar_angle_l'};
fieldsToPlotRight = fieldsToPlotLeft;
plotTitles = fieldsToPlotLeft;
for i = 1 : numel(fieldsToPlotRight)
    if strcmp(fieldsToPlotRight{i}(end-1 : end), '_l')
        fieldsToPlotRight{i} = [fieldsToPlotRight{i}(1 : end-2) '_r'];
        plotTitles{i} = plotTitles{i}(1 : end-2);
    end
end

factorsRight = ones(1, 9);
% specify which fields have to be inverted (e.g. pelvis rotation is
% inverted between left and right side

factorsRight([2 3]) = -1;

models = fieldnames(data.(section));
for i = 1 : numel(models)
    model = models{i};
    figure('Units', 'normalized', 'Position', [0.1 0.1 0.7 0.7]); 
    tiledlayout(3, 3)
    sgtitle(model, 'Interpreter', 'none')

    for f = 1 : numel(fieldsToPlotLeft)
        nexttile(f); hold on;
        % f_getArrayForField returns the required data for the left and right steps
        [tmp_data, ~] = f_getArrayForField(data.(section).(model), fieldsToPlotLeft{f});
        stdshade(tmp_data, alpha, [1 0 0]);
        title(plotTitles{f}, 'Interpreter', 'none');

        % f_getArrayForField returns the required data for the left and right steps
        [~, tmp_data] = f_getArrayForField(data.(section).(model), fieldsToPlotRight{f});
        stdshade(tmp_data * factorsRight(f), alpha, [0 0 1]);
    end

    leg = legend({'', 'LEFT', '', 'RIGHT'}, 'Orientation', 'Horizontal');
    leg.Layout.Tile = 'north';
end
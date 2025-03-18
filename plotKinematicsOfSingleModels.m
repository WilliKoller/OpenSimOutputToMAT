clear;
dataPath = 'C:\Users\willi\ucloud\TreadmillTest\SimulationOutput';
load(fullfile(dataPath, 'dataStruct_ErrorScores_no_trials_removed.mat'));
% load("participantData.mat");
%%

alpha = 0.1;
section = 'IK';
fieldsToPlotLeft = {'pelvis_tilt', 'pelvis_list', 'pelvis_rotation', 'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l', 'knee_angle_l', 'ankle_angle_l', 'calcn_l_Oy'};

fieldsToPlotRight = fieldsToPlotLeft;

plotTitles = fieldsToPlotLeft;
for i = 1 : numel(fieldsToPlotRight)
    if strcmp(fieldsToPlotRight{i}(end-1 : end), '_l')
        fieldsToPlotRight{i} = [fieldsToPlotRight{i}(1 : end-2) '_r'];
        plotTitles{i} = plotTitles{i}(1 : end-2);
    end
end
fieldsToPlotRight{9} = 'calcn_r_Oy';

factorsRight = ones(1, 9);
% specify which fields have to be inverted (e.g. pelvis rotation is
% inverted between left and right side

factorsRight([2 3 9]) = -1;

models = fieldnames(data.(section));
% models = models(contains(models, 'bdt'));
% models = models(contains(models, 'WT'));

for i = 1 : numel(models)
    model = models{i};
    figure('Units', 'normalized', 'Position', [0.1 0.1 0.7 0.7]);
    tiledlayout(3, 3)
    sgtitle(model, 'Interpreter', 'none')

    % participantID = str2double(model(3:4));
    % affectedLeg = particpantData(participantID, 2); % 1 == left; 2 == right

    % plot bdt model
    for f = 1 : numel(fieldsToPlotLeft)
        try
            nexttile(f); hold on;
            % f_getArrayForField returns the required data for the left and right steps
            % if affectedLeg == 1 % left
            [tmp_data, ~] = f_getArrayForField(data.(section).(model), fieldsToPlotLeft{f});
            stdshade(tmp_data, alpha, [1 0 0]);
            %     if f == 1
            %         disp([model ' affected leg is left; number of trials = ' num2str(size(tmp_data,1))]);
            %     end
            % else % right
            % f_getArrayForField returns the required data for the left and right steps
            [~, tmp_data] = f_getArrayForField(data.(section).(model), fieldsToPlotRight{f});
            stdshade(tmp_data * factorsRight(f), alpha, [0 0 1]);
            % if f == 1
            %     disp([model ' affected leg is right; number of trials = ' num2str(size(tmp_data,1))]);
            % end
            % end
            title(plotTitles{f}, 'Interpreter', 'none');
            ylabel('Angle[°]');
            xlabel('Gait cycle[%]');
        end
    end

    % plot torsion tool model
    model = strrep(models{i}, 'bdt_', '');
    for f = 1 : numel(fieldsToPlotLeft)
        nexttile(f); hold on;
        % f_getArrayForField returns the required data for the left and right steps
        if affectedLeg == 1 % left
            [tmp_data, ~] = f_getArrayForField(data.(section).(model), fieldsToPlotLeft{f});
            stdshade(tmp_data, alpha, [1 1 0]);
        else % right
            % f_getArrayForField returns the required data for the left and right steps
            [~, tmp_data] = f_getArrayForField(data.(section).(model), fieldsToPlotRight{f});
            stdshade(tmp_data * factorsRight(f), alpha, [0 1 1]);
        end
    end

    % plot generic model
    model = strrep(models{i}, 'bdt_scaled', 'WT');
    for f = 1 : numel(fieldsToPlotLeft)
        nexttile(f); hold on;
        % f_getArrayForField returns the required data for the left and right steps
        if affectedLeg == 1 % left
            [tmp_data, ~] = f_getArrayForField(data.(section).(model), fieldsToPlotLeft{f});
            stdshade(tmp_data, alpha, [0.4940 0.1840 0.5560]);
        else % right
            % f_getArrayForField returns the required data for the left and right steps
            [~, tmp_data] = f_getArrayForField(data.(section).(model), fieldsToPlotRight{f});
            stdshade(tmp_data * factorsRight(f), alpha, [0.6350 0.0780 0.1840]);
        end
    end

    if affectedLeg == 1
        leg = legend({'', 'LEFT bdt', '', 'left TT', '', 'left Generic'}, 'Orientation', 'Horizontal');
    else
        leg = legend({'', 'RIGHT bdt','', 'right TT','', 'right Generic'}, 'Orientation', 'Horizontal');
    end
    leg.Layout.Tile = 'north';
end
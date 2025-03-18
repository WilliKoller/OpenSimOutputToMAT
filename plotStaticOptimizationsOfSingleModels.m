clear;
dataPath = 'C:\Users\Willi\Documents\UniDataLokal\TorsionToolComparison\Simulationsresults';
load(fullfile(dataPath, 'dataStruct_ErrorScores.mat'));
load('participantData.mat');
%%

alpha = 0.1;

section = 'SO_Activation';
fields = fieldnames(data.SO_Activation.ID01_bdt_scaled_final.T_ID01g1_4_1_right);
fields = fields(1 : 81);
fieldsToPlotLeft = {};
for f = 1 : numel(fields)
    if strcmp(fields{f}(end-1:end), '_l')        %wenn die letzte und vorletzte stelle mit _ und l enden dann macht er es ins array
        fieldsToPlotLeft{end+1} = fields{f};     %end+1 macht das array von 0 auf 1 und füllt es dann mit den muskeln mit der Bedingiung
    end
end
% fieldsToPlotLeft = fields(fields, '_l'))
% fieldsToPlotLeft = {'glmax1_l', 'glmax2_l', 'glmax3_l', 'recfem_l', 'vaslat_l', 'vasmed_l', 'semiten_l', 'semimem_l', 'bflh_l'};

fieldsToPlotRight = strrep(fieldsToPlotLeft, '_l', '_r');

plotTitles = fieldsToPlotLeft;
for i = 1 : numel(fieldsToPlotRight)
    if strcmp(fieldsToPlotRight{i}(end-1 : end), '_l')
        fieldsToPlotRight{i} = [fieldsToPlotRight{i}(1 : end-2) '_r'];
        plotTitles{i} = plotTitles{i}(1 : end-2);
    end
end

factorsRight = ones(1, 40);
% specify which fields have to be inverted (e.g. pelvis rotation is
% inverted between left and right side

%factorsRight([2 3 9]) = -1;

models = fieldnames(data.(section));
models = models(contains(models, 'bdt'));
% models = models(contains(models, 'WT'));

for i = 1 : 2%numel(models)
    model = models{i}; 
    figure('Units', 'normalized', 'Position', [0.1 0.1 0.7 0.7]); 
    tiledlayout(6, 7)
    sgtitle(model, 'Interpreter', 'none')

    participantID = str2double(model(3:4));
    affectedLeg = particpantData(participantID, 2); % 1 == left; 2 == right

    % plot bdt model
    for f = 1 : numel(fieldsToPlotLeft)
        nexttile(f); hold on;
        % f_getArrayForField returns the required data for the left and right steps
        if affectedLeg == 1 % left
            [tmp_data, ~] = f_getArrayForField(data.(section).(model), fieldsToPlotLeft{f});
            stdshade(tmp_data, alpha, [1 0 0]);
        else % right
            % f_getArrayForField returns the required data for the left and right steps
            [~, tmp_data] = f_getArrayForField(data.(section).(model), fieldsToPlotRight{f});
            stdshade(tmp_data * factorsRight(f), alpha, [0 0 1]);
        end
        title(plotTitles{f}, 'Interpreter', 'none');
        ylabel('Activation [0-1]');
        xlabel('Gait cycle [%]');
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


    for f = 1 : numel(fieldsToPlotLeft)
        model_bdt = models{i}; 
        modelGeneric = strrep(models{i}, 'bdt_scaled', 'WT');
        modelTT = strrep(models{i}, '_bdt', '');
        if affectedLeg == 1 % left
            [tmp_data_bdt, ~] = f_getArrayForField(data.(section).(model_bdt), fieldsToPlotLeft{f});
            [tmp_data_tt, ~] = f_getArrayForField(data.(section).(modelTT), fieldsToPlotLeft{f});
            [tmp_data_generic, ~] = f_getArrayForField(data.(section).(modelGeneric), fieldsToPlotLeft{f});
        else % right
            % f_getArrayForField returns the required data for the left and right steps
            [~, tmp_data_bdt] = f_getArrayForField(data.(section).(model_bdt), fieldsToPlotRight{f});
            [~, tmp_data_tt] = f_getArrayForField(data.(section).(modelTT), fieldsToPlotRight{f});
            [~, tmp_data_generic] = f_getArrayForField(data.(section).(modelGeneric), fieldsToPlotRight{f});
        end

        % calculate RMS to generic
        mean_bdt = mean(tmp_data_bdt);
        mean_tt = mean(tmp_data_tt);
        mean_generic = mean(tmp_data_generic);

        rms_bdt(i, f) = sqrt(mean((mean_bdt - mean_generic).^2));
        rms_tt(i, f) = sqrt(mean((mean_tt - mean_generic).^2));
        rms_bdt_tt(i, f) = sqrt(mean((mean_tt - mean_bdt).^2));
    end
end

meanRMS_bdt = mean(rms_bdt, 1);
meanRMS_tt = mean(rms_tt, 1);
meanRMS_bdt_tt = mean(rms_bdt_tt, 1);

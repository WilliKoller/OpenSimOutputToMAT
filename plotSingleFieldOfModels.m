clear;
dataPath = '/Users/niklasnagel/Documents/Simulationsresults';
load(fullfile(dataPath, 'dataStruct_ErrorScores.mat'));
%%
% fieldToCompareLeft = 'knee_angle_l';
% fieldToCompareRight = 'knee_angle_r';
% section = 'IK';

fieldToCompareLeft = 'hip_flexion_l_moment';
fieldToCompareRight = 'hip_flexion_r_moment';
section = 'ID';

% fieldToCompareLeft = 'recfem_r';
% fieldToCompareRight = 'recfem_l';
% section = 'SO';

% fieldToCompareLeft = 'hip_r_on_femur_r_in_femur_r_fx';
% fieldToCompareRight = 'hip_l_on_femur_l_in_femur_l_fx';
% section = 'JRL';

alpha = 0.1;

models = fieldnames(data.(section));
for i = 1 : numel(models)
    model = models{i};
    figure('Position', [488,342,1038,420]); 
    tiledlayout(1, 2)
    sgtitle(model, 'Interpreter', 'none')
    
    % f_getArrayForField returns the required data for the left and right steps
    [tmp_data, ~] = f_getArrayForField(data.(section).(model), fieldToCompareLeft);
    nexttile; hold on;
    stdshade(tmp_data, alpha, [0 0 0]);
    title('left');
    
    % f_getArrayForField returns the required data for the left and right steps
    [~, tmp_data] = f_getArrayForField(data.(section).(model), fieldToCompareRight);
    nexttile; hold on;
    stdshade(tmp_data, alpha, [0 0 0]);
    title('right');

    leg = legend({'std', 'mean'}, 'Orientation', 'Horizontal');
    leg.Layout.Tile = 'north';
end
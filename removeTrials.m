%% load data
clear;
dataPath = './ExampleOutput';
load(fullfile(dataPath, 'dataStruct_ErrorScores_no_trials_removed.mat'));

%% remove trials with ErrorScore > 2 --> IK error or Simulation Error
disp('Remove trials with IK error');
if isfield(data.ErrorScoreList, 'ErrorScore3')
    for i = 1 : size(data.ErrorScoreList.ErrorScore3, 2)
        trial2remove = data.ErrorScoreList.ErrorScore3{i};
        data = removeFromStruct(data, trial2remove);
        disp(['removed ', trial2remove]);
    end
end
%%
disp('Remove trials with simulation error');
if isfield(data.ErrorScoreList, 'ErrorScore4')
    for i = 1 : size(data.ErrorScoreList.ErrorScore4, 2)
        trial2remove = data.ErrorScoreList.ErrorScore4{i};
        data = removeFromStruct(data, trial2remove);
        disp(['removed ', trial2remove]);
    end
end
% %% code to exclude one model (only 2 steps each side)
% disp('exclude TD06 (only 2 steps each side)');
% modelToExclude = 'TD06_generic_final';
% data.IK = rmfield(data.IK, modelToExclude);
% data.ID = rmfield(data.ID, modelToExclude);
% data.SO = rmfield(data.SO, modelToExclude);
% data.SO_Activation = rmfield(data.SO_Activation, modelToExclude);
% data.JRL = rmfield(data.JRL, modelToExclude);
% disp(['removed ' modelToExclude]);

%% overwrite .mat

delete([dataPath '\dataStruct_ErrorScores.mat']);
dataFile = [dataPath '\dataStruct_ErrorScores.mat'];
save(dataFile, 'data', 'Scores_table')
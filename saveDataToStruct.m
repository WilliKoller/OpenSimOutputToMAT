%%
clear;
% outputPath = 'C:\Users\Willi\ucloud\PhD\Study_LongitudinalMSK\SimulationOutput_Sanguex_GRF20Hz';
outputPath = 'C:\Users\Biomechanik\SynologyDrive\AlexP\SimulationOutputCatelli';

modelList = GetSubDirsFirstLevelOnly(outputPath);

mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);

ErrorScores = zeros(4,1);

electromechanicalDelay = 0.1; % in seconds

% markerNamesToIgnoreWhenCheckingForErrors = {'C7'}; % example of markers
% that you want to ignore
markerNamesToIgnoreWhenCheckingForErrors = {};
maxMarkerErrorThreshold = 0.05;
maxMarkerRMSErrorThreshold = 0.03;
ignoreContralateralSide = 0;
markersToCheckBothSides = {'RASI', 'LASI', 'LPSI', 'RPSI'};

maxReserveActivation = 1; % set to -1 to ignore
ignoreContrallateralReserves = 0;
reservesToIgnoreBothSides = {' '}; % do not set to '';
importTrialEvenIfReservesAreHigh = 1; % 0 skips also the import of IK and ID results; 1 reads all data and adds data.SO_Activation.(trial).reserversAreBelowThreshold = 0 flag

%%
for p = 1 : numel(modelList)
    disp(' ');
    probandFolder = [outputPath filesep modelList{p}];
    modelList{p} = strrep(modelList{p}, '-', '_');
    trialList = GetSubDirsFirstLevelOnly(probandFolder);
    for trialNr = 1 : numel(trialList)
        disp(['currently processing:' ' ' modelList{p} ' ' '|' ' ' trialList{trialNr}]);

        currentFolder = [probandFolder filesep trialList{trialNr}];

        clear preframes EMG EMGChannels ratioEmgToFrames emgFrequency;
        load(fullfile(currentFolder, 'settings.mat'));

        cycleIsSorted = 0;
        if isfield(cycle, 'left') && isfield(cycle, 'right')
            if issorted(cycle.left.start) && issorted(cycle.right.start)
                cycleIsSorted = 1;
            end
        elseif isfield(cycle, 'left')
            if issorted(cycle.left.start)
                cycleIsSorted = 1;
            end
        else
            if issorted(cycle.right.start)
                cycleIsSorted = 1;
            end
        end
        if cycleIsSorted

            ikFile = fullfile(currentFolder, 'Output', 'IK', 'IK.mot');
            if isfile(ikFile)


                fid = fopen(fullfile(currentFolder, 'Output', 'IK', 'ikSettings.xml'), 'r');
                f=fread(fid,'*char')';
                fclose(fid);
                i0 = strfind(f, '<marker_file>');
                i1 = strfind(f, '</marker_file>');

                markerFileName = f(i0 + 13 : i1 - 1);

                % markerFileName = strrep(markerFileName, 'C:\Users\Biomechanik\Desktop\Willi\Study_LongitudinalMSK', 'C:\Users\Willi\ucloud\PhD\Study_LongitudinalMSK');

                markerData = load_marker_trc(markerFileName);
                markerNames = fieldnames(markerData);
                heelMarker = markerNames(contains(markerNames, 'HEE'));
                heelMarker = heelMarker{1}(1 : end-2);
                heelMarker = strrep(heelMarker, 'LH', 'H');
                heelMarker = strrep(heelMarker, 'RH', 'H');

                leftHeel = [cell2mat(markerData.(['L' heelMarker '_X'])), cell2mat(markerData.(['L' heelMarker '_Y'])), cell2mat(markerData.(['L' heelMarker '_Z']))];
                rightHeel = [cell2mat(markerData.(['R' heelMarker '_X'])), cell2mat(markerData.(['R' heelMarker '_Y'])), cell2mat(markerData.(['R' heelMarker '_Z']))];

                if isfield(cycle, 'left') && isfield(cycle, 'right')
                    frameZero = min(min(cycle.left.start), min(cycle.right.start)) - 1;
                elseif isfield(cycle, 'left')
                    frameZero = min(cycle.left.start) - 1;
                else
                    frameZero = min(cycle.right.start) - 1;
                end

                if exist('preframes', 'var')
                    preframes = preframes + 1;
                    frameZero = frameZero - floor(preframes); % this is a fix for the SO errors --> simulation started a few frames earlier to avoid activation limit
                end

                if isfile(fullfile(currentFolder, 'Output', 'IK', 'IK_model_marker_locations.sto'))

                    inputMarkerTime = cell2mat(markerData.Time);
                    markerLocationsIK = load_sto_file(fullfile(currentFolder, 'Output', 'IK', 'IK_model_marker_locations.sto'));
                    markerLocationsIKTime = markerLocationsIK.time;
                    offsetFrames = find(inputMarkerTime >= markerLocationsIKTime(1), 1);

                    markerNamesIK = fieldnames(markerLocationsIK);
                    markerNamesIK = markerNamesIK(~contains(markerNamesIK, 'time'));
                    for m = 1 : numel(markerNamesIK)
                        markerNamesIK{m} = markerNamesIK{m}(1 : end-3);
                    end
                    markerNamesIK = unique(markerNamesIK);

                    markerNamesIK = markerNamesIK(~ismember(markerNamesIK, markerNamesToIgnoreWhenCheckingForErrors));

                    markersIdxLeftToIgnore = []; markersIdxRightToIgnore = [];
                    for m = 1 : numel(markerNamesIK)
                        markerName = markerNamesIK{m};
                        if ~ismember(markerName, markersToCheckBothSides)
                            if markerName(1) == 'L'
                                markersIdxLeftToIgnore(end+1) = m;
                            elseif markerName(1) == 'R'
                                markersIdxRightToIgnore(end+1) = m;
                            end
                        end
                        ik_loc = [markerLocationsIK.([markerName '_tx']) markerLocationsIK.([markerName '_ty']) markerLocationsIK.([markerName '_tz'])];

                        % inputMarkerName = [markerNameIK(1:end-2) upper(markerNameIK(end))];
                        input_loc = [cell2mat(markerData.([markerName '_X'])) cell2mat(markerData.([markerName '_Y'])) cell2mat(markerData.([markerName '_Z']))];
                        if range(ik_loc(:, 1) > range(input_loc(:, 1)) * 50 )
                            ik_loc = ik_loc ./ 1000;
                        elseif range(input_loc(:, 1)) > range(ik_loc(:, 1)) * 50
                            input_loc = input_loc ./ 1000;
                        end

                        input_loc = input_loc(offsetFrames : end, :);

                        for t = 1 : size(ik_loc, 1)
                            markerErrors(m, t) = norm(ik_loc(t, :) - input_loc(t, :));
                        end
                    end

                    for t = 1 : size(ik_loc, 1)
                        markerRMSErrors(t) = rms(markerErrors(:, t));
                    end

                    if isfield(cycle, 'left')
                        for j = 1 : size(cycle.left.start, 2)
                            markerErrorsInThisCycle = markerErrors(:, cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                            markerRMSErrorsInThisCycle = markerRMSErrors(:, cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                            timeOfThisCycle = markerLocationsIKTime(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                            tmpMarkerNames = markerNamesIK;
                            if ignoreContralateralSide
                                markerErrorsInThisCycle(markersIdxRightToIgnore, :) = [];
                                tmpMarkerNames(markersIdxRightToIgnore) = [];
                            end
                            [maxValue, ind] = max(markerErrorsInThisCycle, [], 'all');
                            [row, col] = ind2sub(size(markerErrorsInThisCycle), ind);
                            if maxValue > maxMarkerErrorThreshold
                                disp(['Left cycle ' num2str(j) ' not valid, Marker ' tmpMarkerNames{row} ' has error of ' num2str(maxValue) ' at time ' num2str(timeOfThisCycle(col)) ' s']);
                                cycle.left.valid(j) = 0;
                            else
                                [maxValue, ind] = max(markerRMSErrorsInThisCycle, [], 'all');
                                [row, col] = ind2sub(size(markerRMSErrorsInThisCycle), ind);
                                if maxValue > maxMarkerRMSErrorThreshold
                                    disp(['Left cycle ' num2str(j) ' not valid, RMS Marker error is ' num2str(max(markerRMSErrorsInThisCycle)) ' at time ' num2str(timeOfThisCycle(col)) ' s']);
                                    cycle.left.valid(j) = 0;
                                else
                                    cycle.left.valid(j) = 1;
                                end
                            end
                        end
                    end

                    if isfield(cycle, 'right')
                        for j = 1 : size(cycle.right.start, 2)
                            markerErrorsInThisCycle = markerErrors(:, cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                            markerRMSErrorsInThisCycle = markerRMSErrors(:, cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                            timeOfThisCycle = markerLocationsIKTime(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                            tmpMarkerNames = markerNamesIK;
                            if ignoreContralateralSide
                                markerErrorsInThisCycle(markersIdxLeftToIgnore, :) = [];
                                tmpMarkerNames(markersIdxLeftToIgnore) = [];
                            end
                            [maxValue, ind] = max(markerErrorsInThisCycle, [], 'all');
                            [row, col] = ind2sub(size(markerErrorsInThisCycle), ind);
                            if maxValue > maxMarkerErrorThreshold
                                disp(['Right cycle ' num2str(j) ' not valid, Marker ' tmpMarkerNames{row} ' has error of ' num2str(maxValue) ' at time ' num2str(timeOfThisCycle(col)) ' s']);
                                cycle.right.valid(j) = 0;
                            else
                                [maxValue, ind] = max(markerRMSErrorsInThisCycle, [], 'all');
                                [row, col] = ind2sub(size(markerRMSErrorsInThisCycle), ind);
                                if maxValue > maxMarkerRMSErrorThreshold
                                    disp(['Right cycle ' num2str(j) ' not valid, RMS Marker error is ' num2str(max(markerRMSErrorsInThisCycle)) ' at time ' num2str(timeOfThisCycle(col)) ' s']);
                                    cycle.right.valid(j) = 0;
                                else
                                    cycle.right.valid(j) = 1;
                                end
                            end
                        end
                    end
                else
                    % this is just if the model's marker locations have not
                    % been reported in the IK step. Then the whole trial
                    % needs to be treated depending on "ErrorScore"

                    % check for IK_fehler.txt bzw IK_withTorsoOutOfRecommendation.txt
                    % ErrorScore 1: Alles cool
                    % ErrorScore 2: Nur Marker errors am Oberkörper zu hoch
                    % ErrorScore 3: Marker error bei Inverse kinematic zu hoch
                    % ErrorScore 4: Fehler bei Simulation

                    if isfile([currentFolder '\IK_fehler.txt']) && isfile([currentFolder '\IK_withTorsoOutOfRecommendation.txt']) && ~isfile([currentFolder '\IK_noProb.txt'])
                        % beide files vorhanden --> ErrorScore 3 --> bad
                        ErrorScores(3,1) = ErrorScores(3,1)+1;
                        data.ErrorScoreList.ErrorScore3{ErrorScores(3,1)} = currentFolder;
                    elseif  isfile([currentFolder '\IK_withTorsoOutOfRecommendation.txt']) && ~isfile([currentFolder '\IK_fehler.txt'])
                        % nur IK_withTorso... vorhanden --> ErrorScore 2 --> ok
                        ErrorScores(2,1) = ErrorScores(2,1)+1;
                        data.ErrorScoreList.ErrorScore2{ErrorScores(2,1)} = currentFolder;
                    else
                        % nichts vorhanden --> ErrorScore 1 --> perfect
                        ErrorScores(1,1) = ErrorScores(1,1)+1;
                        data.ErrorScoreList.ErrorScore1{ErrorScores(1,1)} = currentFolder;
                    end

                    cycle.right.valid(1 : size(cycle.right.start, 2)) = 1;
                    cycle.left.valid(1 : size(cycle.left.start, 2)) = 1;
                end

                if maxReserveActivation > 0
                    soFile = fullfile(currentFolder, 'Output', 'SO', '_StaticOptimization_activation.sto');
                    if isfile(soFile)
                        so_tempData = load_sto_file(soFile);
                        tempFieldNames = fieldnames(so_tempData);
                        reserves = tempFieldNames(contains(tempFieldNames, '_reserve'));
                        reserves = reserves(~contains(reserves, 'lumbar'));
                        if isfield(cycle, 'left')
                            for j = 1 : size(cycle.left.start, 2)
                                if cycle.left.valid(j)
                                    cycle.left.validMuscleActivation(j) = 1;
                                    for r = 1 : numel(reserves)
                                        if ~contains(reserves{r}, reservesToIgnoreBothSides)
                                            if ~(ignoreContrallateralReserves && contains(reserves{r}, '_r_'))
                                                reserveData = so_tempData.(reserves{r})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                                                % disp(reserves{r});
                                                [maxValue, ind] = max(reserveData);
                                                if maxValue > maxReserveActivation
                                                    if importTrialEvenIfReservesAreHigh == 0
                                                        cycle.left.valid(j) = 0;
                                                    end
                                                    cycle.left.validMuscleActivation(j) = 0;
                                                    disp(['Left cycle ' num2str(j) ' muscle activations not valid, ' reserves{r} ' is higher ( ' num2str(maxValue, 2) ' ) than max threshold ( ' num2str(maxReserveActivation) ' )']);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        if isfield(cycle, 'right')
                            for j = 1 : size(cycle.right.start, 2)
                                if cycle.right.valid(j)
                                    cycle.right.validMuscleActivation(j) = 1;
                                    for r = 1 : numel(reserves)
                                        if ~contains(reserves{r}, reservesToIgnoreBothSides)
                                            if ~(ignoreContrallateralReserves && contains(reserves{r}, '_l_'))
                                                reserveData = so_tempData.(reserves{r})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                                                % disp(reserves{r});
                                                [maxValue, ind] = max(reserveData);
                                                if maxValue > maxReserveActivation
                                                    if importTrialEvenIfReservesAreHigh == 0
                                                        cycle.right.valid(j) = 0;
                                                    end
                                                    cycle.right.validMuscleActivation(j) = 0;
                                                    disp(['Right cycle ' num2str(j) ' muscle activations not valid, ' reserves{r} ' is higher (' num2str(maxValue, 2) ') than max threshold ( = ' num2str(maxReserveActivation) ' )']);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end


                %crop data from left and right leg to stance phase
                tempData = load_sto_file(ikFile);
                tempFieldNames = fieldnames(tempData);

                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        if cycle.left.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructLeft;

                            heelPositionStrike1 = leftHeel(cycle.left.start(j) - frameZero, :);
                            heelPositionStrike2 = leftHeel(cycle.left.end(j) - frameZero - 1, :);

                            % rotate pelvis rotation into walking direction
                            walkingDir = heelPositionStrike2 - heelPositionStrike1;
                            walkingDirAngle = atan2d(walkingDir(3), walkingDir(1));
                            for i = 1 : numel(tempFieldNames)
                                if strcmp(tempFieldNames{i}, 'pelvis_rotation')
                                    tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                                    tempStructLeft.([tempFieldNames{i} 'adjWalkingDir']) = tempStructLeft.(tempFieldNames{i}) + walkingDirAngle;
                                else
                                    tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                                end
                            end

                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = tempStructLeft;
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).walkingDirection = walkingDir;
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).strideLength = norm(walkingDir);
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds = double(cycle.left.end(j) - cycle.left.start(j)) / frequency;
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).walkingSpeed = norm(walkingDir) / data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds;
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).footOffFrame = double(cycle.left.footOff(j) - cycle.left.start(j));
                        end
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        if cycle.right.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructRight;

                            heelPositionStrike1 = rightHeel(cycle.right.start(j) - frameZero, :);
                            heelPositionStrike2 = rightHeel(cycle.right.end(j) - frameZero - 1, :);

                            % rotate pelvis rotation into walking direction
                            walkingDir = heelPositionStrike2 - heelPositionStrike1;
                            walkingDirAngle = atan2d(walkingDir(3), walkingDir(1));

                            for i = 1 : numel(tempFieldNames)
                                if strcmp(tempFieldNames{i}, 'pelvis_rotation')
                                    tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                                    tempStructRight.([tempFieldNames{i} 'adjWalkingDir']) = tempStructRight.(tempFieldNames{i}) - walkingDirAngle;
                                    if mean(tempStructRight.(tempFieldNames{i}) - walkingDirAngle) > 150
                                        tempStructRight.([tempFieldNames{i} 'adjWalkingDir']) = tempStructRight.(tempFieldNames{i}) - walkingDirAngle - 180;
                                    end
                                else
                                    tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                                end
                            end

                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = tempStructRight;
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).walkingDirection = walkingDir;
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).strideLength = norm(walkingDir);
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).durationInSeconds = double(cycle.right.end(j) - cycle.right.start(j)) / frequency;
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).walkingSpeed = norm(walkingDir) / data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).durationInSeconds;
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).footOffFrame = double(cycle.right.footOff(j) - cycle.right.start(j));
                        end
                    end
                end
            end


            idFile = fullfile(currentFolder, 'Output', 'ID', 'inverse_dynamics.sto');
            if isfile(idFile)
                data.ID.(modelList{p}).model_mass = model_mass;
                %crop data from left and right leg to stance phase
                tempData = load_sto_file(idFile);
                tempFieldNames = fieldnames(tempData);

                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        if cycle.left.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructLeft;
                            for i = 1 : numel(tempFieldNames)
                                tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                            end
                            data.ID.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = tempStructLeft;
                            data.ID.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds = double(cycle.left.end(j) - cycle.left.start(j)) / frequency;
                        end
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        if cycle.right.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructRight;
                            for i = 1 : numel(tempFieldNames)
                                tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                            end
                            data.ID.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = tempStructRight;
                            data.ID.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).durationInSeconds = double(cycle.right.end(j) - cycle.right.start(j)) / frequency;
                        end
                    end
                end
            end


            soFile = fullfile(currentFolder, 'Output', 'SO', '_StaticOptimization_force.sto');
            if isfile(soFile)
                data.SO.(modelList{p}).model_mass = model_mass;
                %crop data from left and right leg to stance phase
                tempData = load_sto_file(soFile);
                tempFieldNames = fieldnames(tempData);
                
                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        if cycle.left.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructLeft;
                            for i = 1 : numel(tempFieldNames)
                                tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                            end
                            data.SO.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = tempStructLeft;
                            data.SO.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds = double(cycle.left.end(j) - cycle.left.start(j)) / frequency;
                        end
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        if cycle.right.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructRight;
                            for i = 1 : numel(tempFieldNames)
                                tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                            end
                            data.SO.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = tempStructRight;
                            data.SO.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).durationInSeconds = double(cycle.right.end(j) - cycle.right.start(j)) / frequency;
                        end
                    end
                end
            end

            soFile = fullfile(currentFolder, 'Output', 'SO', '_StaticOptimization_activation.sto');
            if isfile(soFile)
                data.SO_Activation.(modelList{p}).model_mass = model_mass;
                %crop data from left and right leg to stance phase
                if maxReserveActivation > 0 % SO file has already been loaded
                    tempData = so_tempData;
                else % read SO data from file
                    tempData = load_sto_file(soFile);
                end
                tempFieldNames = fieldnames(tempData);
                
                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        if cycle.left.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructLeft;
                            for i = 1 : numel(tempFieldNames)
                                tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                            end
                            data.SO_Activation.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = tempStructLeft;
                            data.SO_Activation.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).reserversAreBelowThreshold = cycle.left.validMuscleActivation(j);
                            data.SO_Activation.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds = double(cycle.left.end(j) - cycle.left.start(j)) / frequency;                            
                        end
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        if cycle.right.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructRight;
                            for i = 1 : numel(tempFieldNames)
                                tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                            end
                            data.SO_Activation.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = tempStructRight;
                            data.SO_Activation.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).reserversAreBelowThreshold = cycle.right.validMuscleActivation(j);
                            data.SO_Activation.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).durationInSeconds = double(cycle.right.end(j) - cycle.right.start(j)) / frequency;
                        end
                    end
                end
            end

            % EMG stuff
            fid = fopen(fullfile(currentFolder, 'Output', 'IK', 'ikSettings.xml'), 'r');
            f=fread(fid,'*char')';
            fclose(fid);
            i0 = strfind(f, '<marker_file>');
            i1 = strfind(f, '</marker_file>');

            markerFileName = f(i0 + 13 : i1 - 1);

            markerFileName = strrep(markerFileName, 'C:\Users\Biomechanik\Desktop\Willi\Study_LongitudinalMSK', 'C:\Users\Willi\ucloud\PhD\Study_LongitudinalMSK');

            [path, ~] = fileparts(markerFileName);
            if isfile(fullfile(path, 'EMG_filtered.sto'))
                emg_filtered = load_sto_file(fullfile(path, 'EMG_filtered.sto'));
                emg_fields = fieldnames(emg_filtered);
                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        if cycle.left.valid(j) || importTrialEvenIfReservesAreHigh
                            startTimeCycle = data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).time(1);
                            startTimeCycleMinusDelay = data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).time(1) - electromechanicalDelay;
                            endTimeCycle = data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).time(end);
                            emgIndizes = and(emg_filtered.time >= startTimeCycle, emg_filtered.time <= endTimeCycle);
                            emgIndizesWithDelay = and(emg_filtered.time >= startTimeCycleMinusDelay, emg_filtered.time <= endTimeCycle);
                            for i = 1 : numel(emg_fields)
                                if strcmp(emg_fields{i}(1:2), 'L_') || strcmp(emg_fields{i}(1:2), 'R_') || strcmp(emg_fields{i}, 'time') || contains(emg_fields{i}, 'EMG')
                                    data.EMG.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).(emg_fields{i}) = emg_filtered.(emg_fields{i})(emgIndizes);
                                    data.EMG_electromechanicalDelay.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).(emg_fields{i}) = emg_filtered.(emg_fields{i})(emgIndizesWithDelay);
                                    if startTimeCycleMinusDelay >= 0
                                        data.EMG_electromechanicalDelay.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).electromechanicalDelay = electromechanicalDelay;
                                    else
                                        data.EMG_electromechanicalDelay.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).electromechanicalDelay = electromechanicalDelay + startTimeCycleMinusDelay;
                                    end
                                    % one could include RMS calculations or
                                    % something similar here
                                end
                            end
                        end
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        if cycle.right.valid(j) || importTrialEvenIfReservesAreHigh
                            startTimeCycle = data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).time(1);
                            startTimeCycleMinusDelay = data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).time(1) - electromechanicalDelay;
                            endTimeCycle = data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).time(end);
                            emgIndizes = and(emg_filtered.time >= startTimeCycle, emg_filtered.time <= endTimeCycle);
                            emgIndizesWithDelay = and(emg_filtered.time >= startTimeCycleMinusDelay, emg_filtered.time <= endTimeCycle);
                            for i = 1 : numel(emg_fields)
                                if strcmp(emg_fields{i}(1:2), 'L_') || strcmp(emg_fields{i}(1:2), 'R_') || strcmp(emg_fields{i}, 'time') || contains(emg_fields{i}, 'EMG')
                                    data.EMG.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).(emg_fields{i}) = emg_filtered.(emg_fields{i})(emgIndizes);
                                    data.EMG_electromechanicalDelay.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).(emg_fields{i}) = emg_filtered.(emg_fields{i})(emgIndizesWithDelay);
                                    if startTimeCycleMinusDelay >= 0
                                        data.EMG_electromechanicalDelay.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).electromechanicalDelay = electromechanicalDelay;
                                    else
                                        data.EMG_electromechanicalDelay.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).electromechanicalDelay = electromechanicalDelay + startTimeCycleMinusDelay;
                                    end
                                    % one could include RMS calculations or
                                    % something similar here
                                end
                            end
                        end
                    end
                end
            end

            jrlFile = fullfile(currentFolder, 'Output', 'JRL', '_JointReaction_ReactionLoads.sto');
            if isfile(jrlFile)
                data.JRL.(modelList{p}).model_mass = model_mass;
                %crop data from left and right leg to stance phase
                tempData = load_sto_file(jrlFile);
                tempFieldNames = fieldnames(tempData);
                
                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        if cycle.left.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructLeft;
                            for i = 1 : numel(tempFieldNames)
                                tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                            end
                            data.JRL.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = tempStructLeft;
                            data.JRL.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds = double(cycle.left.end(j) - cycle.left.start(j)) / frequency;
                        end
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        if cycle.right.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructRight;
                            for i = 1 : numel(tempFieldNames)
                                tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                            end
                            data.JRL.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = tempStructRight;
                            data.JRL.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).durationInSeconds = double(cycle.right.end(j) - cycle.right.start(j)) / frequency;
                        end
                    end
                end
            elseif isfile(ikFile) % jrl file doesn't exist --> error in simulation --> add to ErrorScore 4
                ErrorScores(4,1) = ErrorScores(4,1)+1;
                data.ErrorScoreList.ErrorScore4{ErrorScores(4,1)} = currentFolder;
            end



            bodyKinFile = fullfile(currentFolder, 'Output', 'JRL', '_BodyKinematics_pos_global.sto');
            if isfile(bodyKinFile)
                %crop data from left and right leg to stance phase
                tempData = load_sto_file(bodyKinFile);
                tempFieldNames = fieldnames(tempData);
                
                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        if cycle.left.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructLeft;
                            for i = 1 : numel(tempFieldNames)
                                tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                            end

                            orientYfieldnames = tempFieldNames(contains(tempFieldNames, 'Oy'));
                            walkingDir = data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).walkingDirection;
                            walkingDirAngle = atan2d(walkingDir(3), walkingDir(1));
                            for i = 1 : numel(orientYfieldnames)
                                tempStructLeft.(orientYfieldnames{i}) = tempStructLeft.(orientYfieldnames{i}) + walkingDirAngle;
                            end

                            tempStructLeft = rmfield(tempStructLeft, 'time');
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = mergestructs(data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']), tempStructLeft);
                        end
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        if cycle.right.valid(j) || importTrialEvenIfReservesAreHigh
                            clear tempStructRight;
                            for i = 1 : numel(tempFieldNames)
                                tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                            end

                            orientYfieldnames = tempFieldNames(contains(tempFieldNames, 'Oy'));
                            walkingDir = data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).walkingDirection;
                            walkingDirAngle = atan2d(walkingDir(3), walkingDir(1));
                            for i = 1 : numel(orientYfieldnames)
                                tempStructRight.(orientYfieldnames{i}) = tempStructRight.(orientYfieldnames{i}) - walkingDirAngle;
                            end

                            tempStructRight = rmfield(tempStructRight, 'time');
                            data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = mergestructs(data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']), tempStructRight);
                        end
                    end
                end
            end

        else
            disp('reprocess this trial! cycle variable is unsorted!');
        end
    end
end




%% create table with #trials per Error Score & percentage of total #trials
RowDescription = {'number of trials', 'percentage'};
totalNumberOfTrials = [sum(ErrorScores); 100];
ErrorScore1(1,1) = ErrorScores(1,1);
ErrorScore2(1,1) = ErrorScores(2,1);
ErrorScore3(1,1) = ErrorScores(3,1);
ErrorScore4(1,1) = ErrorScores(4,1);
ErrorScore1(2,1) = ErrorScores(1,1)/totalNumberOfTrials(1,1) * 100;
ErrorScore2(2,1) = ErrorScores(2,1)/totalNumberOfTrials(1,1) * 100;
ErrorScore3(2,1) = ErrorScores(3,1)/totalNumberOfTrials(1,1) * 100;
ErrorScore4(2,1) = ErrorScores(4,1)/totalNumberOfTrials(1,1) * 100;
Scores_table = table(totalNumberOfTrials, ErrorScore1, ErrorScore2, ErrorScore3, ErrorScore4);

%% overwrite existing .mat file
delete([outputPath '\dataStruct_ErrorScores_no_trials_removed.mat']);
dataFile = [outputPath '\dataStruct_ErrorScores_no_trials_removed.mat'];
save(dataFile, 'data', 'Scores_table');
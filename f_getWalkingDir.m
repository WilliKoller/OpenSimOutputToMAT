function [outputArr_left, outputArr_right] = f_getArrayForField(modelStruct, fieldOfInterest, startOfStancePhase, footOffFrameStruct)
    if ~exist('startOfStancePhase', 'var')
        startOfStancePhase = 0;
    elseif ~exist('footOffFrameStruct', 'var')
        startOfStancePhase = 0;
        disp('pass struct with footoff frame numbers');
    end
    
    trialList = fieldnames(modelStruct);
    trialList = trialList(contains(trialList, 'T_'));
    leftTrials = trialList(contains(trialList, 'left'));
    rightTrials = trialList(contains(trialList, 'right'));

    outputArr_left = zeros(numel(leftTrials), 3);
    outputArr_right = zeros(numel(rightTrials), 3);
    
    for j = 1 : numel(leftTrials)
        if startOfStancePhase == 0
            outputArr_left(j, :)  = modelStruct.(leftTrials{j}).(fieldOfInterest);
        else
            outputArr_left(j, :)  = modelStruct.(leftTrials{j}).(fieldOfInterest)(startOfStancePhase : footOffFrameStruct.(leftTrials{j}).footOffFrame);
        end
    end
    for j = 1 : numel(rightTrials)
        if startOfStancePhase == 0
            outputArr_right(j, :)  = modelStruct.(rightTrials{j}).(fieldOfInterest);
        else
            outputArr_right(j, :)  = modelStruct.(rightTrials{j}).(fieldOfInterest)(startOfStancePhase : footOffFrameStruct.(rightTrials{j}).footOffFrame);
        end
    end
end


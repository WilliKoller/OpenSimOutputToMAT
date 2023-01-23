function data = removeFromStruct(data, trial2remove)
%removeFromStruct Summary of this function goes here
%   Detailed explanation goes here
    parts = strsplit(trial2remove, filesep);
    model = parts{end-1};
    trial2remove = parts{end};
    allTrials = fieldnames(data.IK.(model));
    
    for i = 1 : numel(allTrials)
       if contains(allTrials{i}, [trial2remove '_'])
            data.IK.(model) = rmfield(data.IK.(model), allTrials{i});
            if isfield(data.ID.(model), allTrials{i})
                data.ID.(model) = rmfield(data.ID.(model), allTrials{i});
            end
            if isfield(data.SO.(model), allTrials{i})
                data.SO.(model) = rmfield(data.SO.(model), allTrials{i});
            end
            if isfield(data.SO_Activation.(model), allTrials{i})
                data.SO_Activation.(model) = rmfield(data.SO_Activation.(model), allTrials{i});
            end
            if isfield(data.JRL.(model), allTrials{i})
                data.JRL.(model) = rmfield(data.JRL.(model), allTrials{i});
            end
       end
    end
end


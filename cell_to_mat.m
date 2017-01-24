function trialData = cell_to_mat(trialData)


% Convert cell arrays to double for those that need it
switch class(trialData)
    case 'dataset'
vn = fieldnames(trialData);
    case {'table', 'struct'}
vn = fieldnames(trialData);
end

if ismember('fixWindowEntered', vn)
    if iscell(trialData.fixWindowEntered)
        trialData.fixWindowEntered = cell2mat(trialData.fixWindowEntered);
    end
end
if ismember('stopSignalOn', vn)
    if iscell(trialData.stopSignalOn)
        trialData.stopSignalOn = cell2mat(trialData.stopSignalOn);
    end
end
if ismember('fixOn', vn)
    if iscell(trialData.fixOn)
        trialData.fixOn = cell2mat(trialData.fixOn);
    end
end
if ismember('targOn', vn)
    if iscell(trialData.targOn)
        trialData.targOn = cell2mat(trialData.targOn);
    end
end
if ismember('checkerOn', vn)
    if iscell(trialData.checkerOn)
        trialData.checkerOn = cell2mat(trialData.checkerOn);
    end
end
if ismember('responseCueOn', vn)
    if iscell(trialData.responseCueOn)
        trialData.responseCueOn = cell2mat(trialData.responseCueOn);
    end
end
if ismember('responseOnset', vn)
    if iscell(trialData.responseOnset)
        trialData.responseOnset = cell2mat(trialData.responseOnset);
    end
end
if ismember('targ1CheckerProp', vn)
    if iscell(trialData.targ1CheckerProp)
        trialData.targ1CheckerProp = cell2mat(trialData.targ1CheckerProp);
    end
end
if ismember('saccToTargIndex', vn)
    if iscell(trialData.saccToTargIndex)
        trialData.saccToTargIndex = cell2mat(trialData.saccToTargIndex);
    end
end
if ismember('targAngle', vn)
    if iscell(trialData.targAngle)
        trialData.targAngle = cell2mat(trialData.targAngle);
    end
end
if ismember('distAngle', vn)
    if iscell(trialData.distAngle)
        trialData.distAngle = cell2mat(trialData.distAngle);
    end
end
if ismember('targAmp', vn)
    if iscell(trialData.targAmp)
        trialData.targAmp = cell2mat(trialData.targAmp);
    end
end
if ismember('rewardOn', vn)
    if iscell(trialData.rewardOn)
        trialData.rewardOn = cell2mat(trialData.rewardOn(:,1));
    end
end
if ismember('trialDuration', vn)
    if iscell(trialData.trialDuration)
        trialData.trialDuration = cell2mat(trialData.trialDuration(:,1));
    end
end


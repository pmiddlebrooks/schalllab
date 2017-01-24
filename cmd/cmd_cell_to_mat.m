function trialData = cmd_cell_to_mat(trialData)


% Convert cell arrays to double for those that need it
if iscell(trialData.targOn)
    trialData.targOn = cell2mat(trialData.targOn);
end
if iscell(trialData.stopSignalOn)
    trialData.stopSignalOn = cell2mat(trialData.stopSignalOn);
end
if iscell(trialData.responseOnset)
    trialData.responseOnset = cell2mat(trialData.responseOnset);
end
if iscell(trialData.saccToTargIndex)
    trialData.saccToTargIndex = cell2mat(trialData.saccToTargIndex);
end
if iscell(trialData.targAngle)
    trialData.targAngle = cell2mat(trialData.targAngle);
end

function trialData = ccm_cell_to_mat(trialData)


% Convert cell arrays to double for those that need it
if iscell(trialData.stopOnset)
    trialData.stopOnset = cell2mat(trialData.stopOnset);
end
if iscell(trialData.responseCueOnset)
    trialData.responseCueOnset = cell2mat(trialData.responseCueOnset);
end
if iscell(trialData.responseOnset)
    trialData.responseOnset = cell2mat(trialData.responseOnset);
end
if iscell(trialData.target1CheckerProportion)
    trialData.target1CheckerProportion = cell2mat(trialData.target1CheckerProportion);
end
if iscell(trialData.saccadeToTargetIndex)
    trialData.saccadeToTargetIndex = cell2mat(trialData.saccadeToTargetIndex);
end
if iscell(trialData.targetAngle)
    trialData.targetAngle = cell2mat(trialData.targetAngle);
end

function newSSDIndex = staircase(lastStopOutcome, lastSSDIndex, nSSD)

newSSDIndex = lastSSDIndex;

maxStepSize = 2;
iStepSize = randi(maxStepSize);

if strcmp(lastStopOutcome, 'stopCorrect')
    newSSDIndex = min(iStepSize + lastSSDIndex, nSSD);
elseif isempty(lastStopOutcome) ||  strcmp(lastStopOutcome, 'stopIncorrectTarget') || strcmp(lastStopOutcome, 'stopIncorrectDistractor') || ...
        strcmp(lastStopOutcome, 'stopIncorrectPreSSDTarget') || strcmp(lastStopOutcome, 'stopIncorrectPreSSDDistractor') || ...
        strcmp(lastStopOutcome, 'targetHoldAbort') || strcmp(lastStopOutcome, 'distractorHoldAbort')
    newSSDIndex = max(lastSSDIndex - iStepSize, 1);
else
    % do nothing
end



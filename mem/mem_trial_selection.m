function  [trialList] = mem_trial_selection(trialData, outcomeArray, targetHemifield)

% function  [trialList] = mem_trial_selection(subjectID, sessionID, outcomeArray, targetHemifield)
%
% outcomeArray: array of strings indicating the outcomes to include:
%           {'all',
%           'saccToTarget', 
%           'targetHoldAbort', 'distractorHoldAbort',
%           'fixationAbort', 'saccadeAbort'
%
%
% targetHemifield: the location of the CORRECT TARGET 
%           'all', 'right', or 'left'.
%           By default right includes vertical up, left includes vertical
%           down.
%
%
%
%
%
%



nTrial = size(trialData.trialOutcome, 1);
trialLogical = ones(nTrial, 1);

trialData = cell_to_mat(trialData);



% Get list(s) of trials w.r.t. the outcome
if strcmp(outcomeArray, 'all')
    outcomeList = ones(nTrial, 1);
else
    outcomeList = zeros(nTrial, 1);
    for iOutcomeIndex = 1 : length(outcomeArray)
        iOutcome = outcomeArray{iOutcomeIndex};
        
        outcomeList = outcomeList + strcmp(trialData.trialOutcome, iOutcome);
    end
end
trialLogical = trialLogical & outcomeList;




% Get list(s) of trials w.r.t. target Angle Range (useful for 50%
% checkerboard displays wherein a random target is assigned
% Get list(s) of trials w.r.t. the SSDs
targTrial = zeros(nTrial, 1);
targAngle = trialData.targAngle;
if strcmp(targetHemifield, 'all')
    targTrial = ones(nTrial, 1);
elseif strcmp(targetHemifield, 'right')
    targTrial = ((targAngle > 270) & (targAngle <= 360)) | ...
        ((targAngle >= 0) & (targAngle < 90)) | ...
        ((targAngle > -90) & (targAngle < 0)) | ...
        ((targAngle >= -360) & (targAngle < -270));
    trialLogical = trialLogical & targTrial;
elseif strcmp(targetHemifield, 'left')
    targTrial = ((targAngle > 90) & (targAngle <= 270)) | ...
        ((targAngle < -90) & (targAngle > -270));
    trialLogical = trialLogical & targTrial;
end
% find(trialLogical)





trialList = find(trialLogical);

% toc
% disp('angle')













































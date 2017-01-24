function data = ccm_microsaccade_before_RT(subjectID, sessionID)
% data = ccm_microsaccade_before_RT(subjectID, sessionID, plotFlag, figureHandle)
%
% Finds microsaccades that occurred after the go cue but before the
% response
%

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
pSignalArray = ExtraVar.pSignalArray;
ssdArray = ExtraVar.ssdArray;

if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end
nTrial = size(trialData.trialOutcome,1);

microAmp = 2;  % degress below which is considered a microsaccade

% Make logicals of microsaccades for each trial
ms          = cellfun(@(x) x < microAmp, trialData.saccAmp, 'uniformoutput', false);
postCue     = cellfun(@(x,y) x > y, trialData.saccBegin, num2cell(trialData.checkerOn), 'uniformoutput', false);
preRT       = cellfun(@(x,y) x < y, trialData.saccBegin, num2cell(trialData.responseOnset), 'uniformoutput', false);
msValid     = cellfun(@(x,y,z) x & y & z, ms, postCue, preRT, 'uniformoutput', false);
nMS         = cellfun(@sum, msValid);

% Which trials are we interested in analyzing?
% Get default trial selection options
selectOpt           = ccm_trial_selection;
selectOpt.rightCheckerPct = 'collapse';
selectOpt.ssd       = 'none';

selectOpt.targDir  	= 'right';
selectOpt.outcome 	= {'goCorrectTarget', 'targetHoldAbort'};
goTargRightTrial    = ccm_trial_selection(trialData, selectOpt);
selectOpt.outcome 	= {'goCorrectDistractor', 'distractorHoldAbort'};
goDistRightTrial    = ccm_trial_selection(trialData, selectOpt);
selectOpt.targDir  	= 'left';
selectOpt.outcome  	= {'goCorrectTarget', 'targetHoldAbort'};
goTargLeftTrial     = ccm_trial_selection(trialData, selectOpt);
selectOpt.outcome 	= {'goCorrectDistractor', 'distractorHoldAbort'};
goDistLeftTrial     = ccm_trial_selection(trialData, selectOpt);


selectOpt.ssd       = 'collapse';

selectOpt.targDir  	= 'right';
selectOpt.outcome   = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
stopTargRightTrial  = ccm_trial_selection(trialData, selectOpt);
selectOpt.outcome  	=  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
stopDistRightTrial  = ccm_trial_selection(trialData, selectOpt);
selectOpt.targDir  	= 'left';
selectOpt.outcome  	= {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
stopTargLeftTrial   = ccm_trial_selection(trialData, selectOpt);
selectOpt.outcome  	=  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
stopDistLeftTrial   = ccm_trial_selection(trialData, selectOpt);



msTrial = nMS > 0;



goTrial = zeros(nTrial,1);
goTrial([goTargRightTrial;goDistRightTrial;goTargLeftTrial;goDistLeftTrial]) = 1;

stopTrial = zeros(nTrial,1);
stopTrial([stopTargRightTrial;stopDistRightTrial;stopTargLeftTrial;stopDistLeftTrial]) = 1;

rightTrial = zeros(nTrial,1);
rightTrial([goTargRightTrial;goDistRightTrial;stopTargRightTrial;stopDistRightTrial]) = 1;

leftTrial = zeros(nTrial,1);
leftTrial([goTargLeftTrial;goDistLeftTrial;stopTargLeftTrial;stopDistLeftTrial]) = 1;


data.msTrial    = msTrial;
data.goTrial    = goTrial;
data.stopTrial  = stopTrial;
data.leftTrial 	= leftTrial;
data.rightTrial = rightTrial;

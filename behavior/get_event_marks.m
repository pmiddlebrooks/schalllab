function eventMarks = get_event_marks(trialData, trial_type, alignEvent, alignindex, trialList)

markEvent = [];
alignedEvent = [];
alignArray = {'fixationWindowEntered', 'targetOnset', 'stopOnset', 'responseOnset', 'rewardOnset'};

if strcmp(alignEvent, 'fixationWindowEntered')       % Fxiation obtained
    % Add target onset marks w.r.t fixation
    markEvent = cell2mat(trialData.targetOnset(trialList));
    alignedEvent = cell2mat(trialData.fixationWindowEntered(trialList));
elseif strcmp(alignEvent, 'targetOnset')       % Target onset
    % Add saccade onset marks w.r.t target onset
%     if strcmp(trial_type, 'GOCorrect') || strcmp(trial_type, 'NOGOWrong')
if ~isnan(cell2mat(trialData.responseToTargetIndex(trialList)))
    markEvent = cellfun(@(x, y) x(y), trialData.responseOnset(trialList), trialData.responseToTargetIndex(trialList));
    alignedEvent = cell2mat(trialData.targetOnset(trialList));
end
elseif strcmp(alignEvent, 'stopOnset')       % Stop signal onset
    % Add saccade onset marks w.r.t stop signal
    % Only if the trial is a stop trial
if ~isnan(cell2mat(trialData.responseToTargetIndex(trialList)))
    markEvent = cellfun(@(x, y) x(y), trialData.responseOnset(trialList), trialData.responseToTargetIndex(trialList));
    alignedEvent = cell2mat(trialData.stopOnset(trialList));
end
%     if strcmp(trial_type, 'NOGOWrong')
%         eventMarks = alignindex + (Sacc_of_interest(trialList, 1) - StopSignal_(trialList, 1));
%     end;
elseif strcmp(alignEvent, 'responseOnset')   % Saccade
    % Add target onset marks w.r.t saccade onset
    % Only if a saccade was made
if ~isnan(cell2mat(trialData.responseToTargetIndex(trialList)))
    markEvent = cell2mat(trialData.targetOnset(trialList));
    alignedEvent = cellfun(@(x, y) x(y), trialData.responseOnset(trialList), trialData.responseToTargetIndex(trialList));
end
%     if strcmp(trial_type, 'GOCorrect') || strcmp(trial_type, 'NOGOWrong')
%         % For now, just use Target onset- might want to use Stop signal
%         % onset too in the future
%         eventMarks = alignindex + (Target_(trialList, 1) - Sacc_of_interest(trialList, 1));
%     end;
    
   % For now, don't worry about event marks for aligning on reward
% elseif strcmp(align_event, 'Correct_')   % Correct- not sure when this really is
%     % Only if it was a correct trial
%     if strcmp(trial_type, 'GOCorrect') || strcmp(trial_type, 'NOGOCorrect')
%         event_marks = Correct_(trial_list, 1);
%         span = [-400 800];
%     end;
% elseif strcmp(align_event, 'Reward_')   % Correct- not sure when this really is
%     % Only if it was a correct trial
%     if strcmp(trial_type, 'GOCorrect') || strcmp(trial_type, 'NOGOCorrect')
%     event_marks = Reward_(trial_list, 1);
%     span = [-400 800];
%     end;
end;

eventMarks = alignindex + (markEvent - alignedEvent);
function event_marks = get_event_marks(matpath, trial_type, align_event, alignindex, trial_list)

load(matpath)
event_marks = [];


if strcmp(align_event, 'Fixate_')       % Fxiation obtained
    % Add target onset marks w.r.t fixation
    event_marks = alignindex + (Target_(trial_list, 1) - Fixate_(trial_list, 1));
elseif strcmp(align_event, 'Target_')       % Target onset
    % Add saccade onset marks w.r.t target onset
    if strcmp(trial_type, 'GOCorrect') || strcmp(trial_type, 'NOGOWrong')
    event_marks = alignindex + (Sacc_of_interest(trial_list, 1) - Target_(trial_list, 1));
    end;
elseif strcmp(align_event, 'StopSignal_')       % Stop signal onset
    % Add saccade onset marks w.r.t stop signal
    % Only if the trial is a stop trial
    if strcmp(trial_type, 'NOGOWrong')
        event_marks = alignindex + (Sacc_of_interest(trial_list, 1) - StopSignal_(trial_list, 1));
    end;
elseif strcmp(align_event, 'Sacc_of_interest')   % Saccade
    % Add target onset marks w.r.t saccade onset
    % Only if a saccade was made
    if strcmp(trial_type, 'GOCorrect') || strcmp(trial_type, 'NOGOWrong')
        % For now, just use Target onset- might want to use Stop signal
        % onset too in the future
        event_marks = alignindex + (Target_(trial_list, 1) - Sacc_of_interest(trial_list, 1));
    end;
    
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

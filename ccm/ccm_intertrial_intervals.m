function ccm_intertrial_intervals(subjectID, sessionArray)
%% Intertrial interval as a function of trial type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA AND SET VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjectID = 'broca';
% [sessionArray, subjectIDArray] = task_session_array(iSubject, 'ccm', 'behavior1');
% sessionArray = {'bp050n02'};
sessionArray = 'behavior1';
acrossSession = true;
DELETE_ABORTS = false;

goTargInt = nan(length(sessionArray), 1);
goDistInt = nan(length(sessionArray), 1);
stopStopInt = nan(length(sessionArray), 1);
stopTargInt = nan(length(sessionArray), 1);
stopDistInt = nan(length(sessionArray), 1);


switch sessionArray
    case 'concat'
        iFile = fullfile('local_data',subjectID,strcat(subjectID,'RT.mat'));
        load(iFile)  % Loads trialData into workspace (and SessionData)
    case {'behavior1', 'behavior2', 'neural1', 'neural2'}
        iFile = fullfile('local_data',subjectID,strcat(subjectID,'_',sessionArray,'.mat'));
        load(iFile)  % Loads trialData into workspace (and SessionData)
    otherwise
        % Load the data
        [trialData, SessionData, ExtraVar] = load_data(subjectID, sessionArray, [ccm_min_vars, 'sessionTag']);
        pSignalArray    = ExtraVar.pSignalArray;
        targAngleArray	= ExtraVar.targAngleArray;
        ssdArray        = ExtraVar.ssdArray;
        nSSD            = length(ssdArray);
        
        if ~strcmp(SessionData.taskID, 'ccm')
            fprintf('Not a chioce countermanding saccade session, try again\n')
            return
        end
end

% Treat data differently if analyzing across sessions vs collapsed
% sessions:
% excludeTrialTriplet: Find session switch trials so we don't process them as if a new
% session was a continuation of one big session
if acrossSession
    if ismember(sessionArray, {'concat','behavior1', 'behavior2', 'neural1', 'neural2'})
        nSession = max(trialData.sessionTag);
    else
        nSession = 1;
    end
    excludeTrialTriplet = [];
else
    nSession = 1;
    excludeTrialPair = find(diff(trialData.sessionTag) < 0);
    excludeTrialTriplet = [excludeTrialPair; excludeTrialPair-1]; % Exclude last 2 trials of a session as possible beginning trials in triplets
end


for i = 1 : length(sessionArray)
    if acrossSession && ismember(sessionArray, {'concat','behavior1', 'behavior2', 'neural1', 'neural2'})
        td = trialData(trialData.sessionTag == i, :);
    else
        td = trialData;
    end
    
    if DELETE_ABORTS
        selectOpt = ccm_trial_selection;
        selectOpt.outcome = {...
            'goCorrectTarget', 'goCorrectDistractor', ...
            'stopCorrect', ...
            'targetHoldAbort', 'distractorHoldAbort', ...
            'stopIncorrectTarget', 'stopIncorrectDistractor'};
        validTrial = ccm_trial_selection(td, selectOpt);
        td = td(validTrial,:);
    end
    
    
    nTrial = size(td, 1);
    
    % Calculate intertrial intervals for the various trial types:
    Opt = ccm_options;
    
    % Go Target
    Opt.outcome         = {'goCorrectTarget'};
    goTargTrial     	= ccm_trial_selection(td, Opt);
    goTargTrial(goTargTrial == nTrial) = [];
    goTargIntTone        = td.trialOnset(goTargTrial+1) - (td.trialOnset(goTargTrial) + td.toneOn(goTargTrial));
    goTargInt(i)        = nanmean(td.trialOnset(goTargTrial+1) - td.trialOnset(goTargTrial));
    
    % Go Distractor
    Opt.outcome         = {'goCorrectDistractor'};
    goDistTrial            = ccm_trial_selection(td, Opt);
    goDistTrial(goDistTrial == nTrial) = [];
    goDistIntTone          = td.trialOnset(goDistTrial+1) - (td.trialOnset(goDistTrial) + td.toneOn(goDistTrial));
    goDistInt(i)        = nanmean(td.trialOnset(goDistTrial+1) - td.trialOnset(goDistTrial));
    
    % Stop Correct
    Opt.outcome         = {'stopCorrect'};
    stopStopTrial             = ccm_trial_selection(td, Opt);
    stopStopTrial(stopStopTrial == nTrial) = [];
    stopStopIntTone          = td.trialOnset(stopStopTrial+1) - (td.trialOnset(stopStopTrial) + td.toneOn(stopStopTrial));
    stopStopInt(i)        = nanmean(td.trialOnset(stopStopTrial+1) - td.trialOnset(stopStopTrial));
    
    % Stop Target
    Opt.outcome         = {'stopIncorrectTarget'};
    stopTargTrial             = ccm_trial_selection(td, Opt);
    stopTargTrial(stopTargTrial == nTrial) = [];
    stopTargIntTone          = td.trialOnset(stopTargTrial+1) - (td.trialOnset(stopTargTrial) + td.toneOn(stopTargTrial));
    stopTargInt(i)          = nanmean(td.trialOnset(stopTargTrial+1) - td.trialOnset(stopTargTrial));
    
    % Stop Target
    Opt.outcome         = {'stopIncorrectDistractor'};
    stopDistTrial             = ccm_trial_selection(td, Opt);
    stopDistTrial(stopDistTrial == nTrial) = [];
    stopDistIntTone          = td.trialOnset(stopDistTrial+1) - (td.trialOnset(stopDistTrial) + td.toneOn(stopDistTrial));
    stopDistInt(i)          = nanmean(td.trialOnset(stopDistTrial+1) - td.trialOnset(stopDistTrial));
    
    
    
end

fprintf('Intertrial intervals (from end trial tone on trial to next trial start:\n')
fprintf('goTarg:\t\t%2.f\n', mean(goTargIntTone))
fprintf('stopStop:\t%2.f\n', mean(stopStopIntTone))
fprintf('stopTarg:\t%2.f\n', mean(stopTargIntTone))
fprintf('goDist:\t\t%2.f\n', mean(goDistIntTone))
% fprintf('stopDist:\t%2.f\n', mean(stopDistIntTone))

    figITI = figure(287);
    clf
    hold on;
    yData = [mean(goTargIntTone), mean(stopStopIntTone), mean(stopTargIntTone), mean(goDistIntTone)];
    bar(yData)
%     errorbar(1:4, yData, [sem(goTargIntTone) sem(stopStopIntTone) sem(stopTargIntTone) sem(goDistIntTone)], 'linestyle' , 'none', 'color', 'k', 'linewidth' , 2)
%     ylim([.7 .9])
    print(figITI, [local_figure_path, sprintf('%s_%s_ITI', subjectID, sessionArray)],'-dpdf', '-r300')

    
fprintf('\nIntertrial intervals (start of trial to next trial start:\n')
fprintf('goTarg:\t\t%2.f\n', mean(goTargInt))
fprintf('stopStop:\t%2.f\n', mean(stopStopInt))
fprintf('stopTarg:\t%2.f\n', mean(stopTargInt))
fprintf('goDist:\t\t%2.f\n', mean(goDistInt))
% fprintf('stopDist:\t%2.f\n', mean(stopDistInt))

    figDur = figure(288);
    clf
    hold on;
    yData = [mean(goTargInt), mean(stopStopInt), mean(stopTargInt), mean(goDistInt)];
    bar(yData)
%     errorbar(1:4, yData, [sem(goTargInt) sem(stopStopInt) sem(stopTargInt) sem(goDistInt)], 'linestyle' , 'none', 'color', 'k', 'linewidth' , 2)
%     ylim([.7 .9])
    print(figDur, [local_figure_path, sprintf('%s_%s_Duration', subjectID, sessionArray)],'-dpdf', '-r300')


%% Session-long RTs
figureHandle = 66;
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(1, 1, figureHandle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA AND SET VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iSubject = 'broca';
[sessionArray, subjectIDArray] = task_session_array(iSubject, 'ccm', 'behavior2');
for i = 1 : length(sessionArray)
    sessionID = sessionArray{i};
    
    [td, S] = load_data(iSubject,sessionID);
    [allRT, rtOutlierTrial] = truncate_rt(td.rt, 120, 1200, 3);
    
    clf
    ax = axes('units', 'centimeters', 'position', [xAxesPosition(1) yAxesPosition(1) axisWidth axisHeight]);
    plot(ax, allRT, '-k')
    xlim([1 length(allRT)])
    ylim([min(allRT)-20 max(allRT)+20])
    
    savePlot = 'y';
    %     input('save?', 's');
    if strcmp(savePlot, 'y')
        localFigurePath = local_figure_path;
        print(figureHandle,[localFigurePath, sprintf('session_rts_%s', sessionID)],'-dpdf', '-r300')
    end
    
    
    
end
%% Session-long RTs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA AND SET VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iSubject = 'broca';
[sessionArray, subjectIDArray] = task_session_array(iSubject, 'ccm', 'behavior1');
sessionArray = {'bp050n02'};


DELETE_ABORTS = true;

% Set up figures
figureHandle = 96;
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(1, 1, figureHandle);
clf
figureHandle = 97;
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(4, 1, figureHandle);
clf
ax1 = axes('units', 'centimeters', 'position', [xAxesPosition(1) yAxesPosition(1) axisWidth axisHeight]);
cla
hold(ax1, 'on')
ax2 = axes('units', 'centimeters', 'position', [xAxesPosition(1) yAxesPosition(2) axisWidth axisHeight]);
cla
hold(ax2, 'on')
ax3 = axes('units', 'centimeters', 'position', [xAxesPosition(1) yAxesPosition(3) axisWidth axisHeight]);
cla
hold(ax3, 'on')
ax4 = axes('units', 'centimeters', 'position', [xAxesPosition(1) yAxesPosition(4) axisWidth axisHeight]);
cla
hold(ax4, 'on')


for i = 1 : length(sessionArray)
    sessionID = sessionArray{i};
    
    [td, S] = load_data(iSubject,sessionID);
    
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
    
    
    [allRT, rtOutlierTrial] = truncate_rt(td.rt, 120, 1200, 3);
    nTrial = size(td, 1);
    
    TRIAL_LAG = 20;  % moving avergage lag (trials)
    EPOCH_DURATION = 400;
    
    % Initialize moving average vectors
    maRT = nan(nTrial, 1);  % response times
    maPStop = nan(nTrial, 1);  % percentage of stop trials
    maRewRate = nan(nTrial, 1);  % reward rate
    
    for iTrial = TRIAL_LAG : nTrial
        maRT(iTrial)        = nanmean(td.rt(iTrial-TRIAL_LAG+1 : iTrial));
        maPStop(iTrial)     = 100 * (1 - (sum(isnan(td.ssd(iTrial-TRIAL_LAG+1 : iTrial))) / TRIAL_LAG));
        maRewRate(iTrial)        = nansum(cell2mat(td.rewardDuration(iTrial-TRIAL_LAG+1 : iTrial))) / TRIAL_LAG;
    end
    
    
    set(ax1, 'xLim', [0 nTrial])
    set(ax2, 'xLim', [0 nTrial])
    set(ax3, 'xLim', [0 nTrial])
    set(ax4, 'xLim', [0 nTrial])
    plot(ax1, allRT, '-k')
    plot(ax2, maRT, '-b')
    plot(ax3, maPStop, '-r')
    plot(ax4, maRewRate, '-g')
    print(figureHandle,[local_figure_path, sprintf('session_moving_avg_%s', sessionID)],'-dpdf', '-r300')
    
    
    %     for iTrial = 1 : 20 : nTrial
    %         clf
    %         ax = axes('units', 'centimeters', 'position', [xAxesPosition(1) yAxesPosition(1) axisWidth axisHeight]);
    %         cla
    %         hold(ax, 'on')
    %
    %         trialSegment = iTrial : iTrial + EPOCH_DURATION;
    %         rtSegment = allRT(trialSegment);
    %
    %
    %         plot(ax, trialSegment, rtSegment, '-k')
    %         plot(ax, trialSegment, maRT(trialSegment), '-b')
    %         plot(ax, trialSegment, maPStop(trialSegment) + min(allRT), '-r')
    %         plot(ax, trialSegment, maRewRate(trialSegment) + min(allRT), '-g')
    %         xlim([trialSegment(1) trialSegment(end)])
    %         %     ylim([min(allRT)-20 max(allRT)+20])
    %         ylim([0 max(allRT)+20])
    %
    %         savePlot = 'n';
    %         %     input('save?', 's');
    %         if strcmp(savePlot, 'y')
    %             localFigurePath = local_figure_path;
    %             print(figureHandle,[localFigurePath, sprintf('session_rts_%s', sessionID)],'-dpdf', '-r300')
    %         end
    %         pause
    %     end % for iTrial
    
end

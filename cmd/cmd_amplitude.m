function [ssrt, ssrtRight, ssrtLeft, ssdArray, stopRespondProb, nStop] = cmd_amplitude(subjectID, sessionID, plotFlag, figureHandle)
%%
if nargin < 3
    plotFlag = 1;
end
if nargin < 4
    figureHandle = 7575;
end


stopRespondProb = [];
ssrt = [];

% Load the data
[dataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessionID);
% If the file hasn't already been copied to a local directory, do it now
if exist(localDataFile, 'file') ~= 2
    copyfile(dataFile, localDataPath)
end
load(localDataFile);


if ~strcmp(SessionData.taskID, 'cmd')
    fprintf('Not a countermanding session, try again\n')
    return
end


% Convert cells to doubles if necessary
trialData = cell_to_mat(trialData);
nTrial = size(trialData.trialOutcome, 1);
ampArray = unique(trialData.targAmp);
nAmp = length(ampArray);
allRT = trialData.responseOnset - trialData.targOnset;



% Need to do a little SSD value adjusting, due to ms difference and 1-frame
% differences in SSD values
ssdArrayRaw = trialData.stopSignalOn - trialData.responseCueOn;
ssdArray = unique(ssdArrayRaw);
ssdArray = ssdArray(~isnan(ssdArray));
if ~isempty(ssdArray) && DO_STOPS
    a = diff(ssdArray);
    ssdArray(a == 1) = ssdArray(a == 1) + 1;
    ssdArray = unique(ssdArray);
    b = [ssdArray(1); diff(ssdArray)];
    ssdArray(b < 18) = [];
end




if plotFlag
    nRow = nAmp + 1;
    nColumn = 3;
    screenOrSave = 'save';
    %     [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = big_figure(nRow, nColumn, figureHandle, screenOrSave);
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nRow, nColumn, 'portrait', figureHandle);
    stopColor = [1 0 0];
    inhColor = [0 0 0];
    ssdMargin = 20;
    
    % axes names
    axInhFnLeft = (1:nAmp);
    axInhFnBoth = axInhFnLeft(end) + (1:nAmp);
    axInhFnRight = axInhFnBoth(end) + (1:nAmp);
    axInhFnGrand = 545;
    
    for iAmpIndex = 1 : nAmp
        % left targets inhibition function
        ax(axInhFnLeft(iAmpIndex)) = axes('units', 'centimeters', 'position', [xAxesPosition(iAmpIndex, 1) yAxesPosition(iAmpIndex, 1) axisWidth axisHeight]);
        cla
        hold(ax(axInhFnLeft(iAmpIndex)), 'on')
        % both targets inhibition function
        ax(axInhFnBoth(iAmpIndex)) = axes('units', 'centimeters', 'position', [xAxesPosition(iAmpIndex, 2) yAxesPosition(iAmpIndex, 2) axisWidth axisHeight]);
        cla
        hold(ax(axInhFnBoth(iAmpIndex)), 'on')
        % right targets inhibition function
        ax(axInhFnRight(iAmpIndex)) = axes('units', 'centimeters', 'position', [xAxesPosition(iAmpIndex, 3) yAxesPosition(iAmpIndex, 3) axisWidth axisHeight]);
        cla
        hold(ax(axInhFnRight(iAmpIndex)), 'on')
    end
    % both targets inhibition function
    ax(axInhFnGrand) = axes('units', 'centimeters', 'position', [mean(xAxesPosition(nAmp+1, :)) yAxesPosition(nAmp+1, 1) axisWidth axisHeight]);
    cla
    hold(ax(axInhFnGrand), 'on')
end







% ***********************************************************************
% Inhibition Function:
% ***********************************************************************
for iAmpIndex = 1 : nAmp
    iAmp = ampArray(iAmpIndex);
    iAmpTrial = find(trialData.targAmp == iAmp);
    for jSSDIndex = 1 : length(ssdArray)
        %     tic
        jSSD = ssdArray(jSSDIndex);
        
        % stop correct trials
        stopStopTrialRight{iAmpIndex, jSSDIndex} = intersect(iAmpTrial, cmd_trial_selection(trialData,  {'stopCorrect'}, jSSD, 'right'));
        stopStopTrialLeft{iAmpIndex, jSSDIndex} = intersect(iAmpTrial, cmd_trial_selection(trialData,  {'stopCorrect'}, jSSD, 'left'));
        stopStopTrialBoth{iAmpIndex, jSSDIndex} = [stopStopTrialRight{iAmpIndex, jSSDIndex}; stopStopTrialLeft{iAmpIndex, jSSDIndex}];
        
        
        % stop incorrect trials
        stopTargetOutcome = {'targetHoldAbort', 'stopIncorrectTarget', 'stopIncorrectPreSSDTarget'};
        
        
        stopIncorrectTrialRight{iAmpIndex, jSSDIndex} = intersect(iAmpTrial, cmd_trial_selection(trialData,  stopTargetOutcome, jSSD, 'right'));
        stopIncorrectTrialLeft{iAmpIndex, jSSDIndex} = intersect(iAmpTrial, cmd_trial_selection(trialData,  stopTargetOutcome, jSSD, 'left'));
        stopIncorrectTrialBoth{iAmpIndex, jSSDIndex} = [stopIncorrectTrialRight{iAmpIndex, jSSDIndex}; stopIncorrectTrialLeft{iAmpIndex, jSSDIndex}];
        
        
        
        % Inhibition function data points:
        stopRespondProbRight(iAmpIndex, jSSDIndex) = length(stopIncorrectTrialRight{iAmpIndex, jSSDIndex}) / (length(stopStopTrialRight{iAmpIndex, jSSDIndex}) + length(stopIncorrectTrialRight{iAmpIndex, jSSDIndex}));
        stopRespondProbLeft(iAmpIndex, jSSDIndex) = length(stopIncorrectTrialLeft{iAmpIndex, jSSDIndex}) / (length(stopStopTrialLeft{iAmpIndex, jSSDIndex}) + length(stopIncorrectTrialLeft{iAmpIndex, jSSDIndex}));
        stopRespondProbBoth(iAmpIndex, jSSDIndex) = length(stopIncorrectTrialBoth{iAmpIndex, jSSDIndex}) / (length(stopStopTrialBoth{iAmpIndex, jSSDIndex}) + length(stopIncorrectTrialBoth{iAmpIndex, jSSDIndex}));
        
        nStopRight(iAmpIndex, jSSDIndex) = length(stopStopTrialRight{iAmpIndex, jSSDIndex}) + length(stopIncorrectTrialRight{iAmpIndex, jSSDIndex});
        nStopLeft(iAmpIndex, jSSDIndex) = length(stopStopTrialLeft{iAmpIndex, jSSDIndex}) + length(stopIncorrectTrialLeft{iAmpIndex, jSSDIndex});
        nStopBoth(iAmpIndex, jSSDIndex) = length(stopStopTrialBoth{iAmpIndex, jSSDIndex}) + length(stopIncorrectTrialBoth{iAmpIndex, jSSDIndex});
        
        %     toc
    end % iSSDIndex
    
    
    % Need go trials cumulative RTs for calculation of SSRT
    goOutcome = {'goCorrectTarget'};
    goTargetRight{iAmpIndex} = intersect(iAmpTrial, cmd_trial_selection(trialData, goOutcome, 'none', 'right'));
    goTargetLeft{iAmpIndex} = intersect(iAmpTrial, cmd_trial_selection(trialData, goOutcome, 'none', 'left'));
    goTargetBoth{iAmpIndex} = intersect(iAmpTrial, cmd_trial_selection(trialData, goOutcome, 'none', 'all'));
    
    
    goRTRight = allRT(goTargetRight{iAmpIndex});
    goRTLeft = allRT(goTargetLeft{iAmpIndex});
    goRTBoth = allRT(goTargetBoth{iAmpIndex});
    
    
    
    
    
    % Inhibition function data
    keepSSDR = ~isnan(stopRespondProbRight(iAmpIndex,:));
    stopRespondProbRight = stopRespondProbRight(iAmpIndex,keepSSDR);
    nStopRight = nStopRight(iAmpIndex,keepSSDR);
    ssdArrayRight = ssdArray(keepSSDR);
    
    keepSSDL = ~isnan(stopRespondProbLeft(iAmpIndex,:));
    stopRespondProbLeft = stopRespondProbLeft(iAmpIndex,keepSSDL);
    nStopLeft = nStopLeft(iAmpIndex,keepSSDL);
    ssdArrayLeft = ssdArray(keepSSDL);
    
    keepSSDB = ~isnan(stopRespondProbBoth(iAmpIndex,:));
    stopRespondProbBoth = stopRespondProbBoth(iAmpIndex,keepSSDB);
    nStopBoth = nStopBoth(iAmpIndex,keepSSDB);
    ssdArrayBoth = ssdArray(keepSSDB);
    
    
    
    
    
    
    % Inhibition functions
    
    % [fitParametersR, lowestSSE] = Weibull(ssdArrayRight, stopRespondProbRight, nStopRight);
    %         [fitParametersR, lowestSSE] = Weibull_fast(ssdArrayRight, stopRespondProbRight, nStopRight);
    [fitParametersR, lowestSSE] = Weibull(ssdArrayRight, stopRespondProbRight, nStopRight);
    timePointsR = ssdArrayRight(1) : ssdArrayRight(end);
    inhibitionFnRight = weibull_curve(fitParametersR, timePointsR);
    
    
    % [fitParametersL, lowestSSE] = Weibull(ssdArrayLeft, stopRespondProbLeft, nStopLeft);
    %         [fitParametersL, lowestSSE] = Weibull_fast(ssdArrayLeft, stopRespondProbLeft, nStopLeft);
    [fitParametersL, lowestSSE] = Weibull(ssdArrayLeft, stopRespondProbLeft, nStopLeft);
    timePointsL = ssdArrayLeft(1) : ssdArrayLeft(end);
    inhibitionFnLeft = weibull_curve(fitParametersL, timePointsL);
    
    
    % [fitParametersL, lowestSSE] = Weibull(ssdArrayLeft, stopRespondProbLeft, nStopLeft);
    %         [fitParametersL, lowestSSE] = Weibull_fast(ssdArrayLeft, stopRespondProbLeft, nStopLeft);
    [fitParametersB, lowestSSE] = Weibull(ssdArrayBoth, stopRespondProbBoth, nStopBoth);
    timePointsB = ssdArrayBoth(1) : ssdArrayBoth(end);
    inhibitionFnBoth = weibull_curve(fitParametersB, timePointsB);
    
    
    
    % SSRTs
    [ssrtRight both_SSRTs SSRT_c SSRT_r meanSSD] = get_ssrt(ssdArrayRight, stopRespondProbRight, nStopRight', goRTRight, fitParametersR);
    [ssrtLeft both_SSRTs SSRT_c SSRT_r meanSSD] = get_ssrt(ssdArrayLeft, stopRespondProbLeft, nStopLeft', goRTLeft, fitParametersL);
    [ssrtBoth both_SSRTs SSRT_c SSRT_r meanSSD] = get_ssrt(ssdArrayBoth, stopRespondProbBoth, nStopBoth', goRTBoth, fitParametersB);
    fprintf('SSRTs:  Left \tRight \tCollapsed:\n\t%.2f \t%.2f \t%.2f\n', ssrtLeft, ssrtRight, ssrtBoth)
    
    
    if plotFlag
        plot(ax(axInhFnRight(iAmpIndex)), timePointsR, inhibitionFnRight, 'color', inhColor, 'linewidth', 2)
        plot(ax(axInhFnRight(iAmpIndex)), ssdArrayRight, stopRespondProbRight, '.k', 'markersize', 25)
        plot(ax(axInhFnBoth(iAmpIndex)), timePointsB, inhibitionFnBoth, 'color', inhColor, 'linewidth', 2)
        plot(ax(axInhFnBoth(iAmpIndex)), ssdArrayBoth, stopRespondProbBoth, '.k', 'markersize', 25)
        plot(ax(axInhFnLeft(iAmpIndex)), timePointsL, inhibitionFnLeft, 'color', inhColor, 'linewidth', 2)
        plot(ax(axInhFnLeft(iAmpIndex)), ssdArrayLeft, stopRespondProbLeft, '.k', 'markersize', 25)
        
        
        
        
        % inhibition function
        xlim(ax(axInhFnRight(iAmpIndex)), [ssdArray(1)-ssdMargin ssdArray(end)+ssdMargin])
        set(ax(axInhFnRight(iAmpIndex)), 'xtick', ssdArray)
        set(ax(axInhFnRight(iAmpIndex)), 'xtickLabel', ssdArray)
        ylim(ax(axInhFnRight(iAmpIndex)), [0 1]);
        set(get(ax(axInhFnRight(iAmpIndex)), 'ylabel'), 'String', 'p(Respond | stop)')
        
        xlim(ax(axInhFnLeft(iAmpIndex)), [ssdArray(1)-ssdMargin ssdArray(end)+ssdMargin])
        set(ax(axInhFnLeft(iAmpIndex)), 'xtick', ssdArray)
        set(ax(axInhFnLeft(iAmpIndex)), 'xtickLabel', ssdArray)
        ylim(ax(axInhFnLeft(iAmpIndex)), [0 1]);
        set(get(ax(axInhFnLeft(iAmpIndex)), 'ylabel'), 'String', 'p(Respond | stop)')
        
        xlim(ax(axInhFnBoth(iAmpIndex)), [ssdArray(1)-ssdMargin ssdArray(end)+ssdMargin])
        set(ax(axInhFnBoth(iAmpIndex)), 'xtick', ssdArray)
        set(ax(axInhFnBoth(iAmpIndex)), 'xtickLabel', ssdArray)
        ylim(ax(axInhFnBoth(iAmpIndex)), [0 1]);
        set(get(ax(axInhFnBoth(iAmpIndex)), 'ylabel'), 'String', 'p(Respond | stop)')
        
        
    end % if plotflag
    
    
end





% Grand Inhibition behavior
for jSSDIndex = 1 : length(ssdArray)
    jSSD = ssdArray(jSSDIndex);
    
    % stop correct trials
    stopStopTrial{jSSDIndex} = cmd_trial_selection(trialData,  {'stopCorrect'}, jSSD, 'all');
    
    
    % stop incorrect trials
    stopTargetOutcome = {'targetHoldAbort', 'stopIncorrectTarget', 'stopIncorrectPreSSDTarget'};
    
    stopIncorrectTrial{jSSDIndex} = cmd_trial_selection(trialData,  stopTargetOutcome, jSSD, 'all');
    
    % Inhibition function data points:
    stopRespondProb(jSSDIndex) = length(stopIncorrectTrial{jSSDIndex}) / (length(stopStopTrial{jSSDIndex}) + length(stopIncorrectTrial{jSSDIndex}));
    
    nStop(jSSDIndex) = length(stopStopTrial{jSSDIndex}) + length(stopIncorrectTrial{jSSDIndex});
    
end % iSSDIndex
goTarget = cmd_trial_selection(trialData, goOutcome, 'none', 'all');
goRT    = allRT(goTarget);

keepSSD = ~isnan(stopRespondProb);
stopRespondProb = stopRespondProb(keepSSD);
nStop = nStop(keepSSD);
ssdArray = ssdArray(keepSSD);


% [fitParametersG, lowestSSE] = Weibull(ssdArray, stopRespondProb, nStop);
%         [fitParametersG, lowestSSE] = Weibull_fast(ssdArray, stopRespondProb, nStop);
[fitParametersG, lowestSSE] = Weibull(ssdArray, stopRespondProb, nStop);
timePointsG = ssdArray(1) : ssdArray(end);
inhibitionFn = weibull_curve(fitParametersG, timePointsG);

[ssrt both_SSRTs SSRT_c SSRT_r meanSSD] = get_ssrt(ssdArray, stopRespondProb, nStop', goRT, fitParametersG);


if plotFlag
    plot(ax(axInhFnGrand), timePointsG, inhibitionFn, 'color', inhColor, 'linewidth', 2)
    plot(ax(axInhFnGrand), ssdArray, stopRespondProb, '.k', 'markersize', 25)
    
    xlim(ax(axInhFnGrand), [ssdArray(1)-ssdMargin ssdArray(end)+ssdMargin])
    set(ax(axInhFnGrand), 'xtick', ssdArray)
    set(ax(axInhFnGrand), 'xtickLabel', ssdArray)
    ylim(ax(axInhFnGrand), [0 1]);
    set(get(ax(axInhFnGrand), 'ylabel'), 'String', 'p(Respond | stop)')
end







print(gcf, ['~/schalllab/tempfigures/',sessionID, '_Behavior'], '-dpdf')



function gng_rt_distribution(subjectID, sessionID, plotFlag)


% Load the data
[dataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessionID);
% If the file hasn't already been copied to a local directory, do it now
if exist(localDataFile, 'file') ~= 2
    copyfile(dataFile, localDataPath)
end
load(localDataFile);




% Convert cells to doubles if necessary
trialData = cell_to_mat(trialData);



if ~strcmp(SessionData.taskID, 'gng')
    fprintf('Not a go/no-go session, try again /n')
    return
end


% trialOutcome = {'goCorrectTarget', 'stopIncorrectTarget', 'targetHoldAbort'};
nTrial = size(trialData.trialOutcome, 1);
goCorrectColor = [0 0 0];
stopIncorrectColor = [1 0 0];
goProportionArray = unique(trialData.goCheckerProportion);
nBin = 20;



if plotFlag
    [figureHandle, cmWidth, cmHeight] = standard_figure;
    interSpaceX = .3;
    interSpaceY = 1.2;
    leftMargin = 1.5;
    bottomMargin = 19;
    topMargin = 2;
    axisWidth = (cmWidth - 2*leftMargin - interSpaceX) / 2;
    axisHeight = 7;
    
    plotY = cmHeight - topMargin - axisHeight;
    
    ylimArray = [];
    
    %     if errorFlag
    %         xError = fill([1:length(condx_avg) length(condx_avg):-1:1], [(condx_avg - condx_sem) condx_avg(end:-1:1) + condx_sem(end:-1:1)], [1 .4 .4]);
    %         set(xError, 'edgecolor', 'none');
    %         % alpha(xError, .3)
    %         yError = fill([1:length(condy_avg) length(condy_avg):-1:1], [(condy_avg - condy_sem) condy_avg(end:-1:1) + condy_sem(end:-1:1)], [.4 .4 1]);
    %         set(yError, 'edgecolor', 'none');
    %         % alpha(yError, .3)
    %     end
end


for iTarget = 1 : 2
    switch iTarget
        case 1
            angleRange = -89 : 89;
        case 2
            angleRange = [-269 : -91, 91 : 269];
    end
    if plotFlag
        plotX = cmWidth - leftMargin - iTarget * (axisWidth + interSpaceX);
        ax(iTarget) = axes('units', 'centimeters', 'position', [plotX plotY axisWidth axisHeight]);
    end
    
    
    
    
    % Get correct go RT distribution
    goCorrectOutcome =  'goCorrectTarget';
    oddData = find(isnan(trialData.saccToTargIndex) & strcmp(trialData.trialOutcome, goCorrectOutcome) & ismember(trialData.targAngle, angleRange));
    if oddData
        fprintf('%d trials to target %d are listed as %s but don''t have valid saccades to target:\n', length(oddData), iTarget, goCorrectOutcome)
        disp(oddData)
    end
    
    iGoCorrectTrialIndices = strcmp(trialData.trialOutcome, goCorrectOutcome) & ...
        ismember(trialData.targAngle, angleRange);
    iGoCorrectTrialIndices(oddData) = 0;
    
    if sum(iGoCorrectTrialIndices)
        goCorrectRT = trialData.responseOnset(iGoCorrectTrialIndices) - trialData.responseCueOn(iGoCorrectTrialIndices);
        
        timeStep = (max(goCorrectRT) - min(goCorrectRT)) / nBin;
        goCorrectRTBinValues = hist(goCorrectRT, nBin);
        distributionArea = sum(goCorrectRTBinValues * timeStep);
        goCorrectPDF = goCorrectRTBinValues / distributionArea;
        goCorrectBinCenters = min(goCorrectRT)+timeStep/2 : timeStep : max(goCorrectRT)-timeStep/2;
        
        fprintf('Target %d:\t %d GoCorrect trials\t median RT: %.2f\n', iTarget, sum(iGoCorrectTrialIndices), median(goCorrectRT))
    end
    
    
    
    
    % Get incorrect stop RT distribution
    %     goCheckerPercentRange = goProportionArray(goProportionArray < 0.5) * 100;  % No go trials are those with goChekcerProportions less than .5
    goCheckerPercentRange = 'all';  % No go trials are those with goChekcerProportions less than .5
    stopsTargetHold = gng_trial_selection(subjectID, sessionID,  {'stopIncorrectTarget'}, goCheckerPercentRange);
    stopsTargetAbort = gng_trial_selection(subjectID, sessionID,  {'targetHoldAbort'}, goCheckerPercentRange);
    stopIncorrectTargetTrial = union(stopsTargetHold, stopsTargetAbort);
    stopIndices = zeros(nTrial, 1);
    stopIndices(stopIncorrectTargetTrial) = 1;
    
    oddData = find(isnan(trialData.saccToTargIndex) & stopIndices & ismember(trialData.targAngle), angleRange);
    if oddData
        fprintf('%d trials to target %d are listed as %s but don''t have valid saccades to target:\n', length(oddData), iTarget, goCorrectOutcome)
        disp(oddData)
    end
    
    
    
    iStopIncorrectTrialIndices = stopIndices & ...
        ismember(trialData.targAngle, angleRange);
    iStopIncorrectTrialIndices(oddData) = 0;
    
    if sum(iStopIncorrectTrialIndices)
        stopIncorrectRT = trialData.responseOnset(iStopIncorrectTrialIndices) - trialData.responseCueOn(iStopIncorrectTrialIndices);
        
        timeStep = (max(stopIncorrectRT) - min(stopIncorrectRT)) / nBin;
        stopIncorrectRTBinValues = hist(stopIncorrectRT, nBin);
        distributionArea = sum(stopIncorrectRTBinValues * timeStep);
        stopIncorrectPDF = stopIncorrectRTBinValues / distributionArea;
        stopIncorrectBinCenters = min(stopIncorrectRT)+timeStep/2 : timeStep : max(stopIncorrectRT)-timeStep/2;
        
        fprintf('Target %d:\t %d StopIncorrect trials\t median RT: %.2f\n', iTarget, sum(iStopIncorrectTrialIndices), median(stopIncorrectRT))
    end
    
    
    if plotFlag
        hold on
        %             plot(ax(iTarget), goCorrectBinCenters, goCorrectPDF, 'color', goCorrectColor, 'linewidth', 2)
        %             plot(ax(iTarget), stopIncorrectBinCenters, stopIncorrectPDF, 'color', stopIncorrectColor, 'linewidth', 2)
        if sum(iGoCorrectTrialIndices)
            plot(ax(iTarget), goCorrectBinCenters, goCorrectRTBinValues, 'color', goCorrectColor, 'linewidth', 2)
        end
        if sum(iStopIncorrectTrialIndices)
            plot(ax(iTarget), stopIncorrectBinCenters, stopIncorrectRTBinValues, 'color', stopIncorrectColor, 'linewidth', 2)
        end
        yMax = ylim;
        ylimArray = [ylimArray yMax(2)];
%         pause
    end
    
    
end  % iTarget




% Go through again and set all y-axis limits the same: should be able to
% use linkaxes to do this, but couldn't figure it out yet.
%linkaxes(h, 'y')
if plotFlag
    for iTarget = 1 : 2
        if iTarget ~= 2
            set(ax(iTarget), 'yticklabel', [])
        end
        set(ax(iTarget), 'ylim', [0, max(ylimArray)])
        set(ax(iTarget), 'xlim', [0, 1000])
    end
end


delete(localDataFile);




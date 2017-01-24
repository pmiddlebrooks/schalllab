function ccm_rt_distribution(subjectID, sessionID, plotFlag)

useTwoColors = true;

% Load the data
[dataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessionID);
% If the file hasn't already been copied to a local directory, do it now
if exist(localDataFile, 'file') ~= 2
    copyfile(dataFile, localDataPath)
end
load(localDataFile);
trialData = cell_to_mat(trialData);


if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end






% trialOutcome = {'goCorrectTarget', 'stopIncorrectTarget', 'targetHoldAbort'};
nTrial = size(trialData.trialOutcome, 1);
goCorrectColor = [0 0 0];
stopIncorrectColor = [1 0 0];
target1ProportionArray = unique(trialData.targ1CheckerProp);
ssdArray = unique(trialData.stopSignalOn - trialData.responseCueOn);
ssdArray(isnan(ssdArray)) = [];
stopIncorrectRT = [];
goCorrectRT     = [];
nBin = 40;

if useTwoColors
    if length(pSignalArray) == 6
        pSignalArray([2 5]) = [];
    elseif length(pSignalArray) == 7
        pSignalArray([2 4 6]) = [];
    end
end


if plotFlag
    nPlotColumn = 2;
    % Plot as many rows as there are discriminatory levels (use the maximum
    % from the individual targets), plus an extra row for collapsed RTs
    nPlotRow = 1 + max(length(target1ProportionArray(target1ProportionArray >= .51)), length(target1ProportionArray(target1ProportionArray <= .49)));
    figureHandleDist = 1001;
    figureHandleCum = 1002;
    [axisWidthDist, axisHeightDist, xAxesPositionDist, yAxesPositionDist]   = standard_figure(nPlotRow, nPlotColumn, 'portrait', figureHandleDist);
    [axisWidthCum, axisHeightCum, xAxesPositionCum, yAxesPositionCum]       = standard_figure(nPlotRow, nPlotColumn, 'portrait', figureHandleCum);
    
    switch subjectID
        case 'Broca'
            minPlotRT = 100;
            maxPlotRT = 600;
        otherwise
            minPlotRT = 300;
            maxPlotRT = 1200;
    end
    ylimArrayDist = [];
    ylimArrayCum = [];
end


for iTarget = 1 : 2
    switch iTarget
        case 1
            iProportionRange = target1ProportionArray(target1ProportionArray >= .51);
            iProportionRange = flipud(iProportionRange); % Rearrage order from easiest to hardest
            iPlotColumn = 2;
            %             proportionRange = proportionRange(1)  % The hardest
            %             proportionRange = max(target1ProportionArray); % The easiest
        case 2
            iProportionRange = target1ProportionArray(target1ProportionArray <= .49);
            iPlotColumn = 1;
            %             proportionRange = proportionRange(end)  % The hardest
            %             proportionRange = min(target1ProportionArray); % The easiest
    end
    for jProp = 1 : length(iProportionRange)
        jProportion = iProportionRange(jProp);
        
        
        %         if plotFlag
        %             figure(figureHandleDist);
        %             ax(jProp, iTarget) = axes('units', 'centimeters', 'position', [xAxesPositionDist(jProp, iTarget) yAxesPositionDist(jProp, iTarget) axisWidthDist axisHeightDist]);
        %             hold(ax(jProp, iTarget))
        %             figure(figureHandleCum);
        %             ax(jProp, iTarget) = axes('units', 'centimeters', 'position', [xAxesPositionCum(jProp, iTarget) yAxesPositionCum(jProp, iTarget) axisWidthCum axisHeightCum]);
        %             hold(ax(jProp, iTarget))
        %         end
        
        
        
        % Get correct go RT distribution
        goCorrectOutcome =  'goCorrectTarget';
        oddData = find(isnan(trialData.saccToTargIndex) & strcmp(trialData.trialOutcome, goCorrectOutcome) & ismember(trialData.targ1CheckerProp, iProportionRange));
        if oddData
            fprintf('%d trials to target %d are listed as %s but don''t have valid saccades to target:\n', length(oddData), iTarget, goCorrectOutcome)
            disp(oddData)
        end
        
        ijGoCorrectTrialIndices = strcmp(trialData.trialOutcome, goCorrectOutcome) & ...
            ismember(trialData.targ1CheckerProp, jProportion);
        ijGoCorrectTrialIndices(oddData) = 0;
        % Get rid of any trials in which the response to target index was
        % greater than the available reponses from a given trial
        mismatch = (trialData.saccToTargIndex - cellfun(@length, trialData.saccBegin)) > 0;
        ijGoCorrectTrialIndices(mismatch) = 0;
        
        if sum(ijGoCorrectTrialIndices)
            ijGoCorrectRT = trialData.responseOnset(ijGoCorrectTrialIndices) - trialData.responseCueOn(ijGoCorrectTrialIndices);
            ijGoCorrectRT(ijGoCorrectRT < 35) = [];
            goCorrectRT = [goCorrectRT; ijGoCorrectRT];
            
            % Go RT Distribution
            timeStep = (max(ijGoCorrectRT) - min(ijGoCorrectRT)) / nBin;
            goCorrectRTBinValues = hist(ijGoCorrectRT, nBin);
            distributionArea = sum(goCorrectRTBinValues * timeStep);
            goCorrectPDF = goCorrectRTBinValues / distributionArea;
            goCorrectBinCenters = min(ijGoCorrectRT)+timeStep/2 : timeStep : max(ijGoCorrectRT)-timeStep/2;
            
            % Cumulative RT function:
            goRT = sort(ijGoCorrectRT);
            ijPropGoRT = nan(length(min(goRT) : max(goRT)), 1);
            iRTIndex = 1;
            for iRT = min(goRT) : max(goRT)
                ijPropGoRT(iRTIndex) = sum(goRT <= iRT) / length(goRT);
                iRTIndex = iRTIndex + 1;
            end
            
            
            
            fprintf('Target %d:\t %d GoCorrect trials\t median RT: %.2f \t mean RT: %.2f\n', iTarget, sum(ijGoCorrectTrialIndices), median(ijGoCorrectRT), mean(ijGoCorrectRT))
        end
        
        
        % Get incorrect stop RT distribution to the appropriate target (correct
        % target trials)
        ijStopIncorrectTargetTrial = ccm_trial_selection(trialData,  {'stopIncorrectTarget'; 'targetHoldAbort'}, jProportion * 100, 'all', 'all');
        ijStopIndices = zeros(nTrial, 1);
        ijStopIndices(ijStopIncorrectTargetTrial) = 1;
        
        oddData = find(isnan(trialData.saccToTargIndex) & ijStopIndices & ismember(trialData.targ1CheckerProp, iProportionRange));
        if oddData
            fprintf('%d trials to target %d are listed as Incorrect Stops but don''t have valid saccades to target:\n', length(oddData), iTarget)
            disp(oddData)
        end
        
        
        ijStopIncorrectTrialIndices = ijStopIndices & ...
            ismember(trialData.targ1CheckerProp, jProportion);
        ijStopIncorrectTrialIndices(oddData) = 0;
        
        if sum(ijStopIncorrectTrialIndices) > 1
            ijStopIncorrectRT = trialData.responseOnset(ijStopIncorrectTrialIndices) - trialData.responseCueOn(ijStopIncorrectTrialIndices);
            ijStopIncorrectRT(ijStopIncorrectRT < 35) = [];
            stopIncorrectRT = [stopIncorrectRT; ijStopIncorrectRT];
            
            % Stop RT Distribution
            timeStep = (max(ijStopIncorrectRT) - min(ijStopIncorrectRT)) / nBin;
            stopIncorrectRTBinValues = hist(ijStopIncorrectRT, nBin);
            distributionArea = sum(stopIncorrectRTBinValues * timeStep);
            stopIncorrectPDF = stopIncorrectRTBinValues / distributionArea;
            stopIncorrectBinCenters = min(ijStopIncorrectRT)+timeStep/2 : timeStep : max(ijStopIncorrectRT)-timeStep/2;
            
            % Cumulative RT function:
            stopRT = sort(stopIncorrectRT);
            ijPropStopRT = nan(length(min(stopRT) : max(stopRT)), 1);
            iRTIndex = 1;
            for iRT = min(stopRT) : max(stopRT)
                ijPropStopRT(iRTIndex) = sum(stopRT <= iRT) / length(stopRT);
                iRTIndex = iRTIndex + 1;
            end
            
            fprintf('Target %d:\t %d StopIncorrect trials\t median RT: %.2f \t mean RT: %.2f\n', iTarget, sum(ijStopIncorrectTrialIndices), median(ijStopIncorrectRT), mean(ijStopIncorrectRT))
        end
        
        
        if plotFlag
            % RT Distribution
            figure(figureHandleDist)
            ax(jProp, iPlotColumn) = axes('units', 'centimeters', 'position', [xAxesPositionDist(jProp, iPlotColumn) yAxesPositionDist(jProp, iPlotColumn) axisWidthDist axisHeightDist]);
            hold(ax(jProp, iPlotColumn), 'on')
            set(ax(jProp, iPlotColumn), 'xlim', [minPlotRT maxPlotRT])
            if sum(ijGoCorrectTrialIndices) > 1
                plot(ax(jProp, iPlotColumn), goCorrectBinCenters, goCorrectRTBinValues, 'color', goCorrectColor, 'linewidth', 2)
            end
            if sum(ijStopIncorrectTrialIndices) > 1
                plot(ax(jProp, iPlotColumn), stopIncorrectBinCenters, stopIncorrectRTBinValues, 'color', stopIncorrectColor, 'linewidth', 2)
            end
            yLimitDist = ylim;
            ylimArrayDist = [ylimArrayDist yLimitDist(2)];
            
            % Cumulative RT
            figure(figureHandleCum)
            ax(jProp, iPlotColumn) = axes('units', 'centimeters', 'position', [xAxesPositionCum(jProp, iPlotColumn) yAxesPositionCum(jProp, iPlotColumn) axisWidthCum axisHeightCum]);
            hold(ax(jProp, iPlotColumn), 'on')
            set(ax(jProp, iPlotColumn), 'xlim', [minPlotRT maxPlotRT])
            if sum(ijGoCorrectTrialIndices) > 1
                
                plot(ax(jProp, iPlotColumn), min(goRT):max(goRT), ijPropGoRT, 'color', goCorrectColor, 'linewidth', 2)
            end
            if sum(ijStopIncorrectTrialIndices) > 1
                plot(ax(jProp, iPlotColumn), min(stopRT):max(stopRT), ijPropStopRT, 'color', stopIncorrectColor, 'linewidth', 2)
            end
            yLimitCum = ylim;
            ylimArrayCum = [ylimArrayCum yLimitCum(2)];
        end
        
    end % jProp
end  % iTarget






% Go through again and set all y-axis limits the same: should be able to
% use linkaxes to do this, but couldn't figure it out yet.
%linkaxes(h, 'y')
if plotFlag
    figure(figureHandleDist)
    for iPlotColumn = 1 : 2
        for jProp = 1 : nPlotRow-1
            if ax(jProp, iPlotColumn)
                %                 disp('ikay')
                %                 max(ylimArrayDist)
                if iPlotColumn == 2
                    set(ax(jProp, iPlotColumn), 'yticklabel', [])
                end
                ylim(ax(jProp, iPlotColumn), [0 max(ylimArrayDist)])
%                                 set(ax(jProp, iPlotColumn), 'xlim', [minPlotRT maxPlotRT])
%                                 set(ax(jProp, iPlotColumn), 'xlim', [300 1000])
            end
        end
    end
    
    figure(figureHandleCum)
    for iPlotColumn = 1 : 2
        for jProp = 1 : nPlotRow-1
            if ax(jProp, iPlotColumn)
                if iPlotColumn == 2
                    set(ax(jProp, iPlotColumn), 'yticklabel', [])
                end
                ylim(ax(jProp, iPlotColumn), [0 max(ylimArrayCum)])
                %                 set(ax(jProp, iPlotColumn), 'xlim', [minPlotRT maxPlotRT])
            end
        end
    end
end







% Do it again for the RT distributions collapsed across targets.
timeStep = (max(goCorrectRT) - min(goCorrectRT)) / nBin;
goCorrectRTBinValues = hist(goCorrectRT, nBin);
distributionArea = sum(goCorrectRTBinValues * timeStep);
goCorrectPDF = goCorrectRTBinValues / distributionArea;
goCorrectBinCenters = min(goCorrectRT)+timeStep/2 : timeStep : max(goCorrectRT)-timeStep/2;

if sum(ijStopIncorrectTrialIndices)
    timeStep = (max(stopIncorrectRT) - min(stopIncorrectRT)) / nBin;
    stopIncorrectRTBinValues = hist(stopIncorrectRT, nBin);
    distributionArea = sum(stopIncorrectRTBinValues * timeStep);
    stopIncorrectPDF = stopIncorrectRTBinValues / distributionArea;
    stopIncorrectBinCenters = min(stopIncorrectRT)+timeStep/2 : timeStep : max(stopIncorrectRT)-timeStep/2;
end

% if plotFlag
%     plotY = cmHeight - topBottomMargin - axisHeight*2 - interSpaceY;
%     plotX = leftRightMargin;
%     ax(3) = axes('units', 'centimeters', 'position', [plotX plotY axisWidth axisHeight]);
%     hold(ax(3))
%     plot(ax(3), goCorrectBinCenters, goCorrectRTBinValues, 'color', goCorrectColor, 'linewidth', 2)
%     plot(ax(3), mean(goCorrectRT), 5, '.', 'color', goCorrectColor, 'markersize', 10)
%     if sum(ijStopIncorrectTrialIndices)
%         plot(ax(3), stopIncorrectBinCenters, stopIncorrectRTBinValues, 'color', stopIncorrectColor, 'linewidth', 2)
%         plot(ax(3), mean(stopIncorrectRT), 5, '.', 'color', stopIncorrectColor, 'markersize', 10)
%     end
%     set(ax(3), 'xlim', [100, 500])
%     %         pause
% end




% Cumulative RT functions:
goRT = sort(goCorrectRT);
iRTIndex = 1;
for i = min(goRT) : max(goRT)
    propGoRT(iRTIndex) = sum(goRT <= i) / length(goRT);
    iRTIndex = iRTIndex + 1;
end
if sum(ijStopIncorrectTrialIndices)
    stopRT = sort(stopIncorrectRT);
    iRTIndex = 1;
    for i = min(stopRT) : max(stopRT)
        propStopRT(iRTIndex) = sum(stopRT <= i) / length(stopRT);
        iRTIndex = iRTIndex + 1;
    end
    [h,p] = kstest2(goRT, stopRT);
    fprintf('Kolmogorov-Smirnov test: p = %.4f\n', p)
end




if plotFlag
    
    % RT Distribution
    figure(figureHandleDist);
    plotX = (xAxesPositionDist(jProp, 1) + xAxesPositionDist(jProp, 2)) / 2;
    ax(jProp, iPlotColumn) = axes('units', 'centimeters', 'position', [plotX yAxesPositionDist(1+length(iProportionRange), 1) axisWidthDist axisHeightDist]);
    hold(ax(jProp, iPlotColumn), 'on')
    set(ax(jProp, iPlotColumn), 'xlim', [minPlotRT maxPlotRT])
%     set(ax(jProp, iPlotColumn), 'xlim', [300 1000])
    plot(ax(jProp, iPlotColumn), goCorrectBinCenters, goCorrectRTBinValues, 'color', goCorrectColor, 'linewidth', 2)
    plot(ax(jProp, iPlotColumn), mean(goCorrectRT), 5, '.', 'color', goCorrectColor, 'markersize', 10)
    if sum(ijStopIncorrectTrialIndices)
        plot(ax(jProp, iPlotColumn), stopIncorrectBinCenters, stopIncorrectRTBinValues, 'color', stopIncorrectColor, 'linewidth', 2)
        plot(ax(jProp, iPlotColumn), mean(stopIncorrectRT), 5, '.', 'color', stopIncorrectColor, 'markersize', 10)
    end
    %     set(ax(3), 'xlim', [100, maxPlotRT])
    
    % Cumulative RT
    figure(figureHandleCum);
    ax(jProp, iPlotColumn) = axes('units', 'centimeters', 'position', [plotX yAxesPositionDist(1+length(iProportionRange), 1) axisWidthCum axisHeightCum]);
    hold(ax(jProp, iPlotColumn), 'on')
    set(ax(jProp, iPlotColumn), 'xlim', [minPlotRT maxPlotRT])
%         set(ax(jProp, iPlotColumn), 'xlim', [300 1000])

    
    %         plotY = cmHeight - topBottomMargin - axisHeight*2 - interSpaceY;
    %     plotX = leftRightMargin + axisWidth + interSpaceX;
    %     ax(4) = axes('units', 'centimeters', 'position', [plotX plotY axisWidth axisHeight]);
    %     hold(ax(4))
    plot(ax(jProp, iPlotColumn), min(goRT):max(goRT), propGoRT, 'color', goCorrectColor, 'linewidth', 2)
    if sum(ijStopIncorrectTrialIndices) > 1
        plot(ax(jProp, iPlotColumn), min(stopRT):max(stopRT), propStopRT, 'color', stopIncorrectColor, 'linewidth', 2)
    end
    %     set(ax(jProp, iPlotColumn), 'xlim', [100, maxPlotRT])
    %         pause
end
print(figure(figureHandleDist), ['~/schalllab/tempfigures/',sessionID, '_RT_Dist'], '-dpdf')
print(figure(figureHandleCum), ['~/schalllab/tempfigures/',sessionID, '_RT_Cum'], '-dpdf')
% delete(localDataFile);


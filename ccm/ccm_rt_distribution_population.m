function data = ccm_rt_distribution_population(subjectID, sessionSet, options)
%%
if nargin < 3
    options = ccm_options;
    options.plotFlag = 1;
    options.printPlot = 1;
    options.saveName = [];
end

task = 'ccm';
if iscell(sessionSet)
    % If user enters sessionSet, get rid of repeated sessions in case there
    % were neural recordings with multiple units from a single session
    sessionArray = unique(sessionSet);
    subjectIDArray = repmat({subjectID}, length(sessionArray), 1);
else
[sessionArray, subjectIDArray] = task_session_array(subjectID, task, sessionSet);
end
% ****************************************************************************************
% Population CDF
% ****************************************************************************************

switch lower(subjectID)
    case 'joule'
        [td, S, E] =load_data(subjectID, sessionArray{1}, ccm_min_vars);
        pSignalArray = E.pSignalArray;
    case 'human'
        pSignalArray = [.35 .42 .46 .5 .54 .58 .65];
    case 'broca'
        if iscell(sessionSet)
        [td, S, E] =load_data(subjectID, sessionArray{1}, ccm_min_vars);
        pSignalArray = E.pSignalArray;
        else
              switch sessionSet
                 case 'behavior'
                    pSignalArray = [.41 .45 .48 .5 .52 .55 .59];
                 case 'behavior2'
                    pSignalArray = [.43 .45 .47 .53 .55 .57];
                 case 'neural1'
                    pSignalArray = [.41 .44 .47 .53 .56 .59];
                 case 'neural2'
                    pSignalArray = [.42 .44 .46 .54 .56 .58];
                   otherwise
                       if length(pSignalArray) == 6
                           pSignalArray([2 5]) = [];
                       elseif length(pSignalArray) == 7
                           pSignalArray([2 4 6]) = [];
                       end
              end
        end
    case 'xena'
        pSignalArray = [.35 .42 .47 .5 .53 .58 .65];
end
% Remove 50% coherence condition
pSignalArray(pSignalArray == .5) = [];

if mod(length(pSignalArray), 2)
    nSignalPerSide = ceil(length(pSignalArray) / 2);
else
    nSignalPerSide = length(pSignalArray) / 2;
end

nSession = length(sessionArray);
nRow = 3;
nColumn = 2;
cumAx = 1;
distAx = 2;
figureHandle = 65;
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nRow, nColumn, 'portrait', figureHandle);
ax(cumAx) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
ax(distAx) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
hold(ax(cumAx), 'on')
hold(ax(distAx), 'on')
stopColor = [.5 .5 .5];
stopColor = [1 0 0];
goColor = [0 0 0];
switch lower(subjectID)
    case 'human'
        set(ax(cumAx), 'Xlim', [300 1000])
    case 'broca'
        set(ax(cumAx), 'Xlim', [100 1000])
        set(ax(distAx), 'xlim', [100 1000])
    case 'xena'
        set(ax(cumAx), 'Xlim', [150 550])
        set(ax(distAx), 'xlim', [150 600])
end


% Cumulative RT functions:
goTargRT = [];
goDistRT = [];
stopTargRT = [];
stopDistRT = [];

for iSession = 1 : nSession
    
    optChron = ccm_chronometric;
    optChron.plotFlag = 0;
    optChron.collapseTarg = true;
    
    iData = ccm_chronometric(subjectIDArray{iSession}, sessionArray{iSession}, optChron);
    
    for s = 1 : nSignalPerSide
        
        
        goTargRT = [goTargRT; iData.goLeftToTarg{s}; iData.goRightToTarg{s}];
        goDistRT = [goDistRT; iData.goLeftToDist{s}; iData.goRightToDist{s}];
        stopTargRT = [stopTargRT; cell2mat(iData.stopRightToTarg(:,s)); cell2mat(iData.stopLeftToTarg(:,s))];
        stopDistRT = [stopDistRT; cell2mat(iData.stopRightToDist(:,s)); cell2mat(iData.stopLeftToDist(:,s))];
        
        % % Go data
        %     goLeftToTarg =  [goLeftToTarg; cellfun(@nanmean, iData.goLeftToTarg)];
        %     goRightToTarg = [goRightToTarg; cellfun(@nanmean, iData.goRightToTarg)];
        %     goLeftToDist =  [goLeftToDist; cellfun(@nanmean, iData.goLeftToDist)];
        %     goRightToDist = [goRightToDist; cellfun(@nanmean, iData.goRightToDist)];
        %
        %
        %     % Stop data
        %     iStopLeftToTarg = [nanmean(cell2mat(iData.stopLeftToTarg(:,1))), nanmean(cell2mat(iData.stopLeftToTarg(:,2))), nanmean(cell2mat(iData.stopLeftToTarg(:,3))), nanmean(cell2mat(iData.stopLeftToTarg(:,4)))];
        %     iStopRightToTarg = [nanmean(cell2mat(iData.stopRightToTarg(:,1))), nanmean(cell2mat(iData.stopRightToTarg(:,2))), nanmean(cell2mat(iData.stopRightToTarg(:,3))), nanmean(cell2mat(iData.stopRightToTarg(:,4)))];
        %     iStopLeftToDist = [nanmean(cell2mat(iData.stopLeftToDist(:,1))), nanmean(cell2mat(iData.stopLeftToDist(:,2))), nanmean(cell2mat(iData.stopLeftToDist(:,3))), nanmean(cell2mat(iData.stopLeftToDist(:,4)))];
        %     iStopRightToDist = [nanmean(cell2mat(iData.stopRightToDist(:,1))), nanmean(cell2mat(iData.stopRightToDist(:,2))), nanmean(cell2mat(iData.stopRightToDist(:,3))), nanmean(cell2mat(iData.stopRightToDist(:,4)))];
        %
        %     stopLeftToTarg = [stopLeftToTarg; iStopLeftToTarg];
        %     stopRightToTarg = [stopRightToTarg; iStopRightToTarg];
        %     stopLeftToDist = [stopLeftToDist; iStopLeftToDist];
        %     stopRightToDist = [stopRightToDist; iStopRightToDist];
    end
end

%     for i = 1 : nSession
%     for s = 1 : length(signalStrengthLeft)
%         goTargRT = [goTargRT; populationData.goRightToTarg{i}{s}; populationData.goLeftToTarg{i}{s}];
%         goDistRT = [goDistRT; populationData.goRightToDist{i}{s}; populationData.goLeftToDist{i}{s}];
%         stopTargRT = [stopTargRT; populationData.stopRightToTarg{i}{s}; populationData.stopLeftToTarg{i}{s}];
%         stopDistRT = [stopDistRT; populationData.stopRightToDist{i}{s}; populationData.stopLeftToDist{i}{s}];
%     end
% end
goTargRT = sort(goTargRT);
goDistRT = sort(goDistRT);
stopTargRT = sort(stopTargRT);
stopDistRT = sort(stopDistRT);

goTargRT(isnan(goTargRT)) = [];
goDistRT(isnan(goDistRT)) = [];
stopTargRT(isnan(stopTargRT)) = [];
stopDistRT(isnan(stopDistRT)) = [];

iRTIndex = 1;
% for i = min(goRT) : max(goRT)
for iMS = 1:1200
    propGoTargRT(iRTIndex) = sum(goTargRT <= iMS) / length(goTargRT);
    propGoDistRT(iRTIndex) = sum(goDistRT <= iMS) / length(goDistRT);
    propStopTargRT(iRTIndex) = sum(stopTargRT <= iMS) / length(stopTargRT);
    propStopDistRT(iRTIndex) = sum(stopDistRT <= iMS) / length(stopDistRT);
    iRTIndex = iRTIndex + 1;
end
% if sum(ijStopIncorrectTrialIndices)
%     stopRT = sort(stopIncorrectRT);
%     iRTIndex = 1;
%     for i = min(stopRT) : max(stopRT)
%         propStopTargRT(iRTIndex) = sum(stopTargRT <= i) / length(stopTargRT);
%         iRTIndex = iRTIndex + 1;
%     end
% end
%
% for i = 1 : nSession
%     cla
%     plot([populationData.ssdArray{i}(1): populationData.ssdArray{i}(end)], populationData.inhibitionFnGrand{i})
%     plot(populationData.ssdArray{i}, populationData.stopRespondProbGrand{i}, 'ok', 'markeredgecolor', 'k', 'markerfacecolor', 'k', 'markersize', 10)
%  pause
% end
[h,p] = kstest2(propGoTargRT, propGoDistRT);
fprintf('Kolmogorov-Smirnov test:  Go Correct vs Go Incorrect: p = %1.2d\n', p)
[h,p] = kstest2(propStopTargRT, propStopDistRT);
fprintf('Kolmogorov-Smirnov test:  Stop Correct vs Stop Incorrect: p = %1.2d\n', p)
[h,p] = kstest2(propGoTargRT, propStopTargRT);
fprintf('Kolmogorov-Smirnov test:  Go Correct vs Stop Correct: p = %1.2d\n', p)
[h,p] = kstest2(propGoDistRT, propStopDistRT);
fprintf('Kolmogorov-Smirnov test:  Go Incorrect vs Stop Incorrect: p = %1.2d\n', p)


box(ax(cumAx), 'off')
% set(ax(cumAx), 'xlim', [100 1000]);
plot(ax(cumAx), propGoTargRT, 'color', goColor, 'linewidth', 2)
plot(ax(cumAx), propGoDistRT, '--', 'color', goColor, 'linewidth', 2)

plot(ax(cumAx), propStopTargRT, 'color', stopColor, 'linewidth', 2)
plot(ax(cumAx), propStopDistRT, '--', 'color', stopColor, 'linewidth', 2)
legend(ax(cumAx), {'Go Target', 'Go Distractor', 'Stop Target', 'Stop Distractor'}, 'location', 'southeast');



nBin = 40;
goRT = [goTargRT; goDistRT];
% Go RT Distribution
timeStep = (max(goRT) - min(goRT)) / nBin;
goRTBinValues = hist(goRT, nBin);
distributionArea = sum(goRTBinValues * timeStep);
goCorrectPDF = goRTBinValues / distributionArea;
goCorrectBinCenters = min(goRT)+timeStep/2 : timeStep : max(goRT)-timeStep/2;

stopRT = [stopTargRT; stopDistRT];
% Go RT Distribution
timeStep = (max(stopRT) - min(stopRT)) / nBin;
stopRTBinValues = hist(stopRT, nBin);
distributionArea = sum(stopRTBinValues * timeStep);
stopStopPDF = stopRTBinValues / distributionArea;
stopStopBinCenters = min(stopRT)+timeStep/2 : timeStep : max(stopRT)-timeStep/2;

plot(ax(distAx), goCorrectBinCenters, goCorrectPDF, '-', 'color', goColor, 'linewidth', 2)
plot(ax(distAx), stopStopBinCenters, stopStopPDF, '-', 'color', stopColor, 'linewidth', 2)

print(figureHandle,fullfile(local_figure_path, subjectID,['ccm_population_rt_distribution_',options.saveName]),'-dpdf', '-r300')





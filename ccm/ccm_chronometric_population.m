function data = ccm_chronometric_population(subjectID, sessionSet, options)

if nargin < 3
    options = ccm_options;
    options.plotFlag = 1;
    options.printPlot = 1;
    options.saveName = [];
end

%%
% **************************************************************************************************
% Populaiton chronometric
% **************************************************************************************************
fprintf('\n\n\n\n')
disp('*******************************************************************************')
disp('Populaiton chronometric')

chronOpt = ccm_chronometric;
chronOpt.collapseTarg = true;
chronOpt.plotFlag = false;
chronOpt.USE_TWO_COLORS = false;

task = 'ccm';

if iscell(sessionSet)
    % If user enters sessionSet, get rid of repeated sessions in case there
    % were neural recordings with multiple units from a single session
    sessionArray = unique(sessionSet);
    subjectIDArray = repmat({subjectID}, length(sessionArray), 1);
else
    [sessionArray, subjectIDArray] = task_session_array(subjectID, task, sessionSet);
end

figureHandle = 4950;
% end



switch lower(subjectID)
    case 'joule'
        [~, ~, E] =load_data(subjectID, sessionArray{1}, ccm_min_vars);
        pSignalArray = E.pSignalArray;
    case 'human'
        pSignalArray = [.35 .42 .46 .5 .54 .58 .65];
    case 'broca'
[~, ~, E] = load_data(subjectID, sessionArray{end}, ccm_min_vars);
        pSignalArray = E.pSignalArray;
    case 'xena'
        pSignalArray = [.35 .42 .47 .5 .53 .58 .65];
end



nSession = length(sessionArray);



if options.plotFlag
    nRow = 3;
    nColumn = 2;
    rtAx = 1;
    ssdAx = 2;
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nRow, nColumn, 'portrait', figureHandle);
    ax(rtAx) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
    cla
    ax(ssdAx) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2) yAxesPosition(2, 2) axisWidth axisHeight]);
    cla
    %     stopColor = [.5 .5 .5];
    stopColor = [1 0 0];
    goColor = [0 0 0];
    choicePlotXMargin = .015;
    switch subjectID
        case 'Human'
            set(ax(rtAx), 'Ylim', [400 900])
        otherwise
            set(ax(rtAx), 'Ylim', [200 500])
%             set(ax(rtAx), 'Ylim', [250 550])
    end
    hold(ax(rtAx), 'on')
    hold(ax(ssdAx), 'on')
end



goLeftToTarg = [];
goRightToTarg = [];
goLeftToDist = [];
goRightToDist = [];
stopLeftToTarg = [];
stopRightToTarg = [];
stopLeftToDist = [];
stopRightToDist = [];
goTargSessionMean = [];
goDistSessionMean = [];
stopTargSessionMean = [];
stopDistSessionMean = [];

sessionData = [];
ssdData = [];
ssdRTData = [];
for iSession = 1 :  nSession
    
    display(sessionArray{iSession})
    iData = ccm_chronometric(subjectIDArray{iSession}, sessionArray{iSession}, chronOpt);
    % go and stop means for the session collapsed across signal strength
    goTargSessionMean = [goTargSessionMean; nanmean([cell2mat(iData.goLeftToTarg(:)); cell2mat(iData.goRightToTarg(:))])];
    goDistSessionMean = [goDistSessionMean; nanmean([cell2mat(iData.goLeftToDist(:)); cell2mat(iData.goRightToDist(:))])];
    stopTargSessionMean = [stopTargSessionMean; nanmean([cell2mat(iData.stopLeftToTarg(:)); cell2mat(iData.stopRightToTarg(:))])];
    stopDistSessionMean = [stopDistSessionMean; nanmean([cell2mat(iData.stopLeftToDist(:)); cell2mat(iData.stopRightToDist(:))])];
    
    
    
    % Go data
    goLeftToTarg =  [goLeftToTarg; cellfun(@nanmean, iData.goLeftToTarg)];
    goRightToTarg = [goRightToTarg; cellfun(@nanmean, iData.goRightToTarg)];
    goLeftToDist =  [goLeftToDist; cellfun(@nanmean, iData.goLeftToDist)];
    goRightToDist = [goRightToDist; cellfun(@nanmean, iData.goRightToDist)];
    
    
    
    iStopLeftToTarg = nan(1, size(iData.stopLeftToTarg, 2));
    iStopRightToTarg = nan(1, size(iData.stopLeftToTarg, 2));
    iStopLeftToDist = nan(1, size(iData.stopLeftToTarg, 2));
    iStopRightToDist = nan(1, size(iData.stopLeftToTarg, 2));
    for i = 1 : size(iData.stopLeftToTarg, 2)
        
        %     % Stop data
        iStopLeftToTarg(i) = nanmean(cell2mat(iData.stopLeftToTarg(:,i)));
        iStopRightToTarg(i) = nanmean(cell2mat(iData.stopRightToTarg(:,i)));
        iStopLeftToDist(i) = nanmean(cell2mat(iData.stopLeftToDist(:,i)));
        iStopRightToDist(i) = nanmean(cell2mat(iData.stopRightToDist(:,i)));
    end
    % Stop data
    %     iStopLeftToTarg = [nanmean(cell2mat(iData.stopLeftToTarg(:,1))), nanmean(cell2mat(iData.stopLeftToTarg(:,2))), nanmean(cell2mat(iData.stopLeftToTarg(:,3)))];
    %     iStopRightToTarg = [nanmean(cell2mat(iData.stopRightToTarg(:,1))), nanmean(cell2mat(iData.stopRightToTarg(:,2))), nanmean(cell2mat(iData.stopRightToTarg(:,3)))];
    %     iStopLeftToDist = [nanmean(cell2mat(iData.stopLeftToDist(:,1))), nanmean(cell2mat(iData.stopLeftToDist(:,2))), nanmean(cell2mat(iData.stopLeftToDist(:,3)))];
    %     iStopRightToDist = [nanmean(cell2mat(iData.stopRightToDist(:,1))), nanmean(cell2mat(iData.stopRightToDist(:,2))), nanmean(cell2mat(iData.stopRightToDist(:,3)))];
    
    stopLeftToTarg = [stopLeftToTarg; iStopLeftToTarg];
    stopRightToTarg = [stopRightToTarg; iStopRightToTarg];
    stopLeftToDist = [stopLeftToDist; iStopLeftToDist];
    stopRightToDist = [stopRightToDist; iStopRightToDist];
    
    
    % Mean RT within each SSD (collapsed across signal strength)
    for jSSD = 2 : length(iData.ssdArray)
        %     stopLeftTargRT = [stopLeftTargRT; nanmean(cell2mat(cellfun(@(x) x', iData.stopLeftToTarg(jSSD, :), 'uniformoutput', 0)))];
        %     stopLeftDistRT = [stopLeftDistRT; nanmean(cell2mat(cellfun(@(x) x', iData.stopLeftToDist(jSSD, :), 'uniformoutput', 0)))];
        %     stopRightTargRT = [stopRightTargRT; nanmean(cell2mat(cellfun(@(x) x', iData.stopRightToTarg(jSSD, :), 'uniformoutput', 0)))];
        %     stopRightDistRT = [stopRightDistRT; nanmean(cell2mat(cellfun(@(x) x', iData.stopRightToDist(jSSD, :), 'uniformoutput', 0)))];
        %     nanmean([...
        %         cell2mat(cellfun(@(x) x', iData.stopLeftToTarg(jSSD, :), 'uniformoutput', 0)),...
        %         cell2mat(cellfun(@(x) x', iData.stopLeftToDist(jSSD, :), 'uniformoutput', 0)),...
        %         cell2mat(cellfun(@(x) x', iData.stopRightToTarg(jSSD, :), 'uniformoutput', 0)),...
        %         cell2mat(cellfun(@(x) x', iData.stopRightToDist(jSSD, :), 'uniformoutput', 0))])
        jSSDRTAll = [cell2mat(cellfun(@(x) x', iData.stopLeftToTarg(jSSD, :), 'uniformoutput', 0)),...
            cell2mat(cellfun(@(x) x', iData.stopLeftToDist(jSSD, :), 'uniformoutput', 0)),...
            cell2mat(cellfun(@(x) x', iData.stopRightToTarg(jSSD, :), 'uniformoutput', 0)),...
            cell2mat(cellfun(@(x) x', iData.stopRightToDist(jSSD, :), 'uniformoutput', 0))];
        %         if jSSD < 2 &&
        jSSDRT = nanmean(jSSDRTAll);
        if ~isnan(jSSDRT)
            ssdRTData = [ssdRTData; jSSDRT];
            sessionData = [sessionData; iSession];
            ssdData = [ssdData; iData.ssdArray(jSSD)];
        end
    end
    
    if iSession == nSession
        pSignalArrayLeft = iData.pSignalArrayLeft;
        pSignalArrayRight = iData.pSignalArrayRight;
    end
    clear iData
    
end



% Mean RT within each SSD; Does RT increased with SSD?
[p, s] = polyfit(ssdData, ssdRTData, 1);
[y, delta] = polyval(p, ssdData, s);
stats = regstats(ssdData, ssdRTData)
R = corrcoef(ssdData, ssdRTData);
Rsqrd = R(1, 2)^2;
fprintf('\n\n R-squred: %.3f \t  Slope: %.3f \t p-value for regression: %.4f\n\n', Rsqrd, stats.beta(2), stats.tstat.pval(2))
cov(ssdData, ssdRTData);
xVal = min(ssdData) : .1 : max(ssdData);
yVal = p(1) * xVal + p(2);




goLeftTargMean = nanmean(goLeftToTarg, 1);
goLeftTargStd = nanstd(goLeftToTarg, 1);
goLeftTargSem = nanstd(goLeftToTarg, 1)/sqrt(size(goLeftToTarg, 1));
goRightDistMean = nanmean(goRightToDist, 1);
goRightDistStd = nanstd(goRightToDist, 1);
goRightDistSem = nanstd(goRightToDist, 1)/sqrt(size(goRightToDist, 1));

goRightTargMean = nanmean(goRightToTarg, 1);
goRightTargStd = nanstd(goRightToTarg, 1);
goRightTargSem = nanstd(goRightToTarg, 1)/sqrt(size(goRightToTarg, 1));
goLeftDistMean = nanmean(goLeftToDist, 1);
goLeftDistStd = nanstd(goLeftToDist, 1);
goLeftDistSem = nanstd(goLeftToDist, 1)/sqrt(size(goLeftToDist, 1));

stopLeftTargMean = nanmean(stopLeftToTarg, 1);
stopLeftTargStd = nanstd(stopLeftToTarg, 1);
stopLeftTargSem = nanstd(stopLeftToTarg, 1)/sqrt(size(stopLeftToTarg, 1));
stopRightDistMean = nanmean(stopRightToDist, 1);
stopRightDistStd = nanstd(stopRightToDist, 1);
stopRightDistSem = nanstd(stopRightToDist, 1)/sqrt(size(stopRightToDist, 1));

stopRightTargMean = nanmean(stopRightToTarg, 1);
stopRightTargStd = nanstd(stopRightToTarg, 1);
stopRightTargSem = nanstd(stopRightToTarg, 1)/sqrt(size(stopRightToTarg, 1));
stopLeftDistMean = nanmean(stopLeftToDist, 1);
stopLeftDistStd = nanstd(stopLeftToDist, 1);
stopLeftDistSem = nanstd(stopLeftToDist, 1)/sqrt(size(stopLeftToDist, 1));





goTargMean = [goLeftTargMean, goRightTargMean];
goTargStd = [goLeftTargStd, goRightTargStd];
goTargSem = [goLeftTargSem, goRightTargSem];
goTargPrint = [goTargMean; goTargStd; goTargSem];

goDistMean = [goRightDistMean, goLeftDistMean];
goDistStd = [goRightDistStd, goLeftDistStd];
goDistSem = [goRightDistSem, goLeftDistSem];
goDistPrint = [goDistMean; goDistStd; goDistSem];

stopTargMean = [stopLeftTargMean, stopRightTargMean];
stopTargStd = [stopLeftTargStd, stopRightTargStd];
stopTargSem = [stopLeftTargSem, stopRightTargSem];
stopTargPrint = [stopTargMean; stopTargStd; stopTargSem];

stopDistMean = [stopRightDistMean, stopLeftDistMean];
stopDistStd = [stopRightDistStd, stopLeftDistStd];
stopDistSem = [stopRightDistSem, stopLeftDistSem];
stopDistPrint = [stopDistMean; stopDistStd; stopDistSem];








if options.plotFlag
    % PLOT GO TRIALS
    plot(ax(rtAx), pSignalArrayLeft, goLeftTargMean, '-o', 'color', goColor, 'linewidth' , 2, 'markeredgecolor', goColor, 'markerfacecolor', [1 1 1], 'markersize', 10)
    errorbar(ax(rtAx), pSignalArrayLeft ,goLeftTargMean, goLeftTargStd, 'linestyle' , 'none', 'color', goColor, 'linewidth' , 2)
%     errorbar(ax(rtAx), pSignalArrayLeft ,goLeftTargMean, goLeftTargSem, 'linestyle' , 'none', 'color', goColor, 'linewidth' , 2)
    
    plot(ax(rtAx), pSignalArrayLeft, goRightDistMean, 'd', 'markeredgecolor', goColor, 'markersize', 10)
    % errorbar(ax(rtAx), pSignalArray ,goRightDistMean, goRightDistStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 2)
    
    plot(ax(rtAx), pSignalArrayRight, goRightTargMean, '-o', 'color', goColor, 'linewidth' , 2, 'markeredgecolor', goColor, 'markerfacecolor', [1 1 1], 'markersize', 10)
    errorbar(ax(rtAx), pSignalArrayRight ,goRightTargMean, goRightTargStd, 'linestyle' , 'none', 'color', goColor, 'linewidth' , 2)
%     errorbar(ax(rtAx), pSignalArrayRight ,goRightTargMean, goRightTargSem, 'linestyle' , 'none', 'color', goColor, 'linewidth' , 2)
    
    plot(ax(rtAx), pSignalArrayRight, goLeftDistMean, 'd', 'markeredgecolor', goColor, 'markersize', 10)
    % errorbar(ax(rtAx), pSignalArray ,goRightDistMean, goRightDistStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 2)
    
    
    % PLOT STOP TRIALS
    plot(ax(rtAx), pSignalArrayLeft, stopLeftTargMean, '-o', 'color', stopColor, 'linewidth' , 2, 'markeredgecolor', stopColor, 'markerfacecolor', [1 1 1], 'markersize', 10)
    errorbar(ax(rtAx), pSignalArrayLeft ,stopLeftTargMean, stopLeftTargStd, 'linestyle' , 'none', 'color', stopColor, 'linewidth' , 2)
%     errorbar(ax(rtAx), pSignalArrayLeft ,stopLeftTargMean, stopLeftTargSem, 'linestyle' , 'none', 'color', stopColor, 'linewidth' , 2)
    
    plot(ax(rtAx), pSignalArrayLeft, stopRightDistMean, 'd', 'markeredgecolor', stopColor, 'markersize', 10)
    % errorbar(ax(rtAx), pSignalArray ,stopRightDistMean, stopRightDistStd, '.' , 'linestyle' , 'none', 'color', stopColor, 'linewidth' , 2)
    
    plot(ax(rtAx), pSignalArrayRight, stopRightTargMean, '-o', 'color', stopColor, 'linewidth' , 2, 'markeredgecolor', stopColor, 'markerfacecolor', [1 1 1], 'markersize', 10)
    errorbar(ax(rtAx), pSignalArrayRight ,stopRightTargMean, stopRightTargStd, 'linestyle' , 'none', 'color', stopColor, 'linewidth' , 2)
%     errorbar(ax(rtAx), pSignalArrayRight ,stopRightTargMean, stopRightTargSem, 'linestyle' , 'none', 'color', stopColor, 'linewidth' , 2)
    
    plot(ax(rtAx), pSignalArrayRight, stopLeftDistMean, 'd', 'markeredgecolor', stopColor, 'markersize', 10)
    % errorbar(ax(rtAx), pSignalArray ,stopRightDistMean, stopRightDistStd, '.' , 'linestyle' , 'none', 'color', stopColor, 'linewidth' , 2)
    
    
    
    %     set(ax(ssrtAx), 'Xlim', [.38 .62])
    % set(ax(ssrtAx), 'Xlim', [.33 .67])
    set(ax(rtAx), 'Xlim', [pSignalArray(1) - choicePlotXMargin, pSignalArray(end) + choicePlotXMargin])
    plot(ax(rtAx), [.5 .5], ylim, '--k')
    set(ax(rtAx), 'xtick', pSignalArray)
    set(ax(rtAx), 'xtickLabel', pSignalArray*100)
    
    
    % PLOT RTs within each SSD
    plot(ax(ssdAx), ssdData, ssdRTData, '.', 'color', stopColor, 'markersize', 10)
    % regression line
    plot(ax(ssdAx), xVal, yVal, 'color', 'r', 'linewidth', 2)
    % 1:1 line
    get(ax(ssdAx), 'xlim')
    plot(ax(ssdAx), get(ax(ssdAx), 'ylim'), get(ax(ssdAx), 'ylim'), 'color', 'k', 'linewidth', 1)
    
    
    
end











% ANOVA calculations
% ==================================
rtData = [];
groupInh = [];
groupTarg = [];
groupSig = [];
goFlag = 1;
stopFlag = 2;
targFlag = 1;
distFlag = 2;
pSignalArray = [pSignalArrayLeft; pSignalArrayRight];
goTarg = [goLeftToTarg, goRightToTarg];
goDist = [goRightToDist, goLeftToDist];
stopTarg = [stopLeftToTarg, stopRightToTarg];
stopDist = [stopRightToDist, stopLeftToDist];
for i = 1 : size(goTarg, 2)
    rtData = [rtData; goTarg(:,i)];
    groupInh = [groupInh; repmat({'go'}, length(goTarg(:,i)), 1)];
    groupTarg = [groupTarg; repmat({'targ'}, length(goTarg(:,i)), 1)];
    groupSig = [groupSig; repmat(i, length(goTarg(:,i)), 1)];
end
for i = 1 : size(goDist, 2)
    rtData = [rtData; goDist(:,i)];
    groupInh = [groupInh; repmat({'go'}, length(goDist(:,i)), 1)];
    groupTarg = [groupTarg; repmat({'dist'}, length(goDist(:,i)), 1)];
    groupSig = [groupSig; repmat(i, length(goDist(:,i)), 1)];
end
for i = 1 : size(stopTarg, 2)
    rtData = [rtData; stopTarg(:,i)];
    groupInh = [groupInh; repmat({'stop'}, length(stopTarg(:,i)), 1)];
    groupTarg = [groupTarg; repmat({'targ'}, length(stopTarg(:,i)), 1)];
    groupSig = [groupSig; repmat(i, length(stopTarg(:,i)), 1)];
end
for i = 1 : size(stopDist, 2)
    rtData = [rtData; stopDist(:,i)];
    groupInh = [groupInh; repmat({'stop'}, length(stopDist(:,i)), 1)];
    groupTarg = [groupTarg; repmat({'dist'}, length(stopDist(:,i)), 1)];
    groupSig = [groupSig; repmat(i, length(stopDist(:,i)), 1)];
end
% size(rtData)
% size(groupInh)
% size(groupTarg)
% size(groupSig)
% [rtData, groupInh, groupTarg, groupSig]
fprintf('\n\n *******************  RT ANOVA  *******************  \n')
[p,table,stats] = anovan(rtData,{groupInh, groupTarg, groupSig}, 'varnames', {'Stop/Go', 'Targ/Dist', 'Signal'}, 'model', 'full', 'display', 'off');
fprintf('\nStop vs. Go: \t\tp = %.3f\nTarg vs Dist: \tp = %.3f\nSignal Strength: \tp = %.3f\n', p(1), p(2), p(3))
disp(table)

eta2Sig = table{4,2} / (table{4,2} + table{end-1,2})
eta2Targ = table{3,2} / (table{3,2} + table{end-1,2})
eta2InhTarg = table{5,2} / (table{2,2} + table{end-1,2})

fprintf('\n\n *********  RT Multicompare  *********  \n')
c = multcompare(stats, 'dimension', 3, 'display', 'off');
disp(c)
disp(stats)



% TTESTS for overall means calculations
% ==================================

[h,p, ci, stats] = ttest2(goTargSessionMean, stopTargSessionMean, [], [], 'unequal');
stats
rTest = sqrt(stats.tstat^2 / (stats.tstat^2 + stats.df))
fprintf('GoTarg: %.0f (%.0f) --vs-- stopTarg: %.0f (%.0f) \tt-test: t(%.2f) = %.2f \tp = %.5f \t CI: %.2f - %.2f\n', ...
    mean(goTargSessionMean), std(goTargSessionMean), ...
    mean(stopTargSessionMean), std(stopTargSessionMean), ...
    stats.df, stats.tstat, p, ci(1), ci(2));

[h,p, ci, stats] = ttest2(goDistSessionMean, stopDistSessionMean, [], [], 'unequal');
rTest = sqrt(stats.tstat^2 / (stats.tstat^2 + stats.df))
fprintf('GoDist: %.0f (%.0f) --vs-- stopDist: %.0f (%.0f) \tt-test: t(%.2f) = %.2f \tp = %.5f \t CI: %.2f - %.2f\n', ...
    mean(goDistSessionMean), std(goDistSessionMean), ...
    mean(stopDistSessionMean), std(stopDistSessionMean), ...
    stats.df, stats.tstat, p, ci(1), ci(2));









% Collect data in different format for SPSS Repeated measures ANOVA

rtDataSession = [];
groupSession = [];
session = [];
for iSession = 1 : nSession
    rtData = [];
    group = [];
    for i = 1 : size(goTarg, 2)
        rtData = [rtData, goTarg(iSession,i)];
        group = [group, {['goTarg',num2str(i)]}];
    end
    for i = 1 : size(goDist, 2)
        rtData = [rtData, goDist(iSession,i)];
        group = [group, {['goDist',num2str(i)]}];
    end
    for i = 1 : size(stopTarg, 2)
        rtData = [rtData, stopTarg(iSession,i)];
        group = [group, {['stopTarg',num2str(i)]}];
    end
    for i = 1 : size(stopDist, 2)
        rtData = [rtData, stopDist(iSession,i)];
        group = [group, {['stopDist',num2str(i)]}];
    end
    rtDataSession = [rtDataSession; rtData];
    groupSession = [groupSession; group];
end



%  ANOVA calculations
% ==================================
rtData = [];
groupInh = [];
groupTarg = [];
groupSig = [];
goFlag = 1;
stopFlag = 2;
targFlag = 1;
distFlag = 2;
pSignalArray = [pSignalArrayLeft; pSignalArrayRight];
goTarg = [goLeftToTarg, goRightToTarg];
goDist = [goRightToDist, goLeftToDist];
stopTarg = [stopLeftToTarg, stopRightToTarg];
stopDist = [stopRightToDist, stopLeftToDist];
sessionNumber = [];
for i = 1 : size(goTarg, 2)
    rtData = [rtData; goTarg(:,i)];
    groupInh = [groupInh; repmat(goFlag, length(goTarg(:,i)), 1)];
    groupTarg = [groupTarg; repmat(targFlag, length(goTarg(:,i)), 1)];
    groupSig = [groupSig; repmat(i, length(goTarg(:,i)), 1)];
    sessionNumber = [sessionNumber; (1:nSession)'];
    %     sessionNumber = [sessionNumber; i * ones(length(goTarg(:,i)), 1)];
end
for i = 1 : size(goDist, 2)
    rtData = [rtData; goDist(:,i)];
    groupInh = [groupInh; repmat(goFlag, length(goDist(:,i)), 1)];
    groupTarg = [groupTarg; repmat(distFlag, length(goDist(:,i)), 1)];
    groupSig = [groupSig; repmat(i, length(goDist(:,i)), 1)];
    sessionNumber = [sessionNumber; (1:nSession)'];
    %     sessionNumber = [sessionNumber; i * ones(length(goDist(:,i)), 1)];
end
for i = 1 : size(stopTarg, 2)
    rtData = [rtData; stopTarg(:,i)];
    groupInh = [groupInh; repmat(stopFlag, length(stopTarg(:,i)), 1)];
    groupTarg = [groupTarg; repmat(targFlag, length(stopTarg(:,i)), 1)];
    groupSig = [groupSig; repmat(i, length(stopTarg(:,i)), 1)];
    sessionNumber = [sessionNumber; (1:nSession)'];
    %     sessionNumber = [sessionNumber; i * ones(length(stopTarg(:,i)), 1)];
end
for i = 1 : size(stopDist, 2)
    rtData = [rtData; stopDist(:,i)];
    groupInh = [groupInh; repmat(stopFlag, length(stopDist(:,i)), 1)];
    groupTarg = [groupTarg; repmat(distFlag, length(stopDist(:,i)), 1)];
    groupSig = [groupSig; repmat(i, length(stopDist(:,i)), 1)];
    sessionNumber = [sessionNumber; (1:nSession)'];
    %     sessionNumber = [sessionNumber; i * ones(length(stopDist(:,i)), 1)];
end
size(groupSig)
size(sessionNumber)

datas = [rtData, groupInh, groupTarg, groupSig, sessionNumber];

if options.printPlot
    print(figureHandle,fullfile(local_figure_path, subjectID,['ccm_population_chronometric_',options.saveName]),'-dpdf', '-r300')
end




% fprintf('\nGo:   Targ: %.0f (%.0f) (%.0f) \tDist: %.0f (%.0f) (%.0f)', goTargMean, goTargSD, goTargSEM, goDistMean, goDistSD, goDistSEM);
% fprintf('\nStop: Targ: %.0f (%.0f) (%.0f) \tDist: %.0f (%.0f) (%.0f)', stopTargMean, stopTargSD, stopTargSEM, stopDistMean, stopDistSD, stopDistSEM);


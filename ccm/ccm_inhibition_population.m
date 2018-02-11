function data = ccm_inhibition_population(subjectID, sessionSet, options)

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

%% Population SSRT
%******************************************************************************
fprintf('\n\n\n\n')
disp('*******************************************************************************')
disp('Populaiton SSRT')

optInh              = ccm_inhibition;
optInh.plotFlag     = 0;
optInh.include50    = 0;
optInh.collapseTarg    = true;
optInh.USE_TWO_COLORS = false;



figureHandle = 4900;
plotSurface = false;

switch lower(subjectID)
    case 'joule'
        [td, S, E] =load_data(subjectID, sessionArray{1},ccm_min_vars);
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

if optInh.USE_TWO_COLORS
    if length(pSignalArray) == 6
        pSignalArray([2 5]) = [];
    elseif length(pSignalArray) == 7
        pSignalArray([2 4 6]) = [];
    end
end
if ~optInh.include50
    pSignalArray(pSignalArray == .5) = [];
end
nSignalStrength = length(pSignalArray);



nSession = length(sessionArray);

if options.plotFlag
    nRow = 3;
    nColumn = 2;
    axSSRT = 1;
    axInh = 2;
    axInhEach = 3;
    axRTSSD = 4;
    axInhSess = 5;
    axPred = 6;
    figureHandle = 9898;
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nRow, nColumn, 'portrait', figureHandle);
    
    ax(axSSRT) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
    hold(ax(axSSRT), 'on')
    ax(axInh) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
    hold(ax(axInh), 'on')
    ax(axInhEach) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 1) yAxesPosition(2, 1) axisWidth axisHeight]);
    hold(ax(axInhEach), 'on')
    ax(axRTSSD) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2) yAxesPosition(2, 2) axisWidth axisHeight]);
    hold(ax(axRTSSD), 'on')
    ax(axInhSess) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 1) yAxesPosition(3, 1) axisWidth axisHeight]);
    hold(ax(axInhSess), 'on')
    ax(axPred) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 2) yAxesPosition(3, 2) axisWidth axisHeight]);
    hold(ax(axPred), 'on')
    choicePlotXMargin = .015;
    switch subjectID
        case 'Human'
            set(ax(axSSRT), 'Ylim', [0 300])
        otherwise
            set(ax(axSSRT), 'Ylim', [0 150])
    end
end

nSession = length(sessionArray);
% pSignalArray = unique([populationData.pSignalArrayLeft{1}, populationData.pSignalArrayRight{1}])';







% GRAND INHIBITION FN *****************
stopRespondGrand = [];
ssdArrayGrand = [];

% SIGNAL STRENGTH INHIBITION FNS ******
dataSSD             = [];
dataSignalStrength  = [];
dataProbRespond     = [];
dataSession         = [];
dataRTSSD           = [];
rtGo             = nan(nSession, 1);
inhEach             = cell(1, nSignalStrength);
inhEachFit             = cell(1, nSignalStrength);
inhRTSSD            = cell(1, nSignalStrength);
inhRTSSDFit            = cell(1, nSignalStrength);
inhSess           = cell(1, nSignalStrength);
inhSessTime           = cell(1, nSignalStrength);
inhSessEach           = cell(1, nSignalStrength);

% SSRT ********************************
ssrtIntWeight = [];
ssrtMean = [];
ssrtInt = [];
ssrtIntSimple = [];


ssrtGrandIntWeight = [];
ssrtGrandMean = [];
ssrtGrandInt = [];

% OBSERVED VS. PREDICTED RTS  **********
dataStopRespondRT = [];
dataStopRespondRTPredict = [];

% Mean number of canceled stop trials in each session
nStopMean = nan(nSession, 1);
for iSession = 1 : nSession
    
    disp(sessionArray{iSession})
    iData               = ccm_inhibition(subjectIDArray{iSession}, sessionArray{iSession}, optInh);
    iData.ssrtIntegration = cellfun(@nanmean, iData.ssrtIntegration);
    
    % Mean number of stop trials in each session for the case when canceled stop trials number at least 10
    iOver10Trial = iData.nStopStop >= 10;
    nStopMean(iSession) = mean(iData.nStop(iOver10Trial));
    
    
    % GRAND INHIBITION FN ***********************
    ssdArrayGrand       = [ssdArrayGrand; iData.ssdArray];
    stopRespondGrand    = [stopRespondGrand; iData.stopRespondProbGrand];
    rtGo(iSession)     	= nanmean(cell2mat(iData.goTotalRT(:,1)));
    
    % SIGNAL STRENGTH INHIBITION FNS *************
    for j = 1 : nSignalStrength
        dataSSD             = [dataSSD; iData.ssd{j}];
        dataSignalStrength  = [dataSignalStrength; repmat(pSignalArray(j), length(iData.ssd{j}), 1)];
        iProbResp         = iData.stopRespondProb(j,:)';
        iProbResp(isnan(iProbResp)) = [];
        dataProbRespond  	= [dataProbRespond; iProbResp];
        dataSession        	= [dataSession; repmat(iSession, length(iData.ssd{j}), 1)];
        dataRTSSD           = [dataRTSSD; iData.goRTMinusSSD{j}];
        
        
        dataStopRespondRT       = [dataStopRespondRT; iData.stopRespondRT(j,:)'];
        dataStopRespondRTPredict = [dataStopRespondRTPredict; iData.stopRespondRTPredict(j,:)'];
        
        % Individual session inhibition functions
        inhSess{iSession}    = iData.inhibitionFnGrand;
        inhSessTime{iSession} = min(cell2mat(iData.ssd)) : max(cell2mat(iData.ssd));
        inhSessTime{iSession} = min(iData.ssdArray) : max(iData.ssdArray);
        inhSessEach{iSession} = iData.inhibitionFn;
    end
    
    
    
    % SSRT ********************************
    
    ssrtInt             = [ssrtInt; iData.ssrtIntegration'];
    ssrtIntWeight       = [ssrtIntWeight; iData.ssrtIntegrationWeighted'];
    ssrtIntSimple       = [ssrtIntSimple; iData.ssrtIntegrationSimple'];
    ssrtMean            = [ssrtMean; iData.ssrtMean'];
    
    ssrtGrandIntWeight 	= [ssrtGrandIntWeight; iData.ssrtCollapseIntegrationWeighted];
    ssrtGrandMean       = [ssrtGrandMean; iData.ssrtCollapseMean];
    ssrtGrandInt        = [ssrtGrandInt; nanmean(iData.ssrtCollapseIntegration)];
    
    clear iData
    
end
% Need to do a little SSD value adjusting, due to ms difference and 1-frame
% differences in SSD values
ssdArray = unique(dataSSD(~isnan(dataSSD)));
a = diff(ssdArray);
for i = 1 : 5
    changeSSD = ssdArray(a == i);
    for j = 1 : length(changeSSD)
        dataSSD(dataSSD == changeSSD(j)) = dataSSD(dataSSD == changeSSD(j)) + i;
    end
end
ssdArray = unique(dataSSD);



% GRAND INHIBITION FN  *****************************************
[fitParameters, lowestSSE] = Weibull(ssdArrayGrand, stopRespondGrand);
timePoints = min(ssdArrayGrand) : max(ssdArrayGrand);
inhPop = weibull_curve(fitParameters, timePoints);





% SIGNAL STRENGTH INHIBITION FNS  ********************************
for j = 1 : nSignalStrength
    jProbRespond = dataProbRespond(dataSignalStrength == pSignalArray(j));
    inhEach{j} = jProbRespond;
    
    % SSD vs P(Respond)
    jSSD = dataSSD(dataSignalStrength == pSignalArray(j));
    [fitParameters, lowestSSE] = Weibull(jSSD, jProbRespond);
    inhEachFit{j} = weibull_curve(fitParameters, timePoints);
    
    % RT-SSD vs P(Respond)
    jRTSSD = dataRTSSD(dataSignalStrength == pSignalArray(j));
    rtSSD{j} = jRTSSD;
    offsetVal = 1000;
    goSSDTimepoints = offsetVal + (min(-dataRTSSD) : max(-dataRTSSD));
    [fitParameters, lowestSSE] = Weibull(-jRTSSD + offsetVal, jProbRespond);
    inhRTSSDFit{j} = weibull_curve(fitParameters, goSSDTimepoints);
end




% SSRT   ***********************************************************
ssrtAvg = mean(ssrtIntWeight, 1);
ssrtSTD = std(ssrtIntWeight, 1);
ssrtSEM = std(ssrtIntWeight, 1) / sqrt(nSession);

ssrtMeanAvg = mean(ssrtMean, 1);
ssrtIntegrationAvg = mean(ssrtInt, 1);
ssrtIntegrationStd = std(ssrtInt, 1);
ssrtIntegrationSem = ssrtIntegrationStd ./ nSession;

% ssrtPlot = ssrtIntWeight;
ssrtPlot = ssrtIntWeight;
ssrtPlotAvg = mean(ssrtPlot, 1);
ssrtPlotStd = std(ssrtPlot, 1);
ssrtPlotSem = ssrtPlotStd ./ nSession;


% OBSERVED VS. PREDICTED NON-CANCELED RTS  ********************************
% ssdArray = unique(dataSSD)
stopRespondRT = nan(nSignalStrength, length(ssdArray));
stopRespondRTPredict = nan(nSignalStrength, length(ssdArray));
% for i = 1 : round(nSignalStrength/2)
for i = 1 : nSignalStrength
    for j = 1 : length(ssdArray)
        %         dataPoints = (dataSignalStrength == pSignalArray(i) | dataSignalStrength == 1-pSignalArray(i)) & dataSSD == ssdArray(j);
        dataPoints = dataSignalStrength == pSignalArray(i) & dataSSD == ssdArray(j);
        %         nansum(dat`aPoints)
        %         pause
        if nansum(dataPoints) > 2
            stopRespondRT(i,j) = nanmean(dataStopRespondRT(dataPoints));
            stopRespondRTPredict(i,j) = nanmean(dataStopRespondRTPredict(dataPoints));
            dp(i,j) = nansum(dataPoints);
        end
    end
end
stopRespondRTPredict(isnan(stopRespondRT)) = nan;
% Get rid of columns that are all NaN
nanCol = zeros(length(ssdArray), 1);
for j = 1 : length(ssdArray)
    if sum(isnan(stopRespondRT(:,j))) == length(pSignalArray)
        nanCol(j) = j;
    end
end
stopRespondRT(:,find(nanCol)) = [];
stopRespondRTPredict(:,find(nanCol)) = [];
predictSSDs = ssdArray;
predictSSDs(find(nanCol)) = [];




flipPropArray = pSignalArray;
flipPropArray(flipPropArray > .5) = abs(1 - flipPropArray(flipPropArray > .5));

pSignalArrayData = repmat(flipPropArray, size(ssrtPlot, 1), 1);
[p, s] = polyfit(pSignalArrayData(:), ssrtPlot(:), 1);
[y, delta] = polyval(p, pSignalArrayData(:), s);
stats = regstats(pSignalArrayData(:), ssrtPlot(:))
fprintf('p-value for regression: %.4f\n', stats.tstat.pval(2))
R = corrcoef(pSignalArrayData(:), ssrtPlot(:));
Rsqrd = R(1, 2)^2;
cov(pSignalArrayData(:), ssrtPlot(:));
xVal = min(pSignalArrayData(:)) : .001 : max(pSignalArrayData(:));
yVal = p(1) * xVal + p(2);
regressionLineX = min(flipPropArray) : .001 : abs(1 - min(flipPropArray));
regressionLineY = [yVal(1:end-1), fliplr(yVal)];


disp('********************************************')
disp('ANOVA for SSRT: Integration weighted')
group = num2cell(pSignalArray');
group = num2str(pSignalArray');
g = @(x) sprintf('%.2f', x);
group = cellfun(g, num2cell(pSignalArray'), 'uni', false);
% group = pSignalArray';
[p, table, stats] = anova1(ssrtIntWeight, group, 'off')
eta2Sig = table{2,2}/table{end,2}
c = multcompare(stats, 'display', 'off')
ssrtIntWtMean = mean(ssrtGrandIntWeight)
ssrtIntWtStd = std(ssrtGrandIntWeight)
data.stats = stats;
data.table = table;

disp('********************************************')
disp('ANOVA for SSRT: Mean')
[p, table, stats] = anova1(ssrtMean, group, 'off')
eta2Sig = table{2,2}/table{end,2}
c = multcompare(stats, 'display', 'off')
ssrtMeanGrandMean = mean(ssrtGrandMean)
ssrtMeanGrandStd = std(ssrtGrandMean)


disp(' ********************************************')
disp('ANOVA for SSRT: Integration')
[p, table, stats] = anova1(ssrtInt, group, 'off')
eta2Sig = table{2,2}/table{end,2}
c = multcompare(stats, 'display', 'off')
ssrtIntegrationGrandMean = mean(ssrtGrandInt)
ssrtIntegrationGrandStd = std(ssrtGrandInt)


disp('********************************************')
disp('ANOVA for Inhibitions functions: SSD')
% [p,table,stats] = anovan(dataSSD,{dataSignalStrength, dataProbRespond, dataSession}, 'varnames', {'Signal Strength', 'p(Respond)', 'Session'}, 'model', 'full', 'display', 'off')
% [p,table,stats] = anovan(dataProbRespond,{dataSignalStrength, dataSSD, dataSession}, 'varnames', {'Signal Strength', 'SSD', 'Session'}, 'display', 'off')
[p,table,stats] = anovan(dataProbRespond,{dataSignalStrength, dataSSD}, 'varnames', {'Signal Strength', 'SSD'}, 'display', 'off')
% c = multcompare(stats, 'dimension', 1, 'display', 'off')
eta2SSD = table{2,2} / (table{2,2} + table{end-1,2})


disp('********************************************')
disp('ANOVA for Inhibitions functions: RT - SSD')
% [p,table,stats] = anovan(dataRTSSD,{dataSignalStrength, dataProbRespond, dataSession}, 'varnames', {'Signal Strength', 'p(Respond)', 'Session'}, 'model', 'full', 'display', 'off')
% [p,table,stats] = anovan(dataProbRespond,{dataSignalStrength, dataRTSSD, dataSession}, 'continuous', [2], 'varnames', {'Signal Strength', 'RT-SSD', 'Session'}, 'display', 'off')
[p,table,stats] = anovan(dataProbRespond,{dataSignalStrength, dataRTSSD}, 'continuous', [2], 'varnames', {'Signal Strength', 'RT-SSD'}, 'display', 'off')
% c = multcompare(stats, 'dimension', 1, 'display', 'off')
eta2RTSSD = table{2,2} / (table{2,2} + table{end-1,2})









if options.plotFlag
    
    
    % GRAND INHIBITION FN ***********************
    plot(ax(axInh), ssdArrayGrand, stopRespondGrand, 'ok', 'markeredgecolor', 'k', 'markerfacecolor', 'k', 'markersize', 6)
    % plot(ax(inhAx), inhGrandMean, 'k', 'linewidth', 2)
    plot(ax(axInh), timePoints, inhPop, 'k', 'linewidth', 2)
    ylim([0 1])
    switch subjectID
        case 'Human'
            set(ax(axInh), 'Xlim', [0 max(ssdArrayGrand)])
            set(ax(axInhEach), 'Xlim', [0 max(ssdArrayGrand)])
        otherwise
            set(ax(axInh), 'Xlim', [0 max(ssdArrayGrand)])
            set(ax(axInhEach), 'Xlim', [0 max(ssdArrayGrand)])
    end
    
    
    
    
    
    % SSRT  ********************************
    %    plot(ax(axSSRT), pSignalArray, mean(ssrtInt, 1), '-ok', 'markeredgecolor', 'k', 'markerfacecolor', 'k', 'markersize', 10)
    errorbar(ax(axSSRT), pSignalArray ,ssrtPlotAvg, ssrtPlotStd, '.' , 'linestyle' , 'none', 'color', 'k', 'linewidth' , 3)
    plot(ax(axSSRT), pSignalArray, mean(ssrtPlot, 1), 'ok', 'markeredgecolor', 'k', 'markerfacecolor', 'k', 'markersize', 10)
    errorbar(ax(axSSRT), pSignalArray ,ssrtPlotAvg, ssrtPlotSem, '.' , 'linestyle' , 'none', 'color', 'g', 'linewidth' , 2)
    %    plot(ax(axSSRT), pSignalArray, mean(ssrtMean, 1), '-og', 'markeredgecolor', 'g', 'markerfacecolor', 'g', 'markersize', 10)
%      plot(ax(axSSRT), pSignalArray, mean(ssrtInt, 1), '-ob', 'markeredgecolor', 'r', 'markerfacecolor', 'r', 'markersize', 10)
   
    set(ax(axSSRT), 'Xlim', [pSignalArray(1) - choicePlotXMargin, pSignalArray(end) + choicePlotXMargin])
    set(ax(axSSRT), 'xtick', pSignalArray)
    set(ax(axSSRT), 'xtickLabel', pSignalArray*100)
    
    % Regression Line:
    %     plot(ax(axSSRT), regressionLineX, regressionLineY, 'color', 'r', 'linewidth', 2)
    
    
    
    
    
    
    % SIGNAL STRENGTH INHIBITION FNs: SSD ***********************
    minColorGun = .25;
    maxColorGun = 1;
    
    % Separate (earlier) loop for the data points, so they appear under the
    % lines
    for i = 1 : nSignalStrength
        iPct = pSignalArray(i) * 100;
        % Determine color to use for plot based on which checkerboard color
        % proportion being used. Normalize the available color spectrum to do
        % it
        if iPct == 50
            inhColor = [0 0 0];
        elseif iPct < 50
            colorNorm = .5 - pSignalArray(1);
            colorProp = (.5 - pSignalArray(i)) / colorNorm;
            colorGun = minColorGun + (maxColorGun - minColorGun) * colorProp;
            inhColor = [0 colorGun colorGun];
        elseif iPct > 50
            colorNorm = pSignalArray(end) - .5;
            colorProp = (pSignalArray(i) - .5) / colorNorm;
            colorGun = minColorGun + (maxColorGun - minColorGun) * colorProp;
            inhColor = [colorGun 0 colorGun];
        end
        
        
        iProbRespond = dataProbRespond(dataSignalStrength == pSignalArray(i));
        iSSD = dataSSD(dataSignalStrength == pSignalArray(i));
        iRTSSD = dataRTSSD(dataSignalStrength == pSignalArray(i));
        
        %         plot(ax(axRTSSD), -iRTSSD, iProbRespond, '.', 'color', inhColor, 'markersize', 10);
        %         plot(ax(axInhEach), iSSD, iProbRespond, '.', 'color', inhColor, 'markersize', 10)
    end
    
    for i = 1 : nSignalStrength
        iPct = pSignalArray(i) * 100;
        
        
        goSSDTimepoints = offsetVal + (min(-dataRTSSD) : max(-dataRTSSD));
        cMap = ccm_colormap(pSignalArray);
        inhColor = cMap(i,:);
        
        plot(ax(axInhEach), timePoints, inhEachFit{i}, 'color', inhColor, 'linewidth', 2)
        plot(ax(axRTSSD), (goSSDTimepoints - offsetVal), inhRTSSDFit{i}, '-', 'color', inhColor, 'linewidth', 2);
        %          plot(ax(axInhEach), timePoints, inhEach{i}, 'color', inhColor, 'linewidth', 2)
        %          plot(ax(axRTSSD), (goSSDTimepoints - offsetVal), inhRTSSD{i}, '-', 'color', inhColor, 'linewidth', 2);
        
        %          % greyscale version
        %          if iPct == 50
        %          inhColor = [0 0 0];
        %          plot(ax(axInhEach), timePoints, inhEach{i}, 'color', inhColor, 'linewidth', 2)
        %          plot(ax(axRTSSD), (goSSDTimepoints - offsetVal), inhRTSSD{i}, '-', 'color', inhColor, 'linewidth', 2);
        %       elseif iPct < 50
        %          %             inhColor = [i/nSignalStrength, i/nSignalStrength, i/nSignalStrength];
        %          inhColor = [(4-i)/4, (4-i)/4, (4-i)/4];
        %          plot(ax(axInhEach), timePoints, inhEach{i}, 'color', inhColor, 'linewidth', 2)
        %          plot(ax(axRTSSD), (goSSDTimepoints - offsetVal), inhRTSSD{i}, '-', 'color', inhColor, 'linewidth', 2);
        %       elseif iPct > 50
        %          %             inhColor = [1 - i/nSignalStrength, 1 - i/nSignalStrength, 1 - i/nSignalStrength];
        %          inhColor = [1 - (nSignalStrength + 1 - i)/4, 1 - (nSignalStrength + 1 - i)/4, 1 - (nSignalStrength + 1 - i)/4];
        %          plot(ax(axInhEach), timePoints, inhEach{i}, '--', 'color', inhColor, 'linewidth', 2)
        %          plot(ax(axRTSSD), (goSSDTimepoints - offsetVal), inhRTSSD{i}, '--', 'color', inhColor, 'linewidth', 2);
        %       end
        
        
    end
    
    %     xlim(ax(axInhEach), [min(ssdArrayGrand) max(ssdArrayGrand)])
    ylim(ax(axInhEach), [0 1])
    %     xlim(ax(axRTSSD), [min(-dataRTSSD) max(-dataRTSSD)])
    switch subjectID
        case 'Human'
            xlim(ax(axRTSSD), [-600 300])
        otherwise
            xlim(ax(axRTSSD), [-300 200])
    end
    ylim(ax(axRTSSD), [0 1])
    
    
    
    
    % INDIVIDUAL SESSION INHIBITION FUNCITONS
    for i = 1 : nSession
        plot(ax(axInhSess), inhSessTime{i}, inhSess{i}, 'k', 'linewidth', 2)
    end
    ylim([0 1])
    set(ax(axInhSess), 'Xlim', [0 max(ssdArrayGrand)])
    
    
    
    
    % OBSERVED VS. PREDICTED RTS  ********************************
    d = stopRespondRT - stopRespondRTPredict;
    %     d(:,1:2) = [];
    nanmean(d(:))
    nanstd(d(:))
    if plotSurface
        figure(64)
        clf
        obsSurf = surf(ssdArray(find(~nanCol)), pSignalArray, stopRespondRT);
        set(obsSurf,'FaceColor',[1 0 0],'FaceAlpha',0.7);
        hold on
        predSurf = surf(ssdArray(find(~nanCol)), pSignalArray, stopRespondRTPredict);
        set(predSurf,'FaceColor',[.4 .4 .4],'FaceAlpha',0.7);
    end
    stopRespondRT = round(stopRespondRT)
    stopRespondRTPredict = round(stopRespondRTPredict)
    
    
    if options.printPlot
        print(figureHandle,fullfile(local_figure_path, subjectID,['ccm_population_inhibition_',options.saveName]),'-dpdf', '-r300')
    end
end

% Print out plotted values
ssrtPrint = [ssrtPlotAvg; ssrtPlotStd; ssrtPlotSem];


data.ssrtIntWeight = ssrtIntWeight;
data.ssrtGrandIntWeight = ssrtGrandIntWeight;
data.rtGo = rtGo;
data.nStopMean = nStopMean;



% Print out plotted values
% fprintf('\nSSRT:   Targ: %.0f (%.0f) (%.0f)', goRightProbMean, goRightProbStd, goRightProbSem);




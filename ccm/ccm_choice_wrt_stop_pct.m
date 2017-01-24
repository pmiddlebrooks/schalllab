function data = ccm_choice_wrt_stop_pct(subjectID, sessionID, trialLag, timeConst, plotFlag, figureHandle)

% A series of analyses to test how choice varies as a function of stopping
%%

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
targ1PropArray = ExtraVar.pSignalArray;
ssdArray = ExtraVar.ssdArray;

% Flag to determine whether we want to include stop trial analyses for the
% session
DO_STOPS = 1;

if ~strcmp(SessionData.taskID, 'ccm')
   fprintf('Not a choice countermanding session, try again\n')
   return
end




if nargin < 3
   trialLag = 20;
end
if nargin < 4
   timeConst = 5;
end
if nargin < 5
   plotFlag = 1;
end
if nargin < 6
   figureHandle = 4444;
end

if plotFlag
   nRow = 3;
   nColumn = 2;
   screenOrSave = 'save';
   [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = big_figure(nRow, nColumn, figureHandle, screenOrSave);
   clf
   choicePlotXMargin = .03;
   ssdMargin = 20;
   ylimArray = [];
   
   % Colors
   stopSparseColor = [0 0 0];
   stopDenseColor = [1 0 0];
   stopDenseMatchColor = [.8 .3 .3];
   goSparseColor = [0 0 0];
   goDenseColor = [0 0 1];
   goDenseMatchColor = [.3 .3 .8];
   
   % axes names
   axGo = 1;
   axStop = 2;
   axRun = 3;
   
   ax(axGo) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
   cla
   hold(ax(axGo), 'on')
   
   ax(axStop) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
   cla
   hold(ax(axStop), 'on')
   
   ax(axRun) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 1) yAxesPosition(2, 1) axisWidth*nColumn axisHeight]);
   cla
   hold(ax(axRun), 'on')
end



ssdArrayRaw = trialData.stopSignalOn - trialData.responseCueOn;





% ******************************************************************
% PSYCHOMETRIC FN -VS- PCT STOPPING
% ******************************************************************


% Find indices of usable trials, and get rid of the rest:
doNotIncludeTrial = strcmp(trialData.trialOutcome, 'fixationAbort') | strcmp(trialData.trialOutcome, 'noFixation') | ...
   strcmp(trialData.trialOutcome, 'earlySaccade') | strcmp(trialData.trialOutcome, 'choiceStimulusAbort'); %| ...

includeTrial = ~doNotIncludeTrial;

rt = trialData.responseOnset - trialData.responseCueOn;
rt = rt(includeTrial);





% BREAKING DOWN BY RECENT PCT STOPPING
% Of the included (valid) trials, assign zeros to go trials and ones to
% stop trials
validStop = ssdArrayRaw(includeTrial);
validStop(~isnan(validStop)) = 1;
validStop(isnan(validStop)) = 0;

% trialLag = 40;
% timeConst = 20; % trials
stopPctRunMean = nan(length(validStop), 1);
rtRunMean = nan(length(validStop), 1);
trialsPast = [-trialLag : -1];
convFn = 1 - exp(timeConst ./ trialsPast)';

% trialLag2 = 10;
% stopPctRunMean2 = nan(1, length(validTrial));
% trialsPast2 = [-trialLag2 : -1];
% timeConst2 = 5; % trials
% convFn2 = 1 - exp(timeConst2 ./ trialsPast2)';

for iTrial = trialLag : length(validStop)
   stopPctRunMean(iTrial) = sum(validStop(iTrial - trialLag + 1 : iTrial) .* convFn) / sum(convFn);
   iNanRTTrial = isnan(rt(iTrial - trialLag + 1 : iTrial));
   rtRunMean(iTrial) = nansum(rt(iTrial - trialLag + 1 : iTrial) .* convFn) / sum(convFn(~iNanRTTrial));
   %     stopPctRunMean2(iTrial) = sum(validTrial(iTrial - trialLag2 + 1 : iTrial) .* convFn2) / sum(convFn2);
end


stopPctMean = mean(validStop);
rtMean = nanmean(rt);





% sparseStop = stopPctRunMean < stopPctMean;
% denseStop = stopPctRunMean > stopPctMean;

sparseStopTrial = stopPctRunMean < stopPctMean;
denseStopTrial = stopPctRunMean > stopPctMean;



%
% % Expand the usable trials back into the whole trialData size, indexing
% % appropriate trials for each condition
% denseStopTrial = zeros(size(trialData.trialOutcome, 1), 1);
% denseStopTrial(includeTrial) = denseStop;
% sparseStopTrial = zeros(size(trialData.trialOutcome, 1), 1);
% sparseStopTrial(includeTrial) = sparseStop;

% sparseStopTrial = find(sparseStop);
% denseStopTrial = find(denseStop);



sparseRT = rt(sparseStopTrial);
denseTrialRT = sortrows([rt(denseStopTrial), find(denseStopTrial)]);
% size(lateSSDTrialRT)
% nanmean(lateSSDTrialRT(:,1))
% nanmean(sparseSSDRT)
while nanmean(denseTrialRT(:,1)) > nanmean(sparseRT)
   denseTrialRT(end,:) = [];
end
% size(lateSSDTrialRT)
denseMatchStopTrial = denseTrialRT(:,2);





% Get incorrect stop probabilities rightward
goRightProbSparse = zeros(length(targ1PropArray), 1);
nGoSparse = nan(1, length(targ1PropArray));
nGoSparseRight = nan(1, length(targ1PropArray));
goSparseTrial = [];
goSparseRightTrial = [];

goRightProbDense = zeros(length(targ1PropArray), 1);
nGoDense = nan(1, length(targ1PropArray));
nGoDenseRight = nan(1, length(targ1PropArray));
goDenseTrial = [];
goDenseRightTrial = [];

goRightProbDenseMatch = zeros(length(targ1PropArray), 1);
nGoDenseMatch = nan(1, length(targ1PropArray));
nGoDenseMatchRight = nan(1, length(targ1PropArray));
goDenseMatchTrial = [];
goDenseMatchRightTrial = [];

stopRightProbSparse = zeros(length(targ1PropArray), 1);
nStopSparse = nan(1, length(targ1PropArray));
nStopSparseRight = nan(1, length(targ1PropArray));
stopSparseTrial = [];
stopSparseRightTrial = [];

stopRightProbDense = zeros(length(targ1PropArray), 1);
nStopDense = nan(1, length(targ1PropArray));
nStopDenseRight = nan(1, length(targ1PropArray));
stopDenseTrial = [];
stopDenseRightTrial = [];

stopRightProbDenseMatch = zeros(length(targ1PropArray), 1);
nStopDenseMatch = nan(1, length(targ1PropArray));
nStopDenseMatchRight = nan(1, length(targ1PropArray));
stopDenseMatchTrial = [];
stopDenseMatchRightTrial = [];



% Get default trial selection options
selectOpt = ccm_trial_selection;

for iPropIndex = 1 : length(targ1PropArray);
   iPct = targ1PropArray(iPropIndex) * 100;
   
   selectOpt.rightCheckerPct = iPct;
   
   % All go correct trials
   selectOpt.ssd       = 'none';
   selectOpt.targDir       = 'collapse';
   selectOpt.outcome   = {'goCorrectTarget', 'targetHoldAbort','goCorrectDistractor', 'distractorHoldAbort'};
   goTrial = ccm_trial_selection(trialData, selectOpt);
   
   iGoSparseTrial = intersect(find(sparseStopTrial), goTrial);
   goSparseTrial = [goSparseTrial; iGoSparseTrial];
   iGoDenseTrial = intersect(find(denseStopTrial), goTrial);
   goDenseTrial = [goDenseTrial; iGoDenseTrial];
   iGoDenseMatchTrial = intersect(find(denseMatchStopTrial), goTrial);
   goDenseMatchTrial = [goDenseMatchTrial; iGoDenseMatchTrial];
   nGoSparse(iPropIndex) = length(iGoSparseTrial);
   nGoDense(iPropIndex) = length(iGoDenseTrial);
   nGoDenseMatch(iPropIndex) = length(iGoDenseMatchTrial);
   
   % All go correct rightward trials
   selectOpt.targDir    = 'right';
   selectOpt.outcome    = {'goCorrectTarget', 'targetHoldAbort'};
   goTargetRight        = ccm_trial_selection(trialData, selectOpt);
   selectOpt.targDir  	= 'left';
   selectOpt.outcome    = {'goCorrectDistractor', 'distractorHoldAbort'};
   goDistractorRight    = ccm_trial_selection(trialData, selectOpt);
   goRightTrial         = union(goTargetRight, goDistractorRight);
   
   iGoSparseRightTrial = intersect(find(sparseStopTrial), goRightTrial);
   goSparseRightTrial = [goSparseRightTrial; iGoSparseRightTrial];
   iGoDenseRightTrial = intersect(find(denseStopTrial), goRightTrial);
   goDenseRightTrial = [goDenseRightTrial; iGoDenseRightTrial];
   iGoDenseMatchRightTrial = intersect(find(denseMatchStopTrial), goRightTrial);
   goDenseMatchRightTrial = [goDenseMatchRightTrial; iGoDenseMatchRightTrial];
   nGoSparseRight(iPropIndex) = length(iGoSparseRightTrial);
   nGoDenseRight(iPropIndex) = length(iGoDenseRightTrial);
   nGoDenseMatchRight(iPropIndex) = length(iGoDenseMatchRightTrial);
   
   
   goRightProbSparse(iPropIndex) = nansum(nGoSparseRight(:, iPropIndex)) / nansum(nGoSparse(:, iPropIndex)) ;
   goRightProbDense(iPropIndex) = nansum(nGoDenseRight(:, iPropIndex)) / nansum(nGoDense(:, iPropIndex)) ;
   goRightProbDenseMatch(iPropIndex) = nansum(nGoDenseMatchRight(:, iPropIndex)) / nansum(nGoDenseMatch(:, iPropIndex)) ;
   
   
   
   
   
   % All stop incorrect trials
   selectOpt.ssd       = 'collapse';
   selectOpt.outcome   = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget', ...
      'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
   selectOpt.targDir       = 'collapse';
   stopTrial = ccm_trial_selection(trialData, selectOpt);
   
   iStopSparseTrial = intersect(find(sparseStopTrial), stopTrial);
   stopSparseTrial = [stopSparseTrial; iStopSparseTrial];
   iStopDenseTrial = intersect(find(denseStopTrial), stopTrial);
   stopDenseTrial = [stopDenseTrial; iStopDenseTrial];
   iStopDenseMatchTrial = intersect(find(denseMatchStopTrial), stopTrial);
   stopDenseMatchTrial = [stopDenseMatchTrial; iStopDenseMatchTrial];
   nStopSparse(iPropIndex) = length(iStopSparseTrial);
   nStopDense(iPropIndex) = length(iStopDenseTrial);
   nStopDenseMatch(iPropIndex) = length(iStopDenseMatchTrial);
   
   
   % All stop incorrect rightward trials
   selectOpt.outcome   = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
   selectOpt.targDir       = 'right';
   stopTargetRight = ccm_trial_selection(trialData, selectOpt);
   selectOpt.outcome   = {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
   selectOpt.targDir       = 'left';
   stopDistractorRight = ccm_trial_selection(trialData, selectOpt);
   stopRightTrial = union(stopTargetRight, stopDistractorRight);
   
   iStopSparseRightTrial = intersect(find(sparseStopTrial), stopRightTrial);
   stopSparseRightTrial = [stopSparseRightTrial; iStopSparseRightTrial];
   iStopDenseRightTrial = intersect(find(denseStopTrial), stopRightTrial);
   stopDenseRightTrial = [stopDenseRightTrial; iStopDenseRightTrial];
   iStopDenseMatchRightTrial = intersect(find(denseMatchStopTrial), stopRightTrial);
   stopDenseMatchRightTrial = [stopDenseMatchRightTrial; iStopDenseMatchRightTrial];
   nStopSparseRight(iPropIndex) = length(iStopSparseRightTrial);
   nStopDenseRight(iPropIndex) = length(iStopDenseRightTrial);
   nStopDenseMatchRight(iPropIndex) = length(iStopDenseMatchRightTrial);
   
   
   
   stopRightProbSparse(iPropIndex) = nansum(nStopSparseRight(:, iPropIndex)) / nansum(nStopSparse(:, iPropIndex));
   stopRightProbDense(iPropIndex) = nansum(nStopDenseRight(:, iPropIndex)) / nansum(nStopDense(:, iPropIndex));
   stopRightProbDenseMatch(iPropIndex) = nansum(nStopDenseMatchRight(:, iPropIndex)) / nansum(nStopDenseMatch(:, iPropIndex));
   
end




propPoints = targ1PropArray(1) : .001 : targ1PropArray(end);

% Weibull fit the percent stop signal data:
% Go trials:
[sparseParam, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(goRightProbSparse)), goRightProbSparse(~isnan(goRightProbSparse)));
goRightSlopeSparse = sparseParam(2);
goSparsePsychFn = weibull_curve(sparseParam, propPoints);

[denseParam, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(goRightProbDense)), goRightProbDense(~isnan(goRightProbDense)));
goRightSlopeDense = denseParam(2);
goDensePsychFn = weibull_curve(denseParam, propPoints);

[denseMatchParam, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(goRightProbDenseMatch)), goRightProbDenseMatch(~isnan(goRightProbDenseMatch)));
goRightSlopeDenseMatch = denseMatchParam(2);
goDenseMatchPsychFn = weibull_curve(denseMatchParam, propPoints);



% Stop trials:
[sparseParam, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(stopRightProbSparse)), stopRightProbSparse(~isnan(stopRightProbSparse)));
stopRightSlopeSparse = sparseParam(2);
stopSparsePsychFn = weibull_curve(sparseParam, propPoints);

[denseParam, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(stopRightProbDense)), stopRightProbDense(~isnan(stopRightProbDense)));
stopRightSlopeDense = denseParam(2);
stopDensePsychFn = weibull_curve(denseParam, propPoints);

[denseMatchParam, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(stopRightProbDenseMatch)), stopRightProbDenseMatch(~isnan(stopRightProbDenseMatch)));
stopRightSlopeDenseMatch = denseMatchParam(2);
stopDenseMatchPsychFn = weibull_curve(denseMatchParam, propPoints);





if plotFlag
   % Sparse vs Dense
   stopSparseLine = plot(ax(axStop), propPoints, stopSparsePsychFn, '-', 'color', stopSparseColor, 'linewidth', 2);
   plot(ax(axStop), targ1PropArray, stopRightProbSparse, 'o', 'color', stopSparseColor, 'linewidth', 2, 'markerfacecolor', stopSparseColor, 'markeredgecolor', stopSparseColor)
   stopDenseLine = plot(ax(axStop), propPoints, stopDensePsychFn, '-', 'color', stopDenseColor, 'linewidth', 2);
   plot(ax(axStop), targ1PropArray, stopRightProbDense, 'o', 'color', stopDenseColor, 'linewidth', 2, 'markerfacecolor', stopDenseColor, 'markeredgecolor', stopDenseColor)
   stopDenseMatchLine = plot(ax(axStop), propPoints, stopDenseMatchPsychFn, '--', 'color', stopDenseMatchColor, 'linewidth', 2);
   plot(ax(axStop), targ1PropArray, stopRightProbDenseMatch, 'o', 'color', stopDenseMatchColor, 'linewidth', 2, 'markerfacecolor', [1 1 1], 'markeredgecolor', stopDenseMatchColor)
   plot(ax(axStop), [.5 .5], ylim, '--k')
   legend([stopSparseLine, stopDenseLine, stopDenseMatchLine], 'Stop Sparse', 'Stop Dense', 'Stop Dense Match', 'location', 'northwest')
   
   set(ax(axStop), 'xtick', targ1PropArray)
   set(ax(axStop), 'xtickLabel', targ1PropArray)
   set(get(ax(axStop), 'ylabel'), 'String', 'p(Right)')
   set(ax(axStop),'XLim',[targ1PropArray(1) - choicePlotXMargin targ1PropArray(end) + choicePlotXMargin])
   set(ax(axStop),'YLim',[0 1])
   
   
   goSparseLine = plot(ax(axGo), propPoints, goSparsePsychFn, '-', 'color', goSparseColor, 'linewidth', 2);
   plot(ax(axGo), targ1PropArray, goRightProbSparse, 'o', 'color', goSparseColor, 'linewidth', 2, 'markerfacecolor', goSparseColor, 'markeredgecolor', goSparseColor)
   goDenseLine = plot(ax(axGo), propPoints, goDensePsychFn, '-', 'color', goDenseColor, 'linewidth', 2);
   plot(ax(axGo), targ1PropArray, goRightProbDense, 'o', 'color', goDenseColor, 'linewidth', 2, 'markerfacecolor', goDenseColor, 'markeredgecolor', goDenseColor)
   goDenseMatchLine = plot(ax(axGo), propPoints, goDenseMatchPsychFn, '--', 'color', goDenseMatchColor, 'linewidth', 2);
   plot(ax(axGo), targ1PropArray, goRightProbDenseMatch, 'o', 'color', goDenseMatchColor, 'linewidth', 2, 'markerfacecolor', [1 1 1], 'markeredgecolor', goDenseMatchColor)
   plot(ax(axGo), [.5 .5], ylim, '--k')
   legend([goSparseLine, goDenseLine, goDenseMatchLine], 'Go Sparse', 'Go Dense', 'Go Dense Match', 'location', 'northwest')
   
   set(ax(axGo), 'xtick', targ1PropArray)
   set(ax(axGo), 'xtickLabel', targ1PropArray)
   set(get(ax(axGo), 'ylabel'), 'String', 'p(Right)')
   set(ax(axGo),'XLim',[targ1PropArray(1) - choicePlotXMargin targ1PropArray(end) + choicePlotXMargin])
   set(ax(axGo),'YLim',[0 1])
   
   
   
   % Running avg data
   [AX,H1,H2] = plotyy(ax(axRun), 1:length(validStop), stopPctRunMean, 1:length(rt), rtRunMean);
   set(AX(1), 'ylim', [min(stopPctRunMean)-.02 max(stopPctRunMean)]+.02)
   set(AX(2), 'ylim', [min(rtRunMean)-15 max(rtRunMean)]+15)
   set(AX(1), 'xlim', [0 length(validStop)])
   set(AX(2), 'xlim', [0 length(validStop)])
   set(get(AX(1),'Ylabel'),'String','Stop Trial %')
   set(get(AX(2),'Ylabel'),'String','RT (ms)')
   xlabel('Trial #')
   set(H1,'color','r', 'linewidth', 2)
   set(H2,'color','k', 'linewidth', 2)
   ttl = sprintf('Running Avg Stop %% and RT. Tau = %d ms, Lag = %d trials', timeConst, trialLag);
   title(ttl)
   legend('Stop %', 'RT')
   
   
end

% Sparse/Dense data
data.stopPctMean = stopPctMean;
data.stopPctRunMean = stopPctRunMean;

data.goSparseTrial = goSparseTrial;
data.goDenseTrial = goDenseTrial;
data.goDenseMatchTrial = goDenseMatchTrial;
data.goSparseRightTrial = goSparseRightTrial;
data.goDenseRightTrial = goDenseRightTrial;
data.goDenseMatchRightTrial = goDenseMatchRightTrial;
data.stopSparseTrial = stopSparseTrial;
data.stopDenseTrial = stopDenseTrial;
data.stopDenseMatchTrial = stopDenseMatchTrial;
data.stopSparseRightTrial = stopSparseRightTrial;
data.stopDenseRightTrial = stopDenseRightTrial;
data.stopDenseMatchRightTrial = stopDenseMatchRightTrial;

data.goRightProbSparse = goRightProbSparse;
data.goRightProbDense = goRightProbDense;
data.goRightProbDenseMatch = goRightProbDenseMatch;
data.stopRightProbSparse = stopRightProbSparse;
data.stopRightProbDense = stopRightProbDense;
data.stopRightProbDenseMatch = stopRightProbDenseMatch;

data.nGoSparse = nGoSparse;
data.nGoSparseRight = nGoSparseRight;
data.nStopSparse = nStopSparse;
data.nStopSparseRight = nStopSparseRight;
data.nGoDense = nGoDense;
data.nGoDenseMatch = nGoDenseMatch;
data.nGoDenseRight = nGoDenseRight;
data.nGoDenseMatchRight = nGoDenseMatchRight;
data.nStopDense = nStopDense;
data.nStopDenseMatch = nStopDenseMatch;
data.nStopDenseRight = nStopDenseRight;
data.nStopDenseMatchRight = nStopDenseMatchRight;

data.goRightSlopeSparse = goRightSlopeSparse;
data.goRightSlopeDense = goRightSlopeDense;
data.goRightSlopeDenseMatch = goRightSlopeDenseMatch;
data.stopRightSlopeSparse = stopRightSlopeSparse;
data.stopRightSlopeDense = stopRightSlopeDense;
data.stopRightSlopeDenseMatch = stopRightSlopeDenseMatch;



if plotFlag
   print(figure(figureHandle), ['~/schalllab/tempfigures/choiceStopPct_',subjectID, '_', sessionID, '_Tau', num2str(timeConst), '_Lag', num2str(trialLag)], '-dpdf');
end

return

% ******************************************************************
% RT -VS- SSD
% ******************************************************************

[signalStrengthLeft, signalStrengthRight ...
   goLeftToTarg, goRightToTarg, goLeftToDist, goRightToDist ...
   stopLeftToTarg, stopRightToTarg, stopLeftToDist, stopRightToDist] ...
   = ccm_chronometric(subjectID, sessionID, 0);

stopTargRT = [stopLeftToTarg, stopRightToTarg];
stopDistRT = [stopRightToDist, stopLeftToTarg];

if plotFlag
   % SSD vs RT
   ax(RTvSSD) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
   cla
   hold(ax(RTvSSD), 'on')
end


minColorGun = .25;
maxColorGun = 1;
for iPropIndex = 1 : length(targ1PropArray);
   iPercent = targ1PropArray(iPropIndex) * 100;
   
   % Determine color to use for plot based on which checkerboard color
   % proportion being used. Normalize the available color spectrum to do
   % it
   if iPercent == 50
      signalColor = [0 0 0];
   elseif iPercent < 50
      colorNorm = .5 - targ1PropArray(1);
      colorProp = (.5 - targ1PropArray(iPropIndex)) / colorNorm;
      colorGun = minColorGun + (maxColorGun - minColorGun) * colorProp;
      signalColor = [0 colorGun colorGun];
   elseif iPercent > 50
      colorNorm = targ1PropArray(end) - .5;
      colorProp = (targ1PropArray(iPropIndex) - .5) / colorNorm;
      colorGun = minColorGun + (maxColorGun - minColorGun) * colorProp;
      signalColor = [colorGun 0 colorGun];
   end
   
   for jSSDIndex = 1 : length(ssdArray)
      jSSD = ssdArray(jSSDIndex);
      
      rtMean(jSSDIndex, iPropIndex) = mean(stopTargRT{jSSDIndex, iPropIndex});
   end  % SSD
   
   
   iSSDIndices = ~isnan(rtMean(:, iPropIndex));
   rtMean(:, iPropIndex)
   
   
   [p, s] = polyfit(ssdArray(iSSDIndices), rtMean(iSSDIndices, iPropIndex), 1);
   [y, delta] = polyval(p, ssdArray(iSSDIndices), s);
   stats = regstats(ssdArray(iSSDIndices), rtMean(iSSDIndices, iPropIndex))
   % fprintf('p-value for regression: %.4f\n', stats.tstat.pval(2))
   R = corrcoef(ssdArray(iSSDIndices), rtMean(iSSDIndices, iPropIndex));
   Rsqrd = R(1, 2)^2;
   cov(ssdArray(iSSDIndices), rtMean(iSSDIndices, iPropIndex));
   xVal = min(ssdArray(iSSDIndices)) : .001 : max(ssdArray(iSSDIndices));
   yVal = p(1) * xVal + p(2);
   
   if plotFlag
      %             plot(ax(RTvSSD), ssdArray(iSSDIndices), rtMean(iSSDIndices, iPropIndex), '-', 'color', signalColor)
      plot(ax(RTvSSD), xVal, yVal, '-', 'color', signalColor, 'lineWidth', 2)
   end
end  % chekcer Proportion






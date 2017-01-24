function data = ccm_choice_wrt_ssd(subjectID, sessionID, trialLag, timeConst, plotFlag, figureHandle)

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
   
   %Colors
   stopEarlyColor = [0 0 0];
   stopLateColor = [1 0 0];
   stopLateMatchColor = [.8 .3 .3];
   goEarlyColor = [0 0 0];
   goLateColor = [0 0 1];
   goLateMatchColor = [.3 .3 .8];
   
   % axes names
   axRun = 2;
   axStopSSD = 3;
   axGoSSD = 1;
   % SSD vs PSYCH
   ax(axGoSSD) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
   cla
   hold(ax(axGoSSD), 'on')
   
   ax(axStopSSD) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
   cla
   hold(ax(axStopSSD), 'on')
   
   ax(axRun) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 1) yAxesPosition(2, 1) axisWidth*nColumn axisHeight]);
   cla
   hold(ax(axRun), 'on')
end




ssdArrayRaw = trialData.stopSignalOn - trialData.responseCueOn;










% ******************************************************************
% PSYCHOMETRIC FN -VS- SSD
% ******************************************************************

% Find indices of usable trials, and get rid of the rest:
doNotIncludeTrial = strcmp(trialData.trialOutcome, 'fixationAbort') | strcmp(trialData.trialOutcome, 'noFixation') | ...
   strcmp(trialData.trialOutcome, 'earlySaccade') | strcmp(trialData.trialOutcome, 'choiceStimulusAbort'); %| ...

includeTrial = ~doNotIncludeTrial;

rt = trialData.responseOnset - trialData.responseCueOn;
rt = rt(includeTrial);



% Of the included (valid) trials, assign zeros to go trials and ones to
% stop trials
ssdArrayRaw = ssdArrayRaw(includeTrial);




% BREAKING DOWN BY SSD

% trialLag = 40;
% timeConst = 20; % trials
ssdRunMean  = nan(length(ssdArrayRaw), 1);
rtRunMean   = nan(length(ssdArrayRaw), 1);
trialsPast  = [-trialLag : -1];
convFn      = 1 - exp(timeConst ./ trialsPast)';

for iTrial = trialLag : length(ssdArrayRaw)
   iNanSSDTrial = isnan(ssdArrayRaw(iTrial - trialLag + 1 : iTrial));
   ssdRunMean(iTrial) = nansum(ssdArrayRaw(iTrial - trialLag + 1 : iTrial) .* convFn) / sum(convFn(~iNanSSDTrial));
   iNanRTTrial = isnan(rt(iTrial - trialLag + 1 : iTrial));
   rtRunMean(iTrial) = nansum(rt(iTrial - trialLag + 1 : iTrial) .* convFn) / sum(convFn(~iNanRTTrial));
end

ssdMean = nanmean(ssdArrayRaw);
rtMean = nanmean(rt);


% earlySSDTrial = ssdRunMean < ssdMean & includeTrial;
% lateSSDTrial = ssdRunMean > ssdMean & includeTrial;
earlySSDTrial = ssdRunMean < ssdMean;
lateSSDTrial = ssdRunMean > ssdMean;

earlySSDRT = rt(earlySSDTrial);
lateSSDTrialRT = sortrows([rt(lateSSDTrial), find(lateSSDTrial)]);
% size(lateSSDTrialRT)
% nanmean(lateSSDTrialRT(:,1))
% nanmean(earlySSDRT)
while nanmean(lateSSDTrialRT(:,1)) > nanmean(earlySSDRT)
   lateSSDTrialRT(end,:) = [];
end
% size(lateSSDTrialRT)
lateSSDMatchRTTrial = lateSSDTrialRT(:,2);
% nanmean(lateSSDRT)






goRightProbEarly = zeros(length(targ1PropArray), 1);
nGoEarly = nan(1, length(targ1PropArray));
nGoEarlyRight = nan(1, length(targ1PropArray));
goEarlySSDTrial = [];
goEarlySSDRightTrial = [];

goRightProbLate = zeros(length(targ1PropArray), 1);
nGoLate = nan(1, length(targ1PropArray));
nGoLateRight = nan(1, length(targ1PropArray));
goLateSSDTrial = [];
goLateSSDRightTrial = [];

goRightProbLateMatchRT = zeros(length(targ1PropArray), 1);
nGoLateMatchRT = nan(1, length(targ1PropArray));
nGoLateRightMatchRT = nan(1, length(targ1PropArray));
goLateSSDMatchRTTrial = [];
goLateSSDRightMatchRTTrial = [];


stopRightProbEarly = zeros(length(targ1PropArray), 1);
nStopEarly = nan(1, length(targ1PropArray));
nStopEarlyRight = nan(1, length(targ1PropArray));
stopEarlySSDTrial = [];
stopEarlySSDRightTrial = [];

stopRightProbLate = zeros(length(targ1PropArray), 1);
nStopLate = nan(1, length(targ1PropArray));
nStopLateRight = nan(1, length(targ1PropArray));
stopLateSSDTrial = [];
stopLateSSDRightTrial = [];

stopRightProbLateMatchRT = zeros(length(targ1PropArray), 1);
nStopLateMatchRT = nan(1, length(targ1PropArray));
nStopLateRightMatchRT = nan(1, length(targ1PropArray));
stopLateSSDMatchRTTrial = [];
stopLateSSDRightMatchRTTrial = [];


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
   
   iGoEarlySSDTrial           	= intersect(find(earlySSDTrial), goTrial);
   goEarlySSDTrial             	= [goEarlySSDTrial; iGoEarlySSDTrial];
   iGoLateSSDTrial             	= intersect(find(lateSSDTrial), goTrial);
   goLateSSDTrial            	= [goLateSSDTrial; iGoLateSSDTrial];
   iGoLateSSDMatchRTTrial      	= intersect(lateSSDMatchRTTrial, goTrial);
   goLateSSDMatchRTTrial      	= [goLateSSDMatchRTTrial; iGoLateSSDMatchRTTrial];
   nGoEarly(iPropIndex)       	= length(iGoEarlySSDTrial);
   nGoLate(iPropIndex)        	= length(iGoLateSSDTrial);
   nGoLateMatchRT(iPropIndex)   = length(iGoLateSSDMatchRTTrial);
   
   % All go correct rightward trials
   selectOpt.targDir    = 'right';
   selectOpt.outcome    = {'goCorrectTarget', 'targetHoldAbort'};
   goTargetRight        = ccm_trial_selection(trialData, selectOpt);
   selectOpt.targDir  	= 'left';
   selectOpt.outcome    = {'goCorrectDistractor', 'distractorHoldAbort'};
   goDistractorRight    = ccm_trial_selection(trialData, selectOpt);
   goRightTrial         = union(goTargetRight, goDistractorRight);
   
   iGoEarlySSDRightTrial               = intersect(find(earlySSDTrial), goRightTrial);
   goEarlySSDRightTrial                = [goEarlySSDRightTrial; iGoEarlySSDRightTrial];
   iGoLateSSDRightTrial                = intersect(find(lateSSDTrial), goRightTrial);
   goLateSSDRightTrial                 = [goLateSSDRightTrial; iGoLateSSDRightTrial];
   iGoLateSSDRightMatchRTTrial         = intersect(lateSSDMatchRTTrial, goRightTrial);
   goLateSSDRightMatchRTTrial          = [goLateSSDRightMatchRTTrial; iGoLateSSDRightMatchRTTrial];
   nGoEarlyRight(iPropIndex)           = length(iGoEarlySSDRightTrial);
   nGoLateRight(iPropIndex)            = length(iGoLateSSDRightTrial);
   nGoLateRightMatchRT(iPropIndex)     = length(iGoLateSSDRightMatchRTTrial);
   
   goRightProbEarly(iPropIndex)        = nansum(nGoEarlyRight(:, iPropIndex)) / nansum(nGoEarly(:, iPropIndex)) ;
   goRightProbLate(iPropIndex)         = nansum(nGoLateRight(:, iPropIndex)) / nansum(nGoLate(:, iPropIndex)) ;
   goRightProbLateMatchRT(iPropIndex)  = nansum(nGoLateRightMatchRT(:, iPropIndex)) / nansum(nGoLateMatchRT(:, iPropIndex)) ;
   
   
   
   
   
   
   
   % All stop incorrect trials
   selectOpt.ssd       = 'collapse';
   selectOpt.outcome   = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget', ...
      'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
   selectOpt.targDir       = 'collapse';
   stopTrial = ccm_trial_selection(trialData, selectOpt);
   
   
   iStopEarlySSDTrial = intersect(find(earlySSDTrial), stopTrial);
   stopEarlySSDTrial = [stopEarlySSDTrial; iStopEarlySSDTrial];
   iStopLateSSDTrial = intersect(find(lateSSDTrial), stopTrial);
   stopLateSSDTrial = [stopLateSSDTrial; iStopLateSSDTrial];
   iStopLateSSDMatchRTTrial = intersect(lateSSDMatchRTTrial, stopTrial);
   stopLateSSDMatchRTTrial = [stopLateSSDMatchRTTrial; iStopLateSSDMatchRTTrial];
   nStopEarly(iPropIndex) = length(iStopEarlySSDTrial);
   nStopLate(iPropIndex) = length(iStopLateSSDTrial);
   nStopLateMatchRT(iPropIndex) = length(iStopLateSSDMatchRTTrial);
   
   % All stop incorrect rightward trials
   selectOpt.outcome   = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
   selectOpt.targDir       = 'right';
   stopTargetRight = ccm_trial_selection(trialData, selectOpt);
   selectOpt.outcome   = {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
   selectOpt.targDir       = 'left';
   stopDistractorRight = ccm_trial_selection(trialData, selectOpt);
   stopRightTrial = union(stopTargetRight, stopDistractorRight);
   
   
   iStopEarlySSDRightTrial = intersect(find(earlySSDTrial), stopRightTrial);
   stopEarlySSDRightTrial = [stopEarlySSDRightTrial; iStopEarlySSDRightTrial];
   iStopLateSSDRightTrial = intersect(find(lateSSDTrial), stopRightTrial);
   stopLateSSDRightTrial = [stopLateSSDRightTrial; iStopLateSSDRightTrial];
   iStopLateSSDRightMatchRTTrial = intersect(lateSSDMatchRTTrial, stopRightTrial);
   stopLateSSDRightMatchRTTrial = [stopLateSSDRightMatchRTTrial; iStopLateSSDRightMatchRTTrial];
   nStopEarlyRight(iPropIndex) = length(iStopEarlySSDRightTrial);
   nStopLateRight(iPropIndex) = length(iStopLateSSDRightTrial);
   nStopLateRightMatchRT(iPropIndex) = length(iStopLateSSDRightMatchRTTrial);
   
   
   stopRightProbEarly(iPropIndex) = nansum(nStopEarlyRight(:, iPropIndex)) / nansum(nStopEarly(:, iPropIndex)) ;
   stopRightProbLate(iPropIndex) = nansum(nStopLateRight(:, iPropIndex)) / nansum(nStopLate(:, iPropIndex)) ;
   stopRightProbLateMatchRT(iPropIndex) = nansum(nStopLateRightMatchRT(:, iPropIndex)) / nansum(nStopLateMatchRT(:, iPropIndex)) ;
end




propPoints = targ1PropArray(1) : .001 : targ1PropArray(end);




% Weibull fit the early/late SSD data

% Go trials
% [earlyParam, lowestSSE] = Weibull(targ1PropArray, goRightProbEarly);
[earlyParam, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(goRightProbEarly)), goRightProbEarly(~isnan(goRightProbEarly)));
goRightSlopeEarly = earlyParam(2);
goEarlyPsychFn = weibull_curve(earlyParam, propPoints);
% [lateParam, lowestSSE] = Weibull(targ1PropArray, goRightProbLate);
[lateParam, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(goRightProbLate)), goRightProbLate(~isnan(goRightProbLate)));
goRightSlopeLate = lateParam(2);
goLatePsychFn = weibull_curve(lateParam, propPoints);
% [lateParamMatchRT, lowestSSE] = Weibull(targ1PropArray, goRightProbLateMatchRT);
[lateParamMatchRT, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(goRightProbLateMatchRT)), goRightProbLateMatchRT(~isnan(goRightProbLateMatchRT)));
goRightSlopeLateMatch = lateParamMatchRT(2);
goLatePsychFnMatchRT = weibull_curve(lateParamMatchRT, propPoints);


% Stop trials
% [earlyParam, lowestSSE] = Weibull(targ1PropArray, stopRightProbEarly);
[earlyParam, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(stopRightProbEarly)), stopRightProbEarly(~isnan(stopRightProbEarly)));
stopRightSlopeEarly = earlyParam(2);
stopEarlyPsychFn = weibull_curve(earlyParam, propPoints);
% [lateParam, lowestSSE] = Weibull(targ1PropArray, stopRightProbLate);
[lateParam, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(stopRightProbLate)), stopRightProbLate(~isnan(stopRightProbLate)));
stopRightSlopeLate = lateParam(2);
stopLatePsychFn = weibull_curve(lateParam, propPoints);
% [lateParamMatchRT, lowestSSE] = Weibull(targ1PropArray, stopRightProbLateMatchRT);
[lateParamMatchRT, lowestSSE] = psychometric_weibull_fit(targ1PropArray(~isnan(stopRightProbLateMatchRT)), stopRightProbLateMatchRT(~isnan(stopRightProbLateMatchRT)));
stopRightSlopeLateMatch = lateParamMatchRT(2);
stopLatePsychFnMatchRT = weibull_curve(lateParamMatchRT, propPoints);


if plotFlag
   
   
   stopEarlyLine = plot(ax(axStopSSD), propPoints, stopEarlyPsychFn, '-', 'color', stopEarlyColor, 'linewidth', 2);
   plot(ax(axStopSSD), targ1PropArray, stopRightProbEarly, 'o', 'color', stopEarlyColor, 'linewidth', 2, 'markerfacecolor', stopEarlyColor, 'markeredgecolor', stopEarlyColor)
   stopLateLine = plot(ax(axStopSSD), propPoints, stopLatePsychFn, '--', 'color', stopLateColor, 'linewidth', 2);
   plot(ax(axStopSSD), targ1PropArray, stopRightProbLate, 'o', 'color', stopLateColor, 'linewidth', 2, 'markerfacecolor', stopLateColor, 'markeredgecolor', stopLateColor)
   stopLateLineMatchRT = plot(ax(axStopSSD), propPoints, stopLatePsychFnMatchRT, '-', 'color', stopLateColor, 'linewidth', 2);
   plot(ax(axStopSSD), targ1PropArray, stopRightProbLateMatchRT, 'o', 'color', stopLateMatchColor, 'linewidth', 2, 'markerfacecolor', [1 1 1], 'markeredgecolor', stopLateMatchColor)
   
   set(ax(axStopSSD), 'xtick', targ1PropArray)
   set(ax(axStopSSD), 'xtickLabel', targ1PropArray)
   set(get(ax(axStopSSD), 'ylabel'), 'String', 'p(Right)')
   set(ax(axStopSSD),'XLim',[targ1PropArray(1) - choicePlotXMargin targ1PropArray(end) + choicePlotXMargin])
   set(ax(axStopSSD),'YLim',[0 1])
   legend([stopEarlyLine, stopLateLine, stopLateLineMatchRT], 'Stop Early', 'Stop Late', 'Stop Late LM', 'location', 'northwest')
   plot(ax(axStopSSD), [.5 .5], ylim, '--k')
   title(ax(axStopSSD), 'Stop Trials Early vs Late')
   
   
   goEarlyLine = plot(ax(axGoSSD), propPoints, goEarlyPsychFn, '-', 'color', goEarlyColor, 'linewidth', 2);
   plot(ax(axGoSSD), targ1PropArray, goRightProbEarly, 'o', 'color', goEarlyColor, 'linewidth', 2, 'markerfacecolor', goEarlyColor, 'markeredgecolor', goEarlyColor)
   goLateLine = plot(ax(axGoSSD), propPoints, goLatePsychFn, '--', 'color', goLateColor, 'linewidth', 2);
   plot(ax(axGoSSD), targ1PropArray, goRightProbLate, 'o', 'color', goLateColor, 'linewidth', 2, 'markerfacecolor', goLateColor, 'markeredgecolor', goLateColor)
   goLateLineMatchRT = plot(ax(axGoSSD), propPoints, goLatePsychFnMatchRT, '-', 'color', goLateColor, 'linewidth', 2);
   plot(ax(axGoSSD), targ1PropArray, goRightProbLateMatchRT, 'o', 'color', goLateMatchColor, 'linewidth', 2, 'markerfacecolor', [1 1 1], 'markeredgecolor', goLateMatchColor)
   
   set(ax(axGoSSD), 'xtick', targ1PropArray)
   set(ax(axGoSSD), 'xtickLabel', targ1PropArray)
   set(get(ax(axGoSSD), 'ylabel'), 'String', 'p(Right)')
   set(ax(axGoSSD),'XLim',[targ1PropArray(1) - choicePlotXMargin targ1PropArray(end) + choicePlotXMargin])
   set(ax(axGoSSD),'YLim',[0 1])
   legend([goEarlyLine, goLateLine, goLateLineMatchRT], 'Go Early', 'Go Late', 'Go Late LM', 'location', 'northwest')
   plot(ax(axGoSSD), [.5 .5], ylim, '--k')
   title(ax(axGoSSD), 'Go Trials Early vs Late')
   
   
   % Running avg data
   [AX,H1,H2] = plotyy(ax(axRun), 1:length(ssdArrayRaw), ssdRunMean, 1:length(rt), rtRunMean);
   set(AX(1), 'ylim', [min(ssdRunMean)-.02 max(ssdRunMean)]+.02)
   set(AX(2), 'ylim', [min(rtRunMean)-15 max(rtRunMean)]+15)
   set(AX(1), 'xlim', [0 length(ssdArrayRaw)])
   set(AX(2), 'xlim', [0 length(ssdArrayRaw)])
   set(get(AX(1),'Ylabel'),'String','SSD Running Avg')
   set(get(AX(2),'Ylabel'),'String','RT (ms)')
   xlabel('Trial #')
   set(H1,'color','r', 'linewidth', 2)
   set(H2,'color','k', 'linewidth', 2)
   ttl = sprintf('Running Avg SSD and RT. Tau = %d ms, Lag = %d trials', timeConst, trialLag);
   title(ttl)
   legend('SSD', 'RT')
end




% Early/Late data
% Go Trials
data.goEarlySSDTrial = goEarlySSDTrial;
data.goLateSSDTrial = goLateSSDTrial;
data.goLateSSDMatchRTTrial = goLateSSDMatchRTTrial;
data.goEarlySSDRightTrial = goEarlySSDRightTrial;
data.goLateSSDRightTrial = goLateSSDRightTrial;
data.goLateSSDRightMatchRTTrial = goLateSSDRightMatchRTTrial;

data.goRightProbEarly = goRightProbEarly;
data.goRightProbLate = goRightProbLate;
data.goRightProbLateMatchRT = goRightProbLateMatchRT;

data.nGoEarly = nGoEarly;
data.nGoEarlyRight = nGoEarlyRight;
data.nGoLate = nGoLate;
data.nGoLateRight = nGoLateRight;
data.nGoLateMatchRT = nGoLateMatchRT;
data.nGoLateRightMatchRT = nGoLateRightMatchRT;

data.goRightSlopeEarly = goRightSlopeEarly;
data.goRightSlopeLate = goRightSlopeLate;
data.goRightSlopeLateMatch = goRightSlopeLateMatch;

% Stop trials
data.stopEarlySSDTrial = stopEarlySSDTrial;
data.stopLateSSDTrial = stopLateSSDTrial;
data.stopLateSSDMatchRTTrial = stopLateSSDMatchRTTrial;
data.stopEarlySSDRightTrial = stopEarlySSDRightTrial;
data.stopLateSSDRightTrial = stopLateSSDRightTrial;
data.stopLateSSDRightMatchRTTrial = stopLateSSDRightMatchRTTrial;

data.stopRightProbEarly = stopRightProbEarly;
data.stopRightProbLate = stopRightProbLate;
data.stopRightProbLateMatchRT = stopRightProbLateMatchRT;

data.nStopEarly = nStopEarly;
data.nStopEarlyRight = nStopEarlyRight;
data.nStopLate = nStopLate;
data.nStopLateRight = nStopLateRight;
data.nStopLateMatchRT = nStopLateMatchRT;
data.nStopLateRightMatchRT = nStopLateRightMatchRT;

data.stopRightSlopeEarly = stopRightSlopeEarly;
data.stopRightSlopeLate = stopRightSlopeLate;
data.stopRightSlopeLateMatch = stopRightSlopeLateMatch;





if plotFlag
   print(figure(figureHandle), ['~/schalllab/tempfigures/choiceSSSD_',subjectID, '_', sessionID, '_Tau', num2str(timeConst), '_Lag', num2str(trialLag)], '-dpdf');
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






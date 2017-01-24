function data = ccm_n2pc_population(subjectID, varargin)

% Set defaults
task = 'ccm';
plotFlag    = 1;
printPlot   = 0;
figureHandle = 6124;
alignEvent = 'checkerOn';
sessionSet = 'behavior1';
stopSession = true;
filterData = false;
for i = 1 : 2 : length(varargin)
   switch varargin{i}
      case 'plotFlag'
         plotFlag = varargin{i+1};
      case 'printPlot'
         printPlot = varargin{i+1};
      case 'figureHandle'
         figureHandle = varargin{i+1};
      case 'alignEvent'
         alignEvent = varargin{i+1};
      case 'stopSession'
         stopSession = varargin{i+1};
      case 'filterData'
         filterData = varargin{i+1};
      otherwise
   end
end

subjectID = lower(subjectID);

%%%%%%%%%%%      Constants    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
stopHz      = 50;




switch subjectID
   case 'human'
      pSignalArray = [.35 .42 .46 .5 .54 .58 .65];
   case 'broca'
      pSignalArray = [.41 .45 .48 .5 .52 .55 .59];
   case 'xena'
      pSignalArray = [.35 .42 .47 .5 .53 .58 .65];
end
[sessionArray, subjectIDArray] = task_session_array(subjectID, task, sessionSet, stopSession);
% pSignalArray = unique([populationData.signalStrengthLeft{1}, populationData.signalStrengthRight{1}])';




nSession    = length(sessionArray);
nRow        = 2;
nColumn     = 2;
axLGo       = 1;
axRGo       = 2;
axLStop     = 3;
axRStop     = 4;
if printPlot
   [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
else
   [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
end
clf
% [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nRow, nColumn, 'portrait', figureHandle);
ax(axLGo) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
ax(axRGo) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
ax(axLStop) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 1) yAxesPosition(2, 1) axisWidth axisHeight]);
ax(axRStop) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2) yAxesPosition(2, 2) axisWidth axisHeight]);


%Colors
cMap = ccm_colormap(pSignalArray);
leftColor = cMap(1,:) .* .8;
rightColor = cMap(end,:) .* .8;
% rightColor = [1 0 0];
% leftColor = [0 0 1];
stopContraColor = [1 0 0];
stopIpsiColor = [0 0 1];
diffColor = [.4 .4 .4];

hold(ax(axLGo), 'on')
hold(ax(axRGo), 'on')
hold(ax(axLStop), 'on')
hold(ax(axRStop), 'on')




o1RightTargEEG = cell(nSession, 1);
o1LeftTargEEG = cell(nSession, 1);
o2LeftTargEEG = cell(nSession, 1);
o2RightTargEEG = cell(nSession, 1);
rtR             = nan(nSession, 1);
rtL             = nan(nSession, 1);
emptyArray = [];
for iSession = 1 : nSession
   iSessionID = sessionArray{iSession};
   disp(iSessionID)
   
   iData = ccm_n2pc(subjectID, iSessionID, 'alignEvent', alignEvent, 'plotFlag', 0);
   
   if ~isempty(iData)
     
      if filterData
         o1RightTargEEG{iSession} = lowpass(nanmean(iData.o1RightTargEEG, 1), stopHz);
         o1LeftTargEEG{iSession} = lowpass(nanmean(iData.o1LeftTargEEG, 1), stopHz);
         
         o2LeftTargEEG{iSession} = lowpass(nanmean(iData.o2LeftTargEEG, 1), stopHz);
         o2RightTargEEG{iSession} = lowpass(nanmean(iData.o2RightTargEEG, 1), stopHz);
         
      else
         
         o1RightTargEEG{iSession} = nanmean(iData.o1RightTargEEG, 1);
         o1LeftTargEEG{iSession} = nanmean(iData.o1LeftTargEEG, 1);
         
         o2LeftTargEEG{iSession} = nanmean(iData.o2LeftTargEEG, 1);
         o2RightTargEEG{iSession} = nanmean(iData.o2RightTargEEG, 1);
      end
      
      rtR(iSession) = nanmean(iData.rtR);
      rtL(iSession) = nanmean(iData.rtL);
   else
      emptyArray = [emptyArray, iSession];
   end
end
displayRange = iData.displayRange;
alignTime = iData.alignTime;
% Remove sessions that did not collect data;
o1RightTargEEG(emptyArray)  = [];
o1LeftTargEEG(emptyArray)   = [];
o2LeftTargEEG(emptyArray)   = [];
o2RightTargEEG(emptyArray)  = [];
rtR(emptyArray)             = [];
rtL(emptyArray)             = [];

lContraMean = nanmean(cell2mat(cellfun(@(x) nanmean(x, 1), o1RightTargEEG, 'uniformoutput', false)), 1);
lIpsiMean   = nanmean(cell2mat(cellfun(@(x) nanmean(x, 1), o1LeftTargEEG, 'uniformoutput', false)), 1);
olDiff      = lContraMean - lIpsiMean;

rContraMean = nanmean(cell2mat(cellfun(@(x) nanmean(x, 1), o2LeftTargEEG, 'uniformoutput', false)), 1);
rIpsiMean   = nanmean(cell2mat(cellfun(@(x) nanmean(x, 1), o2RightTargEEG, 'uniformoutput', false)), 1);
o2Diff      = rContraMean - rIpsiMean;


rtMeanR = mean(rtR);
rtMeanL = mean(rtL);
% lContraMean = mean(cell2mat(o1RightTargEEG), 1);
% lIpsiMean = mean(cell2mat(o1LeftTargEEG), 1);
% olDiff = lContraMean - lIpsiMean;
%
% rContraMean = mean(cell2mat(o2LeftTargEEG), 1);
% rIpsiMean = mean(cell2mat(o2RightTargEEG), 1);
% orDiff = rContraMean - rIpsiMean;

% xMax = 1.1 * max(rtMeanL, rtMeanR);
xMax = 500;

if plotFlag
   plot(ax(axLGo), olDiff, 'color', diffColor, 'linewidth', 2)
   plot(ax(axLGo), lContraMean, 'color', rightColor, 'linewidth', 2)
   plot(ax(axLGo), lIpsiMean, 'color', leftColor, 'linewidth', 2)
   ylim(ax(axLGo), [-.015 .015])
   plot(ax(axLGo), [alignTime alignTime], ylim(ax(axLGo)))
   xlim(ax(axLGo), [0 xMax])
   if strcmp(alignEvent, 'checkerOn')
      plot(ax(axLGo), [rtMeanR + -displayRange(1)  rtMeanR + -displayRange(1)], ylim(ax(axLGo)), '--', 'color', rightColor)
      plot(ax(axLGo), [rtMeanL + -displayRange(1)  rtMeanL + -displayRange(1)], ylim(ax(axLGo)), '--', 'color', leftColor)
   end
   %       set(ax(axLGo), 'XTickLabel', (-alignTime : 100 : xMax))
   set(ax(axLGo), 'xtick', (0 : 100 : displayRange(2) - displayRange(1)))
   set(ax(axLGo), 'XTickLabel', [displayRange(1) : 100 : displayRange(2)])
   legend(ax(axLGo),  'Contra - Ipsi', 'o1 Right (Contra)', 'o1 Left (Ipsi)')
   set(gca,'YDir','reverse');
   
   plot(ax(axRGo), o2Diff, 'color', diffColor, 'linewidth', 2)
   plot(ax(axRGo), rContraMean, 'color', leftColor, 'linewidth', 2)
   plot(ax(axRGo), rIpsiMean, 'color', rightColor, 'linewidth', 2)
   ylim(ax(axRGo), [-.015 .015])
   plot(ax(axRGo), [alignTime alignTime], ylim(ax(axRGo)))
   xlim(ax(axRGo), [0 xMax])
   if strcmp(alignEvent, 'checkerOn')
      plot(ax(axRGo), [rtMeanR + -displayRange(1)  rtMeanR + -displayRange(1)], ylim(ax(axRGo)), '--', 'color', rightColor)
      plot(ax(axRGo), [rtMeanL + -displayRange(1)  rtMeanL + -displayRange(1)], ylim(ax(axRGo)), '--', 'color', leftColor)
   end
   %       set(ax(axRGo), 'XTickLabel', (-alignTime : 100 : xMax))
   set(ax(axRGo), 'xtick', (0 : 100 : displayRange(2) - displayRange(1)))
   set(ax(axRGo), 'XTickLabel', [displayRange(1) : 100 : displayRange(2)])
   legend(ax(axRGo),  'Contra - Ipsi', 'o2 Left (Contra)', 'o2 Right (Ipsi)')
end


data.o1RightTargEEG = o1RightTargEEG;
data.o1LeftTargEEG = o1LeftTargEEG;
data.o2LeftTargEEG = o2LeftTargEEG;
data.o2RightTargEEG = o2RightTargEEG;




if printPlot
   localFigurePath = local_figure_path;
   if ~stopSession
      stopText = '_no_stop_sessions';
   else
      stopText = [];
   end
   print(figureHandle,[localFigurePath, subjectID, '_ccm_n2pc_pop_', alignEvent, stopText],'-dpdf', '-r300')
end


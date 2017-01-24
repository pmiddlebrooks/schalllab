function checkerboard = ccm_checkerboard_session(subjectID, sessionID, plotFlag)

% The checkerboard consists of zeros and ones. Zeros are right target
% color, ones are left target color
%%
if nargin < 3, plotFlag = 1; end
normalizeData = 0;
blackWhite = 0;

% subjectID = 'pg';
% sessionID = 'Allkeypress';
% subjectID = 'hu';
% sessionID = 'Allsaccade';

% subjectID = 'Broca';
% sessionID = 'bp093n02';
% sessionID = 'bp085n02';



% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
rightSignalPArray = ExtraVar.pSignalArray;

if iscell(trialData.checkerArray)
   % reshape cell contents if necessary
   if sum(cellfun(@(x) size(x, 1)>1, trialData.checkerArray))
      trialData.checkerArray = cellfun(@(x) x', trialData.checkerArray, 'uniformoutput', false);
   end
   trialData.checkerArray = cell2mat(trialData.checkerArray);
end
nCheckerRow = sqrt(length(trialData.checkerArray(1,:)));
nCheckerCol = nCheckerRow;



DO_50 = false;
if DO_50
   % Add a second 50% signal strength to distinguish between targ1 and targ2
   % trials
   if ismember(.5, rightSignalPArray)
      [a,i] = ismember(.5, rightSignalPArray);
      rightSignalPArray = [rightSignalPArray(1 : i) ; .5; rightSignalPArray(i+1 : end)];
   end
else
   if ismember(.5, rightSignalPArray)
      [a,i] = ismember(.5, rightSignalPArray);
      rightSignalPArray(i) = [];
   end
end
nSignal = length(rightSignalPArray);



if plotFlag
   nColumn = nSignal;
   nRow = nColumn;
   squareAxes = true;
   checkerFig = 487;
   [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, checkerFig, squareAxes);
   clf
   
   rowGoTarg = 1;
   rowGoDist = 2;
   rowStopTarg = 3;
   rowStopDist = 4;
   rowStopStop = 5;
   
   if blackWhite
      %                 invGray = 1 - gray;
      %                 colormap(invGray);
      colormap(gray);
   else
      colormap(flipud(cool(64)));
   end
end






selectOpt = ccm_trial_selection;


for i = 1 : length(rightSignalPArray);
   iPct = rightSignalPArray(i) * 100;
   checkerboard(i).rightCheckerP = rightSignalPArray(i);
   
   
   % If it's not 50% or if there's only one 50% condition in
   % targPropArray
   if iPct ~= 50 || (iPct == 50 &&  rightSignalPArray(i-1) ~= 50 &&  rightSignalPArray(i+1) ~= 50)
      targetHemifield = 'all';
      % If it's the first (left target) 50% condition
   elseif iPct == 50 && rightSignalPArray(i-1) ~= 50
      targetHemifield = 'left';
      % If it's the second (right target) 50% condition
   elseif iPct == 50 && rightSignalPArray(i-1) == 50
      targetHemifield = 'right';
   end
   
   %   Go trials:
   selectOpt.ssd       = 'none';
   selectOpt.rightCheckerPct = iPct;
   
   selectOpt.outcome   = {'goCorrectTarget'; 'targetHoldAbort'};
   iGoTargTrial        = ccm_trial_selection(trialData, selectOpt);
   selectOpt.outcome   = {'goCorrectDistractor', 'distractorHoldAbort'};
   iGoDistTrial        = ccm_trial_selection(trialData, selectOpt);
   
   
   %   Stop trials: collapse across ssds
   selectOpt.ssd       = 'collapse';
   
   selectOpt.outcome   = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
   iStopTargTrial      = ccm_trial_selection(trialData, selectOpt);
   selectOpt.outcome   = {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
   iStopDistTrial      = ccm_trial_selection(trialData, selectOpt);
   selectOpt.outcome   = {'stopCorrect'};
   iStopStopTrial      = ccm_trial_selection(trialData, selectOpt);
   
   
   checkerboard(i).goTarg = trialData.checkerArray(iGoTargTrial,:);
   checkerboard(i).goDist = trialData.checkerArray(iGoDistTrial,:);
   checkerboard(i).stopTarg = trialData.checkerArray(iStopTargTrial,:);
   checkerboard(i).stopDist = trialData.checkerArray(iStopDistTrial,:);
   checkerboard(i).stopStop = trialData.checkerArray(iStopStopTrial,:);
   
   if plotFlag
      %         ttl = sprintf('Target:  %s', epochName);
      %         title(ttl)
      iGoTarg = mean(checkerboard(i).goTarg, 1) .* size(colormap, 1);
      iGoDist = mean(checkerboard(i).goDist, 1) .* size(colormap, 1);
      iStopTarg = mean(checkerboard(i).stopTarg, 1) .* size(colormap, 1);
      iStopDist = mean(checkerboard(i).stopDist, 1) .* size(colormap, 1);
      iStopStop = mean(checkerboard(i).stopStop, 1) .* size(colormap, 1);
      
      if normalizeData
         allData = [iGoTarg, iGoDist, iStopTarg, iStopDist, iStopStop];
         iGoTarg = (iGoTarg - min(allData)) ./ (max(allData) - min(allData));
         iGoDist = (iGoDist - min(allData)) ./ (max(allData) - min(allData));
         iStopTarg = (iStopTarg - min(allData)) ./ (max(allData) - min(allData));
         iStopDist = (iStopDist - min(allData)) ./ (max(allData) - min(allData));
         iStopStop = (iStopStop - min(allData)) ./ (max(allData) - min(allData));
      end
      
      % Go Targ
      ax(rowGoTarg, i) = axes('units', 'centimeters', 'position', [xAxesPosition(rowGoTarg, i) yAxesPosition(rowGoTarg, i) axisWidth axisHeight]);
      %                 hold(ax(rowGoTarg, i), 'on')
      cla
      iGoTarg = reshape(iGoTarg, nCheckerRow, nCheckerCol)';
      image(iGoTarg)
      set(gca,'XTickLabel','','XTick',[],'YTickLabel','','YTick',[])
      
      % Go Dist
      ax(rowGoDist, i) = axes('units', 'centimeters', 'position', [xAxesPosition(rowGoDist, i) yAxesPosition(rowGoDist, i) axisWidth axisHeight]);
      %         hold(ax(rowGoDist, i), 'on')
      cla
      iGoDist = reshape(iGoDist, nCheckerRow, nCheckerCol)';
      image(iGoDist)
      set(gca,'XTickLabel','','XTick',[],'YTickLabel','','YTick',[])
      
      % Stop Targ
      ax(rowStopTarg, i) = axes('units', 'centimeters', 'position', [xAxesPosition(rowStopTarg, i) yAxesPosition(rowStopTarg, i) axisWidth axisHeight]);
      %         hold(ax(rowStopTarg, i), 'on')
      cla
      iStopTarg = reshape(iStopTarg, nCheckerRow, nCheckerCol)';
      image(iStopTarg)
      set(gca,'XTickLabel','','XTick',[],'YTickLabel','','YTick',[])
      
      % Stop Dist
      ax(rowStopDist, i) = axes('units', 'centimeters', 'position', [xAxesPosition(rowStopDist, i) yAxesPosition(rowStopDist, i) axisWidth axisHeight]);
      %         hold(ax(rowStopDist, i), 'on')
      cla
      iStopDist = reshape(iStopDist, nCheckerRow, nCheckerCol)';
      image(iStopDist)
      set(gca,'XTickLabel','','XTick',[],'YTickLabel','','YTick',[])
      
      % Stop Stop
      ax(rowStopStop, i) = axes('units', 'centimeters', 'position', [xAxesPosition(rowStopStop, i) yAxesPosition(rowStopStop, i) axisWidth axisHeight]);
      %         hold(ax(rowStopStop, i), 'on')
      cla
      iStopStop = reshape(iStopStop, nCheckerRow, nCheckerCol)';
      image(iStopStop)
      set(gca,'XTickLabel','','XTick',[],'YTickLabel','','YTick',[])
      
      
      % Labeling stuff
      targL = sprintf('%d%%', round(iPct));
      xlabel(ax(rowStopStop, i), targL,  'fontsize', 18, 'fontweight', 'bold')
      
      
      if i == 1
         % Go Targ
         goTargL = sprintf('Go Target');
         ylabel(ax(rowGoTarg, i), goTargL,  'fontsize', 16, 'fontweight', 'bold')
         
         % Go Dist
         goDistL = sprintf('Go Distractor');
         ylabel(ax(rowGoDist, i), goDistL, 'fontsize', 16, 'fontweight', 'bold')
         
         % Stop Targ
         stopTargL = sprintf('Stop Target');
         ylabel(ax(rowStopTarg, i), stopTargL, 'fontsize', 16, 'fontweight', 'bold')
         
         % Stop Dist
         stopDistL = sprintf('Stop Distractor');
         ylabel(ax(rowStopDist, i),stopDistL, 'fontsize', 16, 'fontweight', 'bold')
         
         % Go Targ
         stopStopL = sprintf('Stop Stop');
         ylabel(ax(rowStopStop, i), stopStopL, 'fontsize', 16, 'fontweight', 'bold')
         
      end
      
   end
end

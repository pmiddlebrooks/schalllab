function checkerboard = ccm_checkerboard_population(subjectID, plotFlag)

% The checkerboard consists of zeros and ones. Zeros are right target
% color, ones are left target color
%% Population SSRT
%******************************************************************************
fprintf('\n\n\n\n')
disp('*******************************************************************************')
disp('Populaiton Checkerboard')
if nargin < 2
   plotFlag = 1;
end


task = 'ccm';
% subjectID = 'Human';
% subjectID = 'Xena';
% subjectID = 'Broca';

% switch subjectID
%     case 'Human'
%         signalStrength = [.35 .42 .46 .5 .54 .58 .65];
%     case 'Broca'
%         signalStrength = [.41 .45 .48 .5 .52 .55 .59];
%     case 'Xena'
%         signalStrength = [.35 .42 .47 .5 .53 .58 .65];
% end
% nSignalStrength = length(signalStrength);


normalizeData = 0;
blackWhite = 0;


[sessionArray, subjectIDArray] = task_session_array(subjectID, task);




nSession = length(sessionArray);



% Initialize variables to be filled *****************
signalPArray = [];


% First loop through the data files to determine the signal strength values
% used throughout (it may be consistent across sessions, or may not)
for iSession = 1 : nSession
   
   
   % Load the data
   [trialData, SessionData, ExtraVar] = load_data(subjectID, sessionArray{iSession});
   pSignalArray = ExtraVar.pSignalArray;
   ssdArray = ExtraVar.ssdArray;
   
   DO_50 = true;
   if DO_50
      % Add a second 50% signal strength to distinguish between targ1 and targ2
      % trials
      if ismember(.5, pSignalArray)
         [a,i] = ismember(.5, pSignalArray);
         pSignalArray = [pSignalArray(1 : i) ; .5; pSignalArray(i+1 : end)];
      end
   else
      if ismember(.5, pSignalArray)
         [a,i] = ismember(.5, pSignalArray);
         pSignalArray(i) = [];
      end
   end
   
   
   %     session(iSession) = sessionArray{iSession};
   signalPArray = [signalPArray; pSignalArray(:)];
   
   
end
signalPArray = unique(signalPArray);
nSignal = length(signalPArray);






% Now loop through sessions and build a struct that can readily be analyzed
checkerArray = [];
checkerboard(i).goTargArray = [];
checkerboard(i).goDistArray = [];
checkerboard(i).stopTargArray = [];
checkerboard(i).stopDistArray = [];
checkerboard(i).stopStopArray = [];
for iSession = 1 : nSession
   
   % Load the data
   [trialData, SessionData, ExtraVar] = load_data(subjectID, sessionArray{iSession});
   
   [a, b] = ismember('checkerArray', fieldnames(trialData));
   if ~a
      fprintf('%s not doesn''t have a checkerArray varialbe, skipping it\n', sessionArray{iSession})
      continue
   else
      
      
      sessionArray{iSession}
      if iscell(trialData.checkerArray)
         % reshape cell contents if necessary
         if sum(cellfun(@(x) size(x, 1)>1, trialData.checkerArray))
            trialData.checkerArray = cellfun(@(x) x', trialData.checkerArray, 'uniformoutput', false);
         end
         trialData.checkerArray = cell2mat(trialData.checkerArray);
      end
      nCheckerRow = sqrt(length(trialData.checkerArray(1,:)));
      nCheckerCol = nCheckerRow;
      
      
      checkerArray = [checkerArray; trialData.checkerArray];
      
      
      
      selectOpt = ccm_trial_selection;
      
      
      for i = 1 : length(signalPArray);
         iPct = signalPArray(i) * 100;
         
         % Initialize each trial outcome array to NaN, b/c sessions may or
         % many not have used a given signal strength
         iGoTargArray = nan(1, nCheckerRow * nCheckerCol);
         iGoDistArray = iGoTargArray;
         iStopTargArray = iGoTargArray;
         iStopDistArray = iGoTargArray;
         iStopStopArray = iGoTargArray;
         
         checkerboard(i).rightCheckerP = signalPArray(i);
         
         
         % If it's not 50% or if there's only one 50% condition in
         % targPropArray
         if iPct ~= 50 || (iPct == 50 &&  signalPArray(i-1) ~= 50 &&  signalPArray(i+1) ~= 50)
            selectOpt.targDir = 'collapse';
            % If it's the first (left target) 50% condition
         elseif iPct == 50 && signalPArray(i-1) ~= 50
            selectOpt.targDir = 'left';
            % If it's the second (right target) 50% condition
         elseif iPct == 50 && signalPArray(i-1) == 50
            selectOpt.targDir = 'right';
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
         
         % If there were trials with the given signal strength during this
         % session, use them instead of the default NaN vector
         if ~isempty(iGoTargTrial)
            iGoTargArray = mean(trialData.checkerArray(iGoTargTrial,:), 1);
            iGoDistArray = mean(trialData.checkerArray(iGoDistTrial,:), 1);
            iStopTargArray = mean(trialData.checkerArray(iStopTargTrial,:), 1);
            iStopDistArray = mean(trialData.checkerArray(iStopDistTrial,:), 1);
            iStopStopArray = mean(trialData.checkerArray(iStopStopTrial,:), 1);
         end
         
         checkerboard(i).goTargArray = [checkerboard(i).goTargArray; iGoTargArray];
         checkerboard(i).goDistArray = [checkerboard(i).goDistArray; iGoDistArray];
         checkerboard(i).stopTargArray = [checkerboard(i).stopTargArray; iStopTargArray];
         checkerboard(i).stopDistArray = [checkerboard(i).stopDistArray; iStopDistArray];
         checkerboard(i).stopStopArray = [checkerboard(i).stopStopArray; iStopStopArray];
      end
      
      
      
      
   end % if checkerArray is part of trialData variables
end









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
   
   
   
   
   
   
   for i = 1 : length(signalPArray);
      
      iPct = signalPArray(i) * 100;
      %         ttl = sprintf('Target:  %s', epochName);
      %         title(ttl)
      iGoTarg = nanmean(checkerboard(i).goTargArray, 1) .* size(colormap, 1);
      iGoDist = nanmean(checkerboard(i).goDistArray, 1) .* size(colormap, 1);
      iStopTarg = nanmean(checkerboard(i).stopTargArray, 1) .* size(colormap, 1);
      iStopDist = nanmean(checkerboard(i).stopDistArray, 1) .* size(colormap, 1);
      iStopStop = nanmean(checkerboard(i).stopStopArray, 1) .* size(colormap, 1);
      
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

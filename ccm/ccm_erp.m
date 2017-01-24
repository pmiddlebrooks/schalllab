function Data = ccm_erp(subjectID, sessionID, epochArray, varargin)

%
% function Data = ccm_erp(subjectID, sessionID, varargin)
%
% Single neuron analyses for choice countermanding task. Only plots the
% sdfs. To see rasters, use ccm_single_neuron_rasters, which displays all
% conditions in a given epoch
%
% input:
%   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%
%   varargin: property names and their values:
%           'plotFlag': 0 or 1
%           'printPlot': 0 or 1: If set to 1, this prints the figure in the local_figures folder
%           'unitArray': a single unit, like 'spikeUnit17a', or an array of units, like {'spikeUnit17a', 'spikeUnit17b'}
%
%
% Returns Unit structure with fields:
%
%   Unit.signalStrength(x).(condition).ssd(x).(epoch name)
%
%   condition can be:  goTarg, goDist, stopTarg, stopDist, stopStop
%   ssd(x):  only applies for stop trials, else the field is absent
%   epoch name: fixOn, targOn, checkerOn, etc.
%   nGo
%   nGoRight
%   nStopIncorrect
%   nStopIncorrectRight
%   goRightLogical
%   goRightSignalStrength
%   stopRightLogical
%   stopRightSignalStrength

%%
% subjectID = 'Broca';
% sessionID = 'bp127n02';
% % sessionID = 'bp088n02';
% % sessionID = 'bp089n02';
% sessionID = 'bp090n02';
% sessionID = 'bp091n02';
% sessionID = 'bp092n02';
% sessionID = 'bp093n02';
% % sessionID = 'bp106n01';
% % sessionID = 'bp095n02';

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID, [ccm_min_vars, 'eegData']);
pSignalArray    = ExtraVar.pSignalArray;
ssdArray        = ExtraVar.ssdArray;

[a, b] = ismember('eegData', fieldnames(trialData));
if ~a;
   fprintf('Session %s does not contain eeg data \n', sessionID)
   Data = [];
   return
end


% Set defaults
plotFlag        = 1;
printPlot       = 0; % If set to 1, this collapses data across signal strengths (within each hemifield)
figureHandle    = 320;
electrodeArray  = eeg_electrode_map(subjectID);
filterData      = false;
for i = 1 : 2 : length(varargin)
    switch varargin{i}
        case 'plotFlag'
            plotFlag = varargin{i+1};
        case 'electrodeArray'
            electrodeArray = varargin{i+1};
        case 'printPlot'
            printPlot = varargin{i+1};
       case 'figureHandle'
            figureHandle = varargin{i+1};
      case 'filterData'
         filterData = varargin{i+1};
        otherwise
    end
end
% electrodeArray = {'o2'};
nElectrode  = length(electrodeArray);
% epochArray = {alignEvent};
      nEpoch = length(epochArray);

clear Data
% constants
cropWindow  = -499 : 500;  % used to extract a semi-small portion of signal for each epoch/alignemnt
baseWindow 	= -149 : 0;   % To baseline-shift the eeg signals, relative to event alignment index;
stopHz      = 30;

DO_STOPS    = 1;
MIN_RT      = 120;
MAX_RT      = 1200;
STD_MULTIPLE = 3;
nTrial      = size(trialData.trialOutcome, 1);

% Get rid of trials with outlying RTs
[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTD);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);


% If we're probing the end of the saccade/response, need to add a column of
% that data to trialData dataset:
if sum(ismember(epochArray, 'responseOffset'))
   trialData.responseOffset = nan(nTrial, 1);
   ind = cellfun(@(x,y) find(x == y),num2cell(trialData.responseOnset), trialData.saccBegin,'uni',false);
for i = 1 : nTrial
   if ~isempty(ind{i})
         trialData.responseOffset(i) = trialData.responseOnset(i) + trialData.saccDuration{i}(ind{i});
   end
end


DO_50 = false;
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
nSignal = length(pSignalArray);




% If there weren't stop trials, skip all stop-related analyses
if isempty(ssdArray) || ~DO_STOPS
   Data = [];
   disp('ccm_inhibition.m: No stop trials or stop trial analyses not requested');
   return
end





[minTrialDuration,i] = min(trialData.trialDuration(~strcmp(trialData.trialOutcome, 'fixationAbort')));
[maxTrialDuration,i] = max(trialData.trialDuration(~strcmp(trialData.trialOutcome, 'fixationAbort')));
nCondition  = 5;
eegMax      = zeros(nElectrode, nEpoch, nSignal, nCondition);  % Keep track of maximum eeg values, for setting y-axis limits in plots
eegMin      = zeros(nElectrode, nEpoch, nSignal, nCondition);  % Keep track of maximum eeg values, for setting y-axis limits in plots

for kChanIndex = 1 : nElectrode

   chanName         = electrodeArray{kChanIndex};
   [a, kChannel]    = ismember(chanName, eeg_electrode_map(subjectID));
   
   % If there is not an electrode for that monkey, inform user and skip
   if ~a || length(trialData.eegData(1,:)) < kChanIndex
      fprintf('Electrode %s is not part of set for %s, skipping it \n\n', chanName, subjectID)
      continue
   end
   
   
   Data(kChanIndex).subjectID   = subjectID;
   Data(kChanIndex).sessionID   = sessionID;
   Data(kChanIndex).name        = electrodeArray{kChanIndex};
   Data(kChanIndex).ssdArray    = ssdArray;
   Data(kChanIndex).pSignalArray = pSignalArray;
   
   for iPropIndex = 1 : length(pSignalArray);
      iPct = pSignalArray(iPropIndex) * 100;
      
      
      %   Go trials:
      ssdRange = 'none';
      
      % If it's not 50% or if there's only one 50% condition in
      % targPropArray
      if iPct ~= 50 || (iPct == 50 &&  pSignalArray(iPropIndex-1) ~= 50 &&  pSignalArray(iPropIndex+1) ~= 50)
         targetHemifield = 'all';
         % If it's the first (left target) 50% condition
      elseif iPct == 50 && pSignalArray(iPropIndex-1) ~= 50
         targetHemifield = 'left';
         % If it's the second (right target) 50% condition
      elseif iPct == 50 && pSignalArray(iPropIndex-1) == 50
         targetHemifield = 'right';
      end
      iGoTargTrial = ccm_trial_selection(trialData, {'goCorrectTarget'; 'targetHoldAbort'}, iPct, ssdRange, targetHemifield);
      iGoDistTrial = ccm_trial_selection(trialData, {'goCorrectDistractor', 'distractorHoldAbort'}, iPct, ssdRange, targetHemifield);
      
      for mEpoch = 1 : length(epochArray)
            alignEvent = epochArray{mEpoch};
         if ~strcmp(alignEvent, 'stopSignalOn')  % No stop signals on go trials

            % Go to Target trials
            alignTimeList = trialData.(alignEvent)(iGoTargTrial);
            [targEEG, alignIndex] 	= align_signals(trialData.eegData(iGoTargTrial, kChannel), alignTimeList, cropWindow);
            satTrial                = eeg_reject_saturate(targEEG, 'alignIndex', alignIndex);
            targEEG(satTrial,:)     = [];
            targEEG                 = eeg_baseline_correct(targEEG, baseWindow, alignIndex);
            if filterData
               targEEGMean = lowpass(nanmean(targEEG, 1)', stopHz);
            else
               targEEGMean = nanmean(targEEG, 1);
            end
            
            Data(kChanIndex).signalStrength(iPropIndex).goTarg.(alignEvent).alignTime = alignIndex;
            if ~isempty(targEEG)
               Data(kChanIndex).signalStrength(iPropIndex).goTarg.(alignEvent).eeg = targEEG;
               Data(kChanIndex).signalStrength(iPropIndex).goTarg.(alignEvent).erp = targEEGMean;
               if ~strcmp(alignEvent, 'rewardOn') && ~strcmp(alignEvent, 'responseOnset')
                  eegMax(kChanIndex, mEpoch, iPropIndex, 1) = max(nanmean(targEEG, 1));
                  eegMin(kChanIndex, mEpoch, iPropIndex, 1) = min(nanmean(targEEG, 1));
               end
            else
               Data(kChanIndex).signalStrength(iPropIndex).goTarg.(alignEvent).eeg = [];
               Data(kChanIndex).signalStrength(iPropIndex).goTarg.(alignEvent).erp = [];
               Data(kChanIndex).signalStrength(iPropIndex).goTarg.(alignEvent).alignTime = [];
            end
            
            % Go to Distractor trials
            alignTimeList = trialData.(alignEvent)(iGoDistTrial);
            [distEEG, alignIndex]	= align_signals(trialData.eegData(iGoDistTrial, kChannel), alignTimeList, cropWindow);
            satTrial                = eeg_reject_saturate(distEEG, 'alignIndex', alignIndex);
            distEEG(satTrial,:)     = [];
            distEEG                 = eeg_baseline_correct(distEEG, baseWindow, alignIndex);
            if filterData
               distEEGMean = lowpass(nanmean(distEEG, 1)', stopHz);
            else
               distEEGMean = nanmean(distEEG, 1);
            end
            
            
            
            Data(kChanIndex).signalStrength(iPropIndex).goDist.(alignEvent).alignTime = alignIndex;
            if ~isempty(distEEG)
               Data(kChanIndex).signalStrength(iPropIndex).goDist.(alignEvent).eeg = distEEG;
               Data(kChanIndex).signalStrength(iPropIndex).goDist.(alignEvent).erp = distEEGMean;
            else
               Data(kChanIndex).signalStrength(iPropIndex).goDist.(alignEvent).eeg = [];
               Data(kChanIndex).signalStrength(iPropIndex).goDist.(alignEvent).erp = [];
               Data(kChanIndex).signalStrength(iPropIndex).goDist.(alignEvent).alignTime = [];
            end
         end
      end %mEpoch
      
      
      
      
      
      
      % Stop trials
      for jSSDIndex = 1 : length(ssdArray)
         jSSD = ssdArray(jSSDIndex);
         %         jSSD = 'all';
         
         jStopTargTrial = ccm_trial_selection(trialData,  {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'}, iPct, jSSD, targetHemifield);
         jStopDistTrial = ccm_trial_selection(trialData,  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'}, iPct, jSSD, targetHemifield);
         jStopCorrectTrial = ccm_trial_selection(trialData,  {'stopCorrect'}, iPct, jSSD, targetHemifield);
         
            
         for mEpoch = 1 : length(epochArray)
            alignEvent = epochArray{mEpoch};
            % Stop to Target trials
            alignTimeList = trialData.(alignEvent)(jStopTargTrial);
            [targEEG, alignIndex] 	= align_signals(trialData.eegData(jStopTargTrial, kChannel), alignTimeList, cropWindow);
            satTrial                = eeg_reject_saturate(targEEG, 'alignIndex', alignIndex);
            targEEG(satTrial,:)     = [];
            targEEG                 = eeg_baseline_correct(targEEG, baseWindow, alignIndex);
            if filterData
               targEEGMean = lowpass(nanmean(targEEG, 1)', stopHz);
            else
               targEEGMean = nanmean(targEEG, 1);
            end
            
            Data(kChanIndex).signalStrength(iPropIndex).stopTarg.ssd(jSSDIndex).(alignEvent).alignTime = alignIndex;
            if ~isempty(targEEG)
               Data(kChanIndex).signalStrength(iPropIndex).stopTarg.ssd(jSSDIndex).(alignEvent).eeg = targEEG;
               Data(kChanIndex).signalStrength(iPropIndex).stopTarg.ssd(jSSDIndex).(alignEvent).erp = targEEGMean;
            else
               Data(kChanIndex).signalStrength(iPropIndex).stopTarg.ssd(jSSDIndex).(alignEvent).eeg = [];
               Data(kChanIndex).signalStrength(iPropIndex).stopTarg.ssd(jSSDIndex).(alignEvent).erp = [];
               Data(kChanIndex).signalStrength(iPropIndex).stopTarg.ssd(jSSDIndex).(alignEvent).alignTime = [];
               
            end
            
            % Stop to Distractor trials
            alignTimeList = trialData.(alignEvent)(jStopDistTrial);
            [distEEG, alignIndex] 	= align_signals(trialData.eegData(jStopDistTrial, kChannel), alignTimeList, cropWindow);
            satTrial                = eeg_reject_saturate(distEEG, 'alignIndex', alignIndex);
            distEEG(satTrial,:)     = [];
            distEEG                 = eeg_baseline_correct(distEEG, baseWindow, alignIndex);
            if filterData
               distEEGMean = lowpass(nanmean(distEEG, 1)', stopHz);
            else
               distEEGMean = nanmean(distEEG, 1);
            end
            
            Data(kChanIndex).signalStrength(iPropIndex).stopDist.ssd(jSSDIndex).(alignEvent).alignTime = alignIndex;
            if ~isempty(distEEG)
               Data(kChanIndex).signalStrength(iPropIndex).stopDist.ssd(jSSDIndex).(alignEvent).eeg = distEEG;
               Data(kChanIndex).signalStrength(iPropIndex).stopDist.ssd(jSSDIndex).(alignEvent).erp = distEEGMean;
            else
               Data(kChanIndex).signalStrength(iPropIndex).stopDist.ssd(jSSDIndex).(alignEvent).eeg = [];
               Data(kChanIndex).signalStrength(iPropIndex).stopDist.ssd(jSSDIndex).(alignEvent).erp = [];
               Data(kChanIndex).signalStrength(iPropIndex).stopDist.ssd(jSSDIndex).(alignEvent).alignTime = [];
            end
            
            if ~strcmp(alignEvent, 'responseOnset')  % No response on correct Stop trials
               % Correct Stop trials
               alignTimeList = trialData.(alignEvent)(jStopCorrectTrial);
               [stopEEG, alignIndex] 	= align_signals(trialData.eegData(jStopCorrectTrial, kChannel), alignTimeList, cropWindow);
               satTrial                = eeg_reject_saturate(stopEEG, 'alignIndex', alignIndex);
               stopEEG(satTrial,:)     = [];
               stopEEG                 = eeg_baseline_correct(stopEEG, baseWindow, alignIndex);
            if filterData
               stopEEGMean = lowpass(nanmean(stopEEG, 1)', stopHz);
            else
               stopEEGMean = nanmean(stopEEG, 1);
            end
               
               Data(kChanIndex).signalStrength(iPropIndex).stopStop.ssd(jSSDIndex).(alignEvent).alignTime = alignIndex;
               if ~isempty(stopEEG)
                  Data(kChanIndex).signalStrength(iPropIndex).stopStop.ssd(jSSDIndex).(alignEvent).eeg = stopEEG;
                  Data(kChanIndex).signalStrength(iPropIndex).stopStop.ssd(jSSDIndex).(alignEvent).erp = stopEEGMean;
               else
                  Data(kChanIndex).signalStrength(iPropIndex).stopStop.ssd(jSSDIndex).(alignEvent).eeg = [];
                  Data(kChanIndex).signalStrength(iPropIndex).stopStop.ssd(jSSDIndex).(alignEvent).erp = [];
                  Data(kChanIndex).signalStrength(iPropIndex).stopStop.ssd(jSSDIndex).(alignEvent).alignTime = [];
               end
               
            end % jSSD
         end % mEpoch
         
      end % jSSDIndex
   end %iPropIndex
   
   
   
   
   
   
   
   
   
   
   
end % kUnitIndex
disp('completed gathering the data')



%%

PLOT_ERROR = true;
% epochArray = {'targOn', 'checkerOn', 'stopSignalOn', 'responseOnset', 'rewardOn'};

if plotFlag
   %     colorGrad = length(pSignalArray)/2;
   %     goC = zeros(colorGrad*2, 3);
   %     goC(1:colorGrad, 1:3) = repmat(linspace(0,.6,colorGrad)', 1, 3);
   %     stopC = goC;
   %     goC(1:colorGrad, 3) = linspace(.3,1,colorGrad);
   %     goC(colorGrad+1 : end, 1:3) = flipud(goC(1:colorGrad, :));
   %     stopC(1:colorGrad, 1) = linspace(.3,1,colorGrad);
   %     stopC(colorGrad+1 : end, 1:3) = flipud(stopC(1:colorGrad, :));
   
   colorGrad = length(pSignalArray);
   goC = zeros(colorGrad, 3);
   goC(:,1) = linspace(1,0,colorGrad);
   goC(:,2) = [linspace(0,0.5,colorGrad/2), linspace(0.5,0,colorGrad/2)];
   goC(:,3) = linspace(0,1,colorGrad);
   stopC = goC;
   
   cMap = ccm_colormap(pSignalArray);
   targLineW = 2;
   distLineW = 1;
   for kChanIndex = 1 : nElectrode
   if length(trialData.eegData(1,:)) < kChanIndex
      continue
   end
   
   
   nRow = 3;
      nColumn = nEpoch * 2 + 1;
      figureHandle = figureHandle + 1;
      if printPlot
         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
      else
         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
      end
      clf
      %         boxMargin = .5;
      %         x = xAxesPosition(end, 1);% - boxMargin;
      %         y = yAxesPosition(end, 1);% - boxMargin;sess
      %         w = axisWidth * nColumn/2;
      %         h = axisHeight * nRow;
      %                 rectangle('Position', [x, y, w, h], 'edgecolor', 'b')
      for mEpoch = 1 : nEpoch
         alignEvent = epochArray{mEpoch};
         epochRange = ccm_epoch_range(alignEvent, 'plot');
         
         yLimMax = nanmax(eegMax(kChanIndex, :)) * 1.1;
         yLimMin = nanmin(eegMin(kChanIndex, :)) * 1.1;
         
         % _______  Set up axes  ___________
         % axes names
         axGo = 1;
         axStopGo = 2;
         axStopStop = 3;
         
         % Set up plot axes
         % Left target Go trials
         ax(axGo, mEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(axGo, mEpoch) yAxesPosition(axGo, mEpoch) axisWidth axisHeight]);
         set(ax(axGo, mEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
         cla
         hold(ax(axGo, mEpoch), 'on')
         plot(ax(axGo, mEpoch), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
         title(epochArray{mEpoch})
         
         % Right target Go trials
         ax(axGo, mEpoch+nEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(axGo, mEpoch+nEpoch+1) yAxesPosition(axGo, mEpoch+nEpoch+1) axisWidth axisHeight]);
         set(ax(axGo, mEpoch+nEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
         cla
         hold(ax(axGo, mEpoch+nEpoch), 'on')
         plot(ax(axGo, mEpoch+nEpoch), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
         title(epochArray{mEpoch})
         %             set(ax(axRight, mEpoch), 'Xtick', [0 : 100 : epochRange(end) - epochRange(1)], 'XtickLabel', [epochRange(1) : 100: epochRange(end)])
         
         % Left target Stop Incorrect trials
         ax(axStopGo, mEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(axStopGo, mEpoch) yAxesPosition(axStopGo, mEpoch) axisWidth axisHeight]);
         set(ax(axStopGo, mEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
         cla
         hold(ax(axStopGo, mEpoch), 'on')
         plot(ax(axStopGo, mEpoch), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
         
         % Right target Stop Incorrect trials
         ax(axStopGo, mEpoch+nEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(axStopGo, mEpoch+nEpoch+1) yAxesPosition(axStopGo, mEpoch+nEpoch+1) axisWidth axisHeight]);
         set(ax(axStopGo, mEpoch+nEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
         cla
         hold(ax(axStopGo, mEpoch+nEpoch), 'on')
         plot(ax(axStopGo, mEpoch+nEpoch), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
         
         % Left target Stop Correct trials
         ax(axStopStop, mEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(axStopStop, mEpoch) yAxesPosition(axStopStop, mEpoch) axisWidth axisHeight]);
         set(ax(axStopStop, mEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
         cla
         hold(ax(axStopStop, mEpoch), 'on')
         plot(ax(axStopStop, mEpoch), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
         
         % Right target Stop Correct trials
         ax(axStopStop, mEpoch+nEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(axStopStop, mEpoch+nEpoch+1) yAxesPosition(axStopStop, mEpoch+nEpoch+1) axisWidth axisHeight]);
         set(ax(axStopStop, mEpoch+nEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
         cla
         hold(ax(axStopStop, mEpoch+nEpoch), 'on')
         plot(ax(axStopStop, mEpoch+nEpoch), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
         
         if mEpoch > 1
            set(ax(axGo, mEpoch), 'yticklabel', [])
            set(ax(axGo, mEpoch+nEpoch), 'yticklabel', [])
            set(ax(axStopGo, mEpoch), 'yticklabel', [])
            set(ax(axStopGo, mEpoch+nEpoch), 'yticklabel', [])
            set(ax(axStopStop, mEpoch), 'yticklabel', [])
            set(ax(axStopStop, mEpoch+nEpoch), 'yticklabel', [])
            
            set(ax(axGo, mEpoch), 'ycolor', [1 1 1])
            set(ax(axGo, mEpoch+nEpoch), 'ycolor', [1 1 1])
            set(ax(axStopGo, mEpoch), 'ycolor', [1 1 1])
            set(ax(axStopGo, mEpoch+nEpoch), 'ycolor', [1 1 1])
            set(ax(axStopStop, mEpoch), 'ycolor', [1 1 1])
            set(ax(axStopStop, mEpoch+nEpoch), 'ycolor', [1 1 1])
         end
         
         
         
         
         
         % __________ Loop signal strengths and plot  _________
         
         
         
         
         
         
         % PLOT LEFT TARGET TRIALS
         for i = 1 : nSignal/2
            iPropIndex = nSignal/2 + 1 - i;
            
            % Go trials
            if ~strcmp(alignEvent, 'stopSignalOn')  % No stop signals on go trials
               alignGoTarg = Data(kChanIndex).signalStrength(iPropIndex).goTarg.(alignEvent).alignTime;
               alignGoDist = Data(kChanIndex).signalStrength(iPropIndex).goDist.(alignEvent).alignTime;
               if ~isempty(alignGoTarg)
                  erpGoTarg = Data(kChanIndex).signalStrength(iPropIndex).goTarg.(alignEvent).erp;
                  plot(ax(axGo, mEpoch), epochRange, erpGoTarg(alignGoTarg + epochRange), 'color', cMap(iPropIndex,:), 'linewidth', targLineW)
               end
               if ~isempty(alignGoDist) && PLOT_ERROR
                  erpGoDist = Data(kChanIndex).signalStrength(iPropIndex).goDist.(alignEvent).erp;
                  plot(ax(axGo, mEpoch), epochRange, erpGoDist(alignGoDist + epochRange), '--', 'color', cMap(iPropIndex,:), 'linewidth', distLineW)
               end
            end
            
            
            % Stop signal trials
            stopTargEeg         = cell(1, length(ssdArray));
            stopTargAlign       = cell(1, length(ssdArray));
            stopDistEeg         = cell(1, length(ssdArray));
            stopDistAlign       = cell(1, length(ssdArray));
            stopStopEeg      = cell(1, length(ssdArray));
            stopStopAlign    = cell(1, length(ssdArray));
            for jSSDIndex = 1 : length(ssdArray)
               stopTargEeg{jSSDIndex}   = Data(kChanIndex).signalStrength(iPropIndex).stopTarg.ssd(jSSDIndex).(alignEvent).eeg;
               stopTargAlign{jSSDIndex} = Data(kChanIndex).signalStrength(iPropIndex).stopTarg.ssd(jSSDIndex).(alignEvent).alignTime;
               
               stopDistEeg{jSSDIndex}   = Data(kChanIndex).signalStrength(iPropIndex).stopDist.ssd(jSSDIndex).(alignEvent).eeg;
               stopDistAlign{jSSDIndex} = Data(kChanIndex).signalStrength(iPropIndex).stopDist.ssd(jSSDIndex).(alignEvent).alignTime;
               
               if ~strcmp(alignEvent, 'responseOnset')  % No stop signals on go trials
                  stopStopEeg{jSSDIndex} = Data(kChanIndex).signalStrength(iPropIndex).stopStop.ssd(jSSDIndex).(alignEvent).eeg;
                  stopStopAlign{jSSDIndex} = Data(kChanIndex).signalStrength(iPropIndex).stopStop.ssd(jSSDIndex).(alignEvent).alignTime;
               end
               
            end  % jSSDIndex = 1 : length(ssdArray)
            
            [eegStopTarg, alignStopTarg] = align_raster_sets(stopTargEeg, stopTargAlign);
            if filterData
               erpStopTarg = lowpass(nanmean(eegStopTarg, 1)', stopHz)';
            else
               erpStopTarg = nanmean(eegStopTarg, 1);
            end
            if size(erpStopTarg, 2) == 1, erpStopTarg = []; end;
            
            [eegStopDist, alignStopDist] = align_raster_sets(stopDistEeg, stopDistAlign);
            if filterData
               erpStopDist = lowpass(nanmean(eegStopDist, 1)', stopHz)';
            else
               erpStopDist = nanmean(eegStopDist, 1);
            end
            if size(erpStopDist, 2) == 1, erpStopDist = []; end;
            
            if ~strcmp(alignEvent, 'responseOnset')  % No stop signals on go trials
               [eegStopCorrect, alignStopCorrect] = align_raster_sets(stopStopEeg, stopStopAlign);
             if filterData
               erpStopCorrect = lowpass(nanmean(eegStopCorrect, 1)', stopHz)';
            else
               erpStopCorrect = nanmean(eegStopCorrect, 1);
            end
               if size(erpStopCorrect, 2) == 1, erpStopCorrect = []; end;
            end
            
            
            
            
            if ~isempty(erpStopTarg)
               plot(ax(axStopGo, mEpoch), epochRange, erpStopTarg(alignStopTarg + epochRange), 'color', cMap(iPropIndex,:), 'linewidth', targLineW)
            end
            if PLOT_ERROR && ~isempty(erpStopDist)
               plot(ax(axStopGo, mEpoch), epochRange, erpStopDist(alignStopDist + epochRange), '--', 'color', cMap(iPropIndex,:), 'linewidth', distLineW)
            end
            
            
            if ~strcmp(alignEvent, 'responseOnset')  % No stop signals on go trials
               if ~isempty(erpStopCorrect)
                  plot(ax(axStopStop, mEpoch), epochRange, erpStopCorrect(alignStopCorrect + epochRange), 'color', cMap(iPropIndex,:), 'linewidth', targLineW)
               end
            end
            
            
            
            
            
            
            
            
            % PLOT RIGHT TARGET TRIALS
            
            iPropIndexR = i + nSignal/2;  % Reverse order of plotting to keep color overlays similar between left and right target
            
            
            % Go trials
            if ~strcmp(alignEvent, 'stopSignalOn')  % No stop signals on go trials
               alignGoTarg = Data(kChanIndex).signalStrength(iPropIndexR).goTarg.(alignEvent).alignTime;
               alignGoDist = Data(kChanIndex).signalStrength(iPropIndexR).goDist.(alignEvent).alignTime;
               erpGoTarg = Data(kChanIndex).signalStrength(iPropIndexR).goTarg.(alignEvent).erp;
               erpGoDist = Data(kChanIndex).signalStrength(iPropIndexR).goDist.(alignEvent).erp;
               
               if ~isempty(erpGoTarg)
                  plot(ax(axGo, mEpoch + nEpoch), epochRange, erpGoTarg(alignGoTarg + epochRange), 'color', cMap(iPropIndexR,:), 'linewidth', targLineW)
               end
               if PLOT_ERROR && ~isempty(erpGoDist)
                  plot(ax(axGo, mEpoch + nEpoch), epochRange, erpGoDist(alignGoDist + epochRange), '--', 'color', cMap(iPropIndexR,:), 'linewidth', distLineW)
               end
            end
            
            
            
            % Stop signal trials
            stopTargEeg     = cell(1, length(ssdArray));
            stopTargAlign   = cell(1, length(ssdArray));
            stopDistEeg     = cell(1, length(ssdArray));
            stopDistAlign   = cell(1, length(ssdArray));
            stopStopEeg  = cell(1, length(ssdArray));
            stopStopAlign = cell(1, length(ssdArray));
            for jSSDIndex = 1 : length(ssdArray)
               stopTargEeg{jSSDIndex}   = Data(kChanIndex).signalStrength(iPropIndexR).stopTarg.ssd(jSSDIndex).(alignEvent).eeg;
               stopTargAlign{jSSDIndex} = Data(kChanIndex).signalStrength(iPropIndexR).stopTarg.ssd(jSSDIndex).(alignEvent).alignTime;
               
               if PLOT_ERROR
                  stopDistEeg{jSSDIndex}    = Data(kChanIndex).signalStrength(iPropIndexR).stopDist.ssd(jSSDIndex).(alignEvent).eeg;
                  stopDistAlign{jSSDIndex}  = Data(kChanIndex).signalStrength(iPropIndexR).stopDist.ssd(jSSDIndex).(alignEvent).alignTime;
               end
               
               if ~strcmp(alignEvent, 'responseOnset')  % No stop signals on go trials
                  stopStopEeg{jSSDIndex}         = Data(kChanIndex).signalStrength(iPropIndexR).stopStop.ssd(jSSDIndex).(alignEvent).eeg;
                  stopStopAlign{jSSDIndex}       = Data(kChanIndex).signalStrength(iPropIndexR).stopStop.ssd(jSSDIndex).(alignEvent).alignTime;
                  [eegStopCorrect, alignStopCorrect] = align_raster_sets(stopStopEeg, stopStopAlign);
                  erpStopCorrect                    = nanmean(eegStopCorrect, 1);
               end
               
            end  % jSSDIndex = 1 : length(ssdArray)
            
            
            [eegStopTarg, alignStopTarg]    = align_raster_sets(stopTargEeg, stopTargAlign);
            if filterData
               erpStopTarg = lowpass(nanmean(eegStopTarg, 1)', stopHz)';
            else
               erpStopTarg = nanmean(eegStopTarg, 1);
            end
            if size(erpStopTarg, 2) == 1, erpStopTarg = []; end;
            
            [eegStopDist, alignStopDist]    = align_raster_sets(stopDistEeg, stopDistAlign);
            if filterData
               erpStopDist = lowpass(nanmean(eegStopDist, 1)', stopHz)';
            else
               erpStopDist = nanmean(eegStopDist, 1);
            end
            if size(erpStopDist, 2) == 1, erpStopDist = []; end;
            
            if ~strcmp(alignEvent, 'responseOnset')  % No stop signals on go trials
               [eegStopCorrect, alignStopCorrect]   = align_raster_sets(stopStopEeg, stopStopAlign);
             if filterData
               erpStopCorrect = lowpass(nanmean(eegStopCorrect, 1)', stopHz)';
            else
               erpStopCorrect = nanmean(eegStopCorrect, 1);
            end
              if size(erpStopCorrect, 2) == 1, erpStopCorrect = []; end;
            end
            
            
            if ~isempty(erpStopTarg)
               plot(ax(axStopGo, mEpoch + nEpoch), epochRange, erpStopTarg(alignStopTarg + epochRange), 'color', cMap(iPropIndexR,:), 'linewidth', targLineW)
            end
            if PLOT_ERROR && ~isempty(erpStopDist)
               plot(ax(axStopGo, mEpoch + nEpoch), epochRange, erpStopDist(alignStopDist + epochRange), '--', 'color', cMap(iPropIndexR,:), 'linewidth', distLineW)
            end
            if ~strcmp(alignEvent, 'responseOnset')  % No stop signals on go trials
               if ~isempty(erpStopCorrect)
                  plot(ax(axStopStop, mEpoch + nEpoch), epochRange, erpStopCorrect(alignStopCorrect + epochRange), 'color', cMap(iPropIndexR,:), 'linewidth', targLineW)
               end
            end
         end % iPropIndex
         %                       if mEpoch == 3
         %                 colormap(goC);
         %                 legend(ax(axGo, 1), '
         
      end % mEpoch
      
      %                             legend(ax(axGo, 1), {num2cell(pSignalArray'), num2str(pSignalArray')})
      
      %         colorbar('peer', ax(axGo, 1), 'location', 'west')
      %         colorbar('peer', ax(axStopGo, 1), 'location', 'west')
      h=axes('Position', [0 0 1 1], 'Visible', 'Off');
      titleString = sprintf('%s \t %s', sessionID, Data(kChanIndex).name);
      text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')
      if printPlot
         localFigurePath = local_figure_path;
         print(figureHandle,[localFigurePath, sessionID, '_', unitArray{kChanIndex}, '_ccm_single_neuron'],'-dpdf', '-r300')
      end
   end % kUnitIndex
   
   
end % plotFlag





end

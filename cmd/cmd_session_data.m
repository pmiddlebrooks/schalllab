function Data = cmd_session_data(subjectID, sessionID, options)%dataType, varargin)

%
% function Data = cmd_session_data(subjectID, sessionID, options)%dataType, varargin)
%
% Single neuron analyses for choice countermanding task. Only plots the
% sdfs. To see rasters, use ccm_single_neuron_rasters, which displays all
% conditions in a given epoch
%
% input:
%   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%
%   options: A structure with various ways to select/organize data: If
%   ccm_session_data.m is called without input arguments, the default
%   options structure is returned. options has the following fields with
%   possible values (default listed first):
%
%    options.dataType = 'neuron', 'lfp', 'erp';
%
%    options.figureHandle   = 1000;
%    options.printPlot      = false, true;
%    options.plotFlag       = true, false;
%    options.doStops        = true, false;
%    options.filterData 	= false, true;
%    options.stopHz         = 50, <any number, above which signal is filtered;
%    options.normalize      = false, true;
%
%
% Returns Data structure with fields:
%
%   Data.angle(x).(condition).ssd(x).(epoch name).(signal)

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
if nargin < 3
   options.dataType = 'neuron';
   
   options.figureHandle     = 1000;
   options.printPlot        = false;
   options.plotFlag         = true;
   options.doStops          = true;
   options.filterData       = false;
   options.stopHz           = 50;
   options.normalize        = false;
   options.targAngle        = 'each';
   
   %    options.xxxx = xxxx;
   %    options.xxxx = xxxx;
   %    options.xxxx = xxxx;
   if nargin == 0
      Data = options;
      return
   end
end
doStops = options.doStops;
normalize = options.normalize;



clear Data
Data    = [];


% Load the data
switch options.dataType
    case 'neuron'
        addVar = 'spikeData';
    case 'lfp'
        addVar = 'lfpData';
    case 'erp'
        addVar = 'eegData';
end
variables = [ccm_min_vars, addVar, 'trialDuration'];

[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID, variables);
ssdArray        = ExtraVar.ssdArray;
nSSD            = length(ssdArray);
targAngleArray	= ExtraVar.targAngleArray;
nTarg           = length(targAngleArray);

if ~strcmp(SessionData.taskID, 'cmd')
   fprintf('Not a simple countermanding session, try again\n')
   return
end



% CONSTANTS
MIN_RT          = 120;
MAX_RT          = 1200;
STD_MULTIPLE    = 3;

Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;


cropWindow  = -499 : 500;  % used to extract a semi-small portion of signal for each epoch/alignemnt
baseWindow 	= -149 : 0;   % To baseline-shift the eeg signals, relative to event alignment index;








% Set defaults
dataType = options.dataType;
switch dataType
   case 'neuron'
      if ~isfield(SessionData, 'spikeUnitArray')
         spikeUnitArray = {};
         for i = 1 : size(trialData.trialOutcome.spikeData, 2)
            spikeUnitArray = [spikeUnitArray, num2str(i)];
            %    SessionData.spikeUnitArray = num2str([1:size(trialData.trialOutcome.spikeData,2)]);
         end
         SessionData.spikeUnitArray = spikeUnitArray;
      end
      dataArray     = SessionData.spikeUnitArray;
   case 'lfp'
      dataArray 	= num2cell(SessionData.lfpChannel);
   case 'erp'
      dataArray     = eeg_electrode_map(subjectID);
end
% Make sure user input a dataType that was recorded during the session
dataTypePossible = {'neuron', 'lfp', 'erp'};
if ~sum(strcmp(dataType, dataTypePossible))
   fprintf('%s Is not a valid data type \n', dataType)
   return
end
if isempty(dataArray)
   fprintf('Session %s apparently does not contain %s data \n', sessionID, dataType)
   return
end

if strcmp(options.targAngle, 'collapse')
   nTarg = 1;
end




% epochArray = {'fixWindowEntered', 'targOn', 'checkerOn', 'stopSignalOn', 'responseOnset', 'rewardOn'};
epochArray = {'fixWindowEntered', 'targOn', 'stopSignalOn', 'responseOnset', 'rewardOn'};
nEpoch = length(epochArray);
nOutcome = 3; % Used to find the maximum signal levels for normalization if desired





% How many units were recorded?
nUnit = size(dataArray, 2);



% Get rid of trials with outlying RTs
[allRT, rtOutlierTrial] = truncate_rt(trialData.rt, MIN_RT, MAX_RT, STD_MULTIPLE);
trialData.rt(rtOutlierTrial) = nan;




% If there weren't stop trials, skip all stop-related analyses
if isempty(ssdArray) || ~doStops
   disp('cmd_session_data.m: No stop trials or stop trial analyses not requested');
end





[minTrialDuration,i] = min(trialData.trialDuration(~strcmp(trialData.trialOutcome, 'fixationAbort')));
[maxTrialDuration,i] = max(trialData.trialDuration(~strcmp(trialData.trialOutcome, 'fixationAbort')));

for kDataIndex = 1 : nUnit
   switch dataType
      case 'neuron'
         [a, kUnit] = ismember(dataArray{kDataIndex}, SessionData.spikeUnitArray);
      case 'lfp'
         [a, kUnit] = ismember(dataArray{kDataIndex}, SessionData.lfpChannel);
      case 'erp'
         [a, kUnit] = ismember(dataArray{kDataIndex}, eeg_electrode_map(subjectID));
   end
   
   yMax = zeros(nEpoch, nTarg, nSSD+1, nOutcome);  % Keep track of maximum signal values, for setting y-axis limits in plots
   yMin = zeros(nEpoch, nTarg, nSSD+1, nOutcome);  % Keep track of maximum signal values, for setting y-axis limits in plots
   
   Data(kDataIndex).subjectID       = subjectID;
   Data(kDataIndex).sessionID       = sessionID;
   Data(kDataIndex).name            = dataArray{kDataIndex};
   Data(kDataIndex).ssdArray        = ssdArray;
   Data(kDataIndex).targAngleArray  = targAngleArray;
   
   
   
   
   % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % LOOP THROUGH ANGLES
   % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for iTarg = 1 : nTarg;
      
      
      if strcmp(options.targAngle, 'collapse')
         selectOpt.targAngle = 'all';
      elseif strcmp(options.targAngle, 'each')
         selectOpt.targAngle = targAngleArray(iTarg);
      else
         selectOpt.targAngle = options.targAngle;
      end
      
      
      
      
      
      % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %              GO TRIALS
      % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      selectOpt.ssdRange = 'none';
      
      % Select appropriate trials to analyzie
      selectOpt.outcomeArray = {'goCorrectTarget'; 'goCorrect'; 'targetHoldAbort'};
      iGoTargTrial = cmd_trial_selection(trialData, selectOpt);
      
      
      for mEpoch = 1 : length(epochArray)
         mEpochName = epochArray{mEpoch};
         if ~strcmp(mEpochName, 'stopSignalOn')  % No stop signals on go trials
            
            alignListGoTarg = trialData.(mEpochName)(iGoTargTrial);
            Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).alignTimeList = alignListGoTarg;   % Keep track of trial-by-trial alignemnt
            
            
            switch dataType
               case 'neuron'
                  % Go to Target trials
                  [alignedRasters, alignmentIndex] = spike_to_raster(trialData.spikeData(iGoTargTrial, kUnit), alignListGoTarg);
                  Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).alignTime = alignmentIndex;
                  sdf = spike_density_function(alignedRasters, Kernel);
                  if ~isempty(sdf); yMax(mEpoch, iTarg, 1, 1) = nanmax(nanmean(sdf, 1)); end;
                  Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).raster = alignedRasters;
                  Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).sdf = sdf;
                  Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).sdfMean = nanmean(sdf, 1);
                  
                  
               case 'lfp'
                  % Go to Target trials
                  [targLFP, alignIndex] 	= align_signals(trialData.lfpData(iGoTargTrial, kUnit), alignListGoTarg, cropWindow);
                  satTrial                = signal_reject_saturate(targLFP, 'alignIndex', alignIndex);
                  targLFP(satTrial,:)     = [];
                  targLFP                 = signal_baseline_correct(targLFP, baseWindow, alignIndex);
                  if filterData
                     targLFPMean = lowpass(nanmean(targLFP, 1)', options.stopHz);
                  else
                     targLFPMean = nanmean(targLFP, 1);
                  end
                  
                  Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).alignTime = alignIndex;
                  Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).lfp = targLFP;
                  Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).lfpMean = targLFPMean;
                  if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                     if ~isempty(targLFP)
                        yMax(mEpoch, iTarg, 1, 1) = nanmax(nanmean(targLFP, 1));
                        yMin(mEpoch, iTarg, 1, 1) = nanmin(nanmean(targLFP, 1));
                     end
                  end
                  
                  
                  
                  
               case 'erp'
                  % Go to Target trials
                  [targEEG, alignIndex] 	= align_signals(trialData.eegData(iGoTargTrial, kUnit), alignListGoTarg, cropWindow);
                  satTrial                = signal_reject_saturate(targEEG, 'alignIndex', alignIndex);
                  targEEG(satTrial,:)     = [];
                  targEEG                 = signal_baseline_correct(targEEG, baseWindow, alignIndex);
                  if filterData
                     targEEGMean = lowpass(nanmean(targEEG, 1)', options.stopHz);
                  else
                     targEEGMean = nanmean(targEEG, 1);
                  end
                  
                  Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).alignTime = alignIndex;
                  Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).eeg = targEEG;
                  Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).erp = targEEGMean;
                  if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                     if ~isempty(targEEG)
                        yMax(kDataIndex, mEpoch, iTarg, 1) = nanmax(nanmean(targEEG, 1));
                        yMin(kDataIndex, mEpoch, iTarg, 1) = nanmin(nanmean(targEEG, 1));
                     end
                  end
                  
                  
            end % switch dataType
         end % ~strcmp(mEpochName, 'stopSignalOn')
      end % mEpoch
      
      
      
      
      
      
      
      % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Stop trials
      % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      
      for jSSDIndex = 1 : nSSD
         jSSD = ssdArray(jSSDIndex);
         selectOpt.ssdRange = jSSD;
         
         selectOpt.outcomeArray = {'stopIncorrectTarget'; 'stopIncorrect'; 'targetHoldAbort'; 'stopIncorrectPreSSDTarget'};
         jStopTargTrial = cmd_trial_selection(trialData, selectOpt);
         selectOpt.outcomeArray = {'stopCorrect'};
         jStopCorrectTrial = cmd_trial_selection(trialData, selectOpt);
         
         for mEpoch = 1 : length(epochArray)
            mEpochName = epochArray{mEpoch};
            
            alignListStopTarg = trialData.(mEpochName)(jStopTargTrial);
            Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).alignTimeList = alignListStopTarg;   % Keep track of trial-by-trial alignemnt
            alignListStopCorrect = trialData.(mEpochName)(jStopCorrectTrial);
            Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).alignTimeList = alignListStopCorrect;   % Keep track of trial-by-trial alignemnt
            
            
            
            switch dataType
               case 'neuron'
                  % Stop to Target trials
                  [alignedRasters, alignmentIndex] = spike_to_raster(trialData.spikeData(jStopTargTrial, kUnit), alignListStopTarg);
                  Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).alignTime = alignmentIndex;
                  
                  sdf = spike_density_function(alignedRasters, Kernel);
                  if ~isempty(sdf); yMax(mEpoch, iTarg, jSSDIndex+1, 2) = nanmax(nanmean(sdf, 1)); end;
                  Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).raster = alignedRasters;
                  Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).sdf = sdf;
                  Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).sdfMean = nanmean(sdf, 1);
                  
                  if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
                     % Stop to Target trials
                     [alignedRasters, alignmentIndex] = spike_to_raster(trialData.spikeData(jStopCorrectTrial, kUnit), alignListStopCorrect);
                     Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).alignTime = alignmentIndex;
                     
                     sdf = spike_density_function(alignedRasters, Kernel);
                     if ~isempty(sdf); yMax(mEpoch, iTarg, jSSDIndex+1, 3) = nanmax(nanmean(sdf, 1)); end;
                     Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).raster = alignedRasters;
                     Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).sdf = sdf;
                     Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).sdfMean = nanmean(sdf, 1);
                     
                  end
                  
                  
                  
                  
               case 'lfp'
                  
                  % Stop to Target trials
                  [targLFP, alignIndex] 	= align_signals(trialData.lfpData(jStopTargTrial, kUnit), alignListStopTarg, cropWindow);
                  satTrial                = signal_reject_saturate(targLFP, 'alignIndex', alignIndex);
                  targLFP(satTrial,:)     = [];
                  targLFP                 = signal_baseline_correct(targLFP, baseWindow, alignIndex);
                  if filterData
                     targLFPMean = lowpass(nanmean(targLFP, 1)', options.stopHz);
                  else
                     targLFPMean = nanmean(targLFP, 1);
                  end
                  
                  Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).alignTime = alignIndex;
                  Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).lfp = targLFP;
                  Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).lfpMean = targLFPMean;
                  if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                     % if ~isempty(targLFP)
                     %                         yMax(kDataIndex, mEpoch, iPropIndex, jSSDIndex + 1, 3) = nanmax(nanmean(targLFP, 1));
                     %                         yMin(kDataIndex, mEpoch, iPropIndex, jSSDIndex + 1, 3) = nanmin(nanmean(targLFP, 1));
                     % end
                  end
                  
                  
                  
                  if ~strcmp(mEpochName, 'responseOnset')  % No response on correct Stop trials
                     % Correct Stop trials
                     [stopLFP, alignIndex] 	= align_signals(trialData.lfpData(jStopCorrectTrial, kUnit), alignListStopCorrect, cropWindow);
                     satTrial                = signal_reject_saturate(stopLFP, 'alignIndex', alignIndex);
                     stopLFP(satTrial,:)     = [];
                     stopLFP                 = signal_baseline_correct(stopLFP, baseWindow, alignIndex);
                     if filterData
                        stopLFPMean = lowpass(nanmean(stopLFP, 1)', options.stopHz);
                     else
                        stopLFPMean = nanmean(stopLFP, 1);
                     end
                     
                     Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).alignTime = alignIndex;
                     Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).lfp = stopLFP;
                     Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).lfpMean = stopLFPMean;
                     if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                        if ~isempty(stopLFP)
                           yMax(kDataIndex, mEpoch, iTarg, jSSDIndex + 1, 3) = nanmax(nanmean(stopLFP, 1));
                           yMin(kDataIndex, mEpoch, iTarg, jSSDIndex + 1, 3) = nanmin(nanmean(stopLFP, 1));
                        end
                     end
                     
                  end % jSSD
                  
                  
                  
               case 'erp'
                  % Stop to Target trials
                  [targEEG, alignIndex] 	= align_signals(trialData.eegData(jStopTargTrial, kUnit), alignListStopTarg, cropWindow);
                  satTrial                = signal_reject_saturate(targEEG, 'alignIndex', alignIndex);
                  targEEG(satTrial,:)     = [];
                  targEEG                 = signal_baseline_correct(targEEG, baseWindow, alignIndex);
                  if filterData
                     targEEGMean = lowpass(nanmean(targEEG, 1)', options.stopHz);
                  else
                     targEEGMean = nanmean(targEEG, 1);
                  end
                  
                  Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).alignTime = alignIndex;
                  Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).eeg = targEEG;
                  Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).erp = targEEGMean;
                  if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                     % if ~isempty(targEEG)
                     %                         yMax(kDataIndex, mEpoch, iPropIndex, jSSDIndex + 1, 3) = nanmax(nanmean(targEEG, 1));
                     %                         yMin(kDataIndex, mEpoch, iPropIndex, jSSDIndex + 1, 3) = nanmin(nanmean(targEEG, 1));
                     %                      end
                  end
                  
                  
                  
                  
                  if ~strcmp(mEpochName, 'responseOnset')  % No response on correct Stop trials
                     % Correct Stop trials
                     [stopEEG, alignIndex] 	= align_signals(trialData.eegData(jStopCorrectTrial, kUnit), alignListStopCorrect, cropWindow);
                     satTrial                = signal_reject_saturate(stopEEG, 'alignIndex', alignIndex);
                     stopEEG(satTrial,:)     = [];
                     stopEEG                 = signal_baseline_correct(stopEEG, baseWindow, alignIndex);
                     if filterData
                        stopEEGMean = lowpass(nanmean(stopEEG, 1)', options.stopHz);
                     else
                        stopEEGMean = nanmean(stopEEG, 1);
                     end
                     
                     Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).alignTime = alignIndex;
                     Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).eeg = stopEEG;
                     Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).erp = stopEEGMean;
                     if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                        if ~isempty(stopEEG)
                           yMax(mEpoch, iTarg, jSSDIndex + 1, 3) = nanmax(nanmean(stopEEG, 1));
                           yMin(mEpoch, iTarg, jSSDIndex + 1, 3) = nanmin(nanmean(stopEEG, 1));
                        end
                     end
                  end
                  
            end % switch dataType
         end % mEpoch
         
      end % jSSD
      
      
   end %iPropIndex
   
   
   
   
   
   
   
   
   a = yMax(:,:,1,:);   % the go trials
   b = yMax(:,2:end,:);  % all the stop trials
   c = yMin(:,:,1,:);   % go trials
   d = yMin(:,2:end,:);  % stop trials
   
   Data(kDataIndex).yMax = max([a(:); nanmean(b(:))]);
   Data(kDataIndex).yMin = min([c(:); nanmean(d(:))]);
   %    Data(kDataIndex).yMax = max([yMax(:,:,1,:); nanmean(nanmax(yMax(:,:,2:end,:)))]);
   %    Data(kDataIndex).yMin = max([yMin(:,:,1,:); nanmean(nanmax(yMin(:,:,2:end,:)))]);
   %    Data(kDataIndex).yMin = nanmin(yMin(:));
   
end % kDataIndex










% If we want the signal to be normalized to the maximum value for each
% unit/channel, loop back through and normalize here
if normalize
   for kDataIndex = 1 : nUnit
      for iTarg = 1 : nSignal;
         for mEpoch = 1 : length(epochArray)
            mEpochName = epochArray{mEpoch};
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Go trials
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~strcmp(mEpochName, 'stopSignalOn')  % No stop signals on go trials
               
               switch dataType
                  case 'neuron'
                     % Go to Target trials
                     Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).sdf = ...
                        Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).sdf ./ Data(kDataIndex).yMax;
                     Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).sdfMean = ...
                        nanmean(Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).sdf, 1);
                     
                     
                  case 'lfp'
                     % Go to Target trials
                     Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).lfp = ...
                        Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).lfp ./ Data(kDataIndex).yMax;
                     Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).lfpMean = ...
                        nanmean(Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).lfp, 1);
                     
                     
                  case 'erp'
                     % Go to Target trials
                     Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).eeg = ...
                        Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).eeg ./ Data(kDataIndex).yMax;
                     Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).lfpMean = ...
                        nanmean(Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).eeg, 1);
                     
               end % switch dataType
            end % ~strcmp(mEpochName, 'stopSignalOn')
            
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Stop trials
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for jSSDIndex = 1 : nSSD
               
               switch dataType
                  case 'neuron'
                     % Stop to Target trials
                     Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).sdf = ...
                        Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).sdf ./ Data(kDataIndex).yMax;
                     Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).sdfMean = ...
                        nanmean(Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).sdf, 1);
                     
                     
                     if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
                        
                        % Stop Correct trials
                        Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).sdf = ...
                           Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).sdf ./ Data(kDataIndex).yMax;
                        Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).sdfMean = ...
                           nanmean(Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).sdf, 1);
                        
                     end % ~strcmp(mEpochName, 'responseOnset')
                     
                     
                  case 'lfp'
                     
                     % Stop to Target trials
                     Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).lfp = ...
                        Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).lfp ./ Data(kDataIndex).yMax;
                     Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).lfpMean = ...
                        nanmean(Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).lfp, 1);
                     
                     
                     if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
                        
                        % Stop Correct trials
                        Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).lfp = ...
                           Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).lfp ./ Data(kDataIndex).yMax;
                        Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).lfpMean = ...
                           nanmean(Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).lfp, 1);
                        
                     end % ~strcmp(mEpochName, 'responseOnset')
                     
                     
                  case 'erp'
                     
                     % Stop to Target trials
                     Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).eeg = ...
                        Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).eeg ./ Data(kDataIndex).yMax;
                     Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).erp = ...
                        nanmean(Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).eeg, 1);
                     
                     if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
                        
                        % Stop Correct trials
                        Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).eeg = ...
                           Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).eeg ./ Data(kDataIndex).yMax;
                        Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).erp = ...
                           nanmean(Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).eeg, 1);
                        
                     end % ~strcmp(mEpochName, 'responseOnset')
                     
               end % switch dataType
            end % jSSD
         end % mEpoch
      end % for mEpoch = 1 : length(epochArray)
      Data(kDataIndex).yMax = 1.1;
      switch dataType
         case 'neuron'
            Data(kDataIndex).yMin = 0;
         case {'lfp', 'erp'}
            Data(kDataIndex).yMin = -1.1;
      end
   end % for kDataIndex = 1 : nUnit
   
   
end % if normalize



disp('completed data collection')

Data(1).dataArray       = dataArray;
Data(1).ssdArray        = ssdArray;
Data(1).targAngleArray  = targAngleArray;
Data(1).sessionID       = sessionID;
Data(1).subjectID       = subjectID;



if options.plotFlag
   cmd_session_data_plot(Data, options)
   
end



function Data = ccm_session_data(subjectID, sessionID, Opt)
%
% function Data = ccm_session_data(subjectID, sessionID, Opt)
%
% Creates a processed data struct for other functions to use, and
% optionally plots session neurophysiology data
%
%
% Single neuron analyses for choice countermanding task. Only plots the
% sdfs. To see rasters, use ccm_single_neuron_rasters, which displays all
% conditions in a given epoch
%
% input:
%   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%
%   Opt: A structure with various ways to select/organize data:
%   The following fields of Opt are relevant for ccm_session_data, with
%   possible values (default listed first):
%
%    Opt.trialData: if Options structure contains trialData and
%    SessionData, don't need to load the data
%   Opt.SessionData: ditto
%    Opt.dataType = 'neuron', 'lfp', 'erp';
%
%    Opt.figureHandle   = 1000;
%    Opt.printPlot      = false, true;
%    Opt.plotFlag       = true, false;
%    Opt.collapseSignal = false, true;
%    Opt.collapseTarg 	= false, true;
%    Opt.doStops        = true, false;
%    Opt.filterData 	= false, true;
%    Opt.stopHz         = 50, <any number, above which signal is filtered;
%    Opt.normalize      = false, true;
%    Opt.howProcess      = how to step through the list of units
%                                 'each' to plot all,
%                                 'step' (default): step through to see one
%                                 plot at a time, pausing between
%                                 'print' (default): step through each plot
%                                 individually, printing each to file
%                                   'pop'- plot the population of units
%                                   from the session alog the depth of the
%                                   electrode
%    Opt.unitArray      = units you want to analyze (default gets filled
%                                   with all available).
%                                 {'spikeUnit17a'}, input a cell of units for a list of individaul units
%    Opt.baselineCorrect = false, true; Baseline correct analog signals?
%    Opt.hemisphere = which hemsiphere were the data recorded from? left or
%    right?
%
%
%
% Returns Data structure with fields:
%
%   Data.(epoch).colorCoherence(x).(outcome).ssd(x)
%
%   outcome can be:  goTarg, goDist, stopTarg, stopDist, stopStop
%   ssd(x):  only applies for stop trials, else the field is absent
%   epoch name: fixOn, targOn, checkerOn, etc.

%%
clear Data

if nargin < 3
    Opt = ccm_options;
    Opt.trialData = [];
end

% User may have input the data already, so you don't need to load it
if iscell(sessionID)
    if nargin < 3
        Opt = ccm_options;
        Opt.multiUnit = false;
        Opt.printFlag = 0;
        Opt.plotFlag = 0;
        Opt.doStops = 1;
        Opt.collapseSignal = 0;
        Opt.collapseTarg = 1;
    end
    Opt.unitArray = sessionID(2);
    sessionID = sessionID{1};
end

% Arrays for looping through data
epochArray          = {'fixWindowEntered', 'targOn', 'checkerOn', 'stopSignalOn', 'responseOnset', 'toneOn', 'rewardOn'};
goOutcomeArray      = {'goTarg', 'goDist'};
stopOutcomeArray    = {'stopTarg', 'stopDist', 'stopStop'};
nEpoch              = length(epochArray);
nOutcome            = length(goOutcomeArray) + length(stopOutcomeArray);


% tic
% Load (or get the input) data
if ~isempty(Opt.trialData)
    trialData = Opt.trialData;
    SessionData = Opt.SessionData;
    ExtraVar = Opt.ExtraVar;
else
    if strcmp(Opt.dataType, 'neuron')
        variables = [ccm_min_vars, Opt.unitArray];
        [trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID, variables, Opt.multiUnit);
        if strcmp(Opt.unitArray, 'spikeData')
            Opt.unitArray = SessionData.spikeUnitArray;
        end
    else
        variables = [ccm_min_vars, 'lfpData'];
        [trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID, variables);
    end
    pSignalArray    = ExtraVar.pSignalArray;
    targAngleArray	= ExtraVar.targAngleArray;
    ssdArray        = ExtraVar.ssdArray;
    nSSD            = length(ssdArray);
end
% fprintf('load time: %.2f\n', toc)

if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a chioce countermanding saccade session, try again\n')
    return
end



% CONSTANTS
MIN_RT          = 120;
MAX_RT          = 1200;
STD_MULTIPLE    = 3;

Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;
% Kernel.method = 'gaussian';
% Kernel.sigma = 10;


cropWindow  = -1000 : 2000;  % used to extract a semi-small portion of signal for each epoch/alignemnt
baseWindow 	= -149 : 0;   % To baseline-shift the eeg signals, relative to event alignment index;





% How many units were recorded?
nUnit = length(Opt.unitArray);




% Make sure user input a dataType that was recorded during the session
dataTypePossible = {'neuron', 'lfp', 'erp'};
if ~sum(strcmp(Opt.dataType, dataTypePossible))
    fprintf('%s Is not a valid data type \n', Opt.dataType)
    return
end
if isempty(Opt.unitArray)
    fprintf('Session %s apparently does not contain %s data \n', sessionID, Opt.dataType)
    return
end


% If collapsing into all left and all right need to note here that there are "2" angles to deal with
% (important for calling ccm_trial_selection.m)
leftTargInd = (targAngleArray < -89) & (targAngleArray > -270) | ...
    (targAngleArray > 90) & (targAngleArray < 269);
rightTargInd = ~leftTargInd;
if Opt.collapseTarg
    nTargPair = 1;
else
    nTargPair = sum(rightTargInd);
    % do nothing, all target angles will be considered separately
end




% Get rid of trials with outlying RTs
[allRT, rtOutlierTrial] = truncate_rt(trialData.rt, MIN_RT, MAX_RT, STD_MULTIPLE);
trialData.rt(rtOutlierTrial) = nan;

trialData = ccm_delete_nan_rt(trialData);


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
% If collapsing data across signal strength, adjust the pSignalArray here
if Opt.collapseSignal
    nSignal = 2;
end





% If there weren't stop trials, skip all stop-related analyses
if isempty(ssdArray) || ~Opt.doStops
    %     disp('ccm_inhibition.m: No stop trials or stop trial analyses not requested');
    Opt.doStops = false;
end





% Get ssrt for the session (to be used for latency matching below)
optInh                  = ccm_options;
optInh.collapseTarg     = true;
optInh.printPlot        = false;
optInh.plotFlag         = false;
dataInh                 = ccm_inhibition(subjectID, sessionID, optInh);
% Which ssrt estimate should we use for latency matching? For now, use a
% collapsed, across color coherence value. Later might want to use
% estimates within each color coherence.
ssrt = dataInh.ssrtCollapseIntegrationWeighted;









% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              LOOP THROUGH UNITS/CHANNELS RECORDED
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for kDataInd = 1 : nUnit
    switch Opt.dataType
        case 'neuron'
        case 'lfp'
                            [a, kUnit] = ismember(chNum(kDataInd), SessionData.lfpChannel);
        case 'erp'
                            [a, kUnit] = ismember(Opt.unitArray{kDataInd}, eeg_electrode_map(subjectID));
    end
    
    
    % Get default trial selection Opt
    selectOpt = ccm_options;
    
    % Loop through all right targets (or collapse them if desired) and
    % account for all target pairs if the session had more than one target
    % pair
    for jTarg = 1 : nTargPair
        
        
        
        
        yMax = zeros(nEpoch, nSignal, nSSD+1, nOutcome);  % Keep track of maximum signal values, for setting y-axis limits in plots
        yMin = zeros(nEpoch, nSignal, nSSD+1, nOutcome);  % Keep track of maximum signal values, for setting y-axis limits in plots
        
        Data(kDataInd, jTarg).subjectID = subjectID;
        Data(kDataInd, jTarg).sessionID = sessionID;
        Data(kDataInd, jTarg).name = Opt.unitArray{kDataInd};
        Data(kDataInd, jTarg).ssdArray = ssdArray;
        Data(kDataInd, jTarg).pSignalArray = pSignalArray;
        
        
        
        
        
        
        
        
        
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                    LOOP THROUGH COLOR COHERENCE VALUES
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for iColor = 1 : nSignal;
            
            
            % If we're collapsing over signal strength or we actually only have
            % 2 levels of signal, determine which iPct (signal) to use this
            % loop iteration
            if iColor == 1 && nSignal == 2
                iPct = pSignalArray(pSignalArray < .5).* 100;
            elseif iColor == 2 && nSignal == 2
                iPct = pSignalArray(pSignalArray > .5).* 100;
            else
                iPct = pSignalArray(iColor).* 100;
            end
            selectOpt.rightCheckerPct = iPct;
            
            
            
            
            
            % If collapsing into all left and all right or all up/all down,
            % need to note here that there are "2" angles to deal with
            % (important for calling ccm_trial_selection.m)
            if Opt.collapseTarg && iPct(1) > 50
                jAngle = 'right';
            elseif Opt.collapseTarg && iPct(1) < 50
                jAngle = 'left';
            else
                if iPct(1) > 50
                    rightTargArray = targAngleArray(rightTargInd);
                    jAngle = rightTargArray(jTarg);
                elseif iPct(1) < 50
                    leftTargArray = targAngleArray(leftTargInd);
                    jAngle = leftTargArray(jTarg);
                end
            end
            selectOpt.targDir = jAngle;
            
            
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                    LOOP THROUGH EPOCHS
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for mEpoch = 1 : length(epochArray)
                mEpochName = epochArray{mEpoch};
                
                
                
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %              GO TRIALS
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                selectOpt.ssd = 'none';
                for g = 1 : length(goOutcomeArray)
                    % Select appropriate trials to analyzie
                    switch goOutcomeArray{g}
                        case 'goTarg'
                            selectOpt.outcome     = {'goCorrectTarget', 'targetHoldAbort', 'toneAbort'};
                            if iPct < 50
                                selectOpt.responseDir = 'left';
                            elseif iPct > 50
                                selectOpt.responseDir = 'right';
                            end
                        case 'goDist'
                            selectOpt.outcome     = {'goCorrectDistractor', 'distractorHoldAbort', 'toneAbort'};
                            if iPct < 50
                                selectOpt.responseDir = 'right';
                            elseif iPct > 50
                                selectOpt.responseDir = 'left';
                            end
                    end
                    iGoTrial = ccm_trial_selection(trialData, selectOpt);
                    
                    % Save go Target trial list for this
                    % color/epoch/outcome, to be used for latency matching
                    % below
                    if mEpoch == 1 && strcmp(goOutcomeArray{g}, 'goTarg')
                        iGoTargTrial = iGoTrial;
                    end
                    
                    if ~strcmp(mEpochName, 'stopSignalOn')  % No stop signals on go trials
                        
                        alignListGo = trialData.(mEpochName)(iGoTrial);
                        Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).alignTimeList = alignListGo;   % Keep track of trial-by-trial alignemnt
                        
                        
                        
                        switch Opt.dataType
                            
                            
                            case 'neuron'
                                % Go to TargetDistractor trials
                                [alignedRasters, alignmentIndex] = spike_to_raster(trialData.(Opt.unitArray{kDataInd})(iGoTrial), alignListGo);
                                %                                     if kUnit > 21
                                %                                         alignmentIndex = alignmentIndex - 98;
                                %                                     end
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).alignTime = alignmentIndex;
                                sdf = spike_density_function(alignedRasters, Kernel);
                                if ~isempty(sdf); yMax(mEpoch, iColor, 1, 1) = nanmax(nanmean(sdf, 1)); end;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).raster = alignedRasters;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).sdf = sdf;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).sdfMean = nanmean(sdf, 1);
                                
                                %                                     fprintf('epoch: %s\talign: %d\n', mEpochName, alignmentIndex)
                                clear alignedRasters sdf
                                
                            case 'lfp'
                                % Go  trials
                                [goLFP, alignmentIndex] 	= align_signals(trialData.lfpData(iGoTrial, kUnit), alignListGo, cropWindow);
                                satTrial                = signal_reject_saturate(goLFP, 'alignIndex', alignmentIndex);
                                iGoTrial(satTrial) = [];
                                goLFP(satTrial,:)     = [];
                                if Opt.baselineCorrect
                                    goLFP                 = signal_baseline_correct(goLFP, baseWindow, alignmentIndex);
                                end
                                if Opt.filterData
                                    goLFPMean = lowpass(nanmean(goLFP, 1)', Opt.stopHz);
                                else
                                    goLFPMean = nanmean(goLFP, 1);
                                end
                                
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).alignTime = alignmentIndex;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).lfp = goLFP;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).lfpMean = goLFPMean;
                                if ~isempty(goLFP)
                                    yMax(mEpoch, iColor, 1, 1) = nanmax(nanmean(goLFP, 1));
                                    yMin(mEpoch, iColor, 1, 1) = nanmin(nanmean(goLFP, 1));
                                end
                                clear goLFP
                                
                                
                                
                                
                                
                                %                         case 'erp'
                                %                             % Go to Target trials
                                %                             [targEEG, alignmentIndex] 	= align_signals(trialData.eegData(iGoTrial, kUnit), alignListGo, cropWindow);
                                %                             satTrial                = signal_reject_saturate(targEEG, 'alignmentIndex', alignmentIndex);
                                %                             iGoTrial(satTrial) = [];
                                %                             targEEG(satTrial,:)     = [];
                                %                             targEEG                 = signal_baseline_correct(targEEG, baseWindow, alignmentIndex);
                                %                             if Opt.filterData
                                %                                 targEEGMean = lowpass(nanmean(targEEG, 1)', Opt.stopHz);
                                %                             else
                                %                                 targEEGMean = nanmean(targEEG, 1);
                                %                             end
                                %
                                %                             Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).alignTime = alignmentIndex;
                                %                             Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).eeg = targEEG;
                                %                             Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).erp = targEEGMean;
                                %                             if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                                %                                 if ~isempty(targEEG)
                                %                                     yMax(kDataInd, mEpoch, iColor, 1) = nanmax(nanmean(targEEG, 1));
                                %                                     yMin(kDataInd, mEpoch, iColor, 1) = nanmin(nanmean(targEEG, 1));
                                %                                 end
                                %                             end
                                %
                        end % switch Opt.dataType
                    end % ~strcmp(mEpochName, 'stopSignalOn')
                    
                    
                    % add go RTs here
                    if strcmp(mEpochName, 'checkerOn')
                        Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).rt = trialData.rt(iGoTrial);
                    end
                    
                end % GO TRIAL OUTCOMES
                
                
                
                
                
                
                
                
                
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % STOP TRIALS   and lantency matched GO trials (to target
                % only - no latency matched choice error data collected
                % here
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                if Opt.doStops
                    for jSSDIndex = 1 : nSSD
                        jSSD = ssdArray(jSSDIndex);
                        selectOpt.ssd       = jSSD;
                        
                        for s = 1 : length(stopOutcomeArray)
                            % Select appropriate trials to analyzie
                            switch stopOutcomeArray{s}
                                case 'stopTarg'
                                    selectOpt.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget', 'toneAbort'};
                                    if iPct < 50
                                        selectOpt.responseDir = 'left';
                                    elseif iPct > 50
                                        selectOpt.responseDir = 'right';
                                    end
                                case 'stopDist'
                                    selectOpt.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor', 'toneAbort'};
                                    if iPct < 50
                                        selectOpt.responseDir = 'right';
                                    elseif iPct > 50
                                        selectOpt.responseDir = 'left';
                                    end
                                case 'stopStop'
                                    selectOpt.outcome       = {'stopCorrect', 'toneAbort'};
                                    selectOpt.responseDir = 'none';
                            end
                            jStopTrial = ccm_trial_selection(trialData, selectOpt);
                            
                            
                            alignListStop = trialData.(mEpochName)(jStopTrial);
                            
                            % Save the list of alignement for stopStop
                            % checkerOn, to use below for aligning stopStop
                            % "responseOnset" to checker Onset plut mean Go
                            % Slow RT
                            if strcmp(stopOutcomeArray{s}, 'stopStop') && strcmp(mEpochName, 'checkerOn')
                                stopStopCheckerTrial = jStopTrial;
                                stopStopCheckerOnList = alignListStop;
                            end
                            
                            Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).alignTimeList = alignListStop;   % Keep track of trial-by-trial alignemnt
                            
                            
                            
                            switch Opt.dataType
                                case 'neuron'
                                    
                                    
                                    if ~(strcmp(stopOutcomeArray{s}, 'stopStop') && strcmp(mEpochName, 'responseOnset'))  % No stop signals on go trials
                                        %  Stop trials
                                        [alignedRasters, alignmentIndex] = spike_to_raster(trialData.(Opt.unitArray{kDataInd})(jStopTrial), alignListStop);
                                        Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).alignTime = alignmentIndex;
                                        
                                        sdf = spike_density_function(alignedRasters, Kernel);
                                        if ~isempty(sdf); yMax(mEpoch, iColor, jSSDIndex+1, 5) = nanmax(nanmean(sdf, 1)); end;
                                        Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).raster = alignedRasters;
                                        Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).sdf = sdf;
                                        Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).sdfMean = nanmean(sdf, 1);
                                        
                                    end
                                    clear alignedRasters sdf
                                    
                                    
                                    
                                    
                                case 'lfp'
                                    
                                    
                                    if ~(strcmp(stopOutcomeArray{s}, 'stopStop') && strcmp(mEpochName, 'responseOnset'))  % No stop signals on go trials
                                        %  Stop trials
                                        [stopLFP, alignmentIndex] 	= align_signals(trialData.lfpData(jStopTrial, kUnit), alignListStop, cropWindow);
                                        satTrial                = signal_reject_saturate(stopLFP, 'alignmentIndex', alignmentIndex);
                                        stopLFP(satTrial,:)     = [];
                                        if Opt.baselineCorrect
                                            stopLFP                 = signal_baseline_correct(stopLFP, baseWindow, alignmentIndex);
                                        end
                                        if Opt.filterData
                                            stopLFPMean = lowpass(nanmean(stopLFP, 1)', Opt.stopHz);
                                        else
                                            stopLFPMean = nanmean(stopLFP, 1);
                                        end
                                        
                                        Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).alignTime = alignmentIndex;
                                        Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).lfp = stopLFP;
                                        Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).lfpMean = stopLFPMean;
                                        if ~isempty(stopLFP)
                                            yMax(kDataInd, mEpoch, iColor, jSSDIndex + 1, 5) = nanmax(nanmean(stopLFP, 1));
                                            yMin(kDataInd, mEpoch, iColor, jSSDIndex + 1, 5) = nanmin(nanmean(stopLFP, 1));
                                        end
                                        
                                    end
                                    clear stopLFP
                                    
                                    %
                                    %                             case 'erp'
                                    %                                 if ~strcmp(mEpochName, 'responseOnset')  % No response on correct Stop trials
                                    %                                     % Correct Stop trials
                                    %                                     [stopEEG, alignmentIndex] 	= align_signals(trialData.eegData(jStopTrial, kUnit), alignListstopStop, cropWindow);
                                    %                                     satTrial                = signal_reject_saturate(stopEEG, 'alignmentIndex', alignmentIndex);
                                    %                                     stopEEG(satTrial,:)     = [];
                                    %                                     stopEEG                 = signal_baseline_correct(stopEEG, baseWindow, alignmentIndex);
                                    %                                     if Opt.filterData
                                    %                                         stopEEGMean = lowpass(nanmean(stopEEG, 1)', Opt.stopHz);
                                    %                                     else
                                    %                                         stopEEGMean = nanmean(stopEEG, 1);
                                    %                                     end
                                    %
                                    %                                     Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).alignTime = alignmentIndex;
                                    %                                     Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).eeg = stopEEG;
                                    %                                     Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).erp = stopEEGMean;
                                    %                                     if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                                    %                                         if ~isempty(stopEEG)
                                    %                                             yMax(kDataInd, mEpoch, iColor, jSSDIndex + 1, 5) = nanmax(nanmean(stopEEG, 1));
                                    %                                             yMin(kDataInd, mEpoch, iColor, jSSDIndex + 1, 5) = nanmin(nanmean(stopEEG, 1));
                                    %                                         end
                                    %                                     end
                                    %                                 end
                                    
                            end % switch Opt.dataType
                            
                            
                            % Move stop rts here
                            if strcmp(mEpochName, 'checkerOn') && ~strcmp(stopOutcomeArray(s), 'stopStop')
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).rt = trialData.rt(jStopTrial);
                            end
                            
                        end % stopOutcomeArray
                        
                        
                        
                        
                        % Get latency matched Go trials
                        mOpt.matchMethod = 'ssrt';
                        mOpt.ssrt = ssrt;
                        mOpt.ssd = ssdArray(jSSDIndex);
                        mOpt.stopRT = trialData.rt(jStopTrial);
                        [goFastTrial, goSlowTrial] = latency_match(trialData.rt(iGoTargTrial), mOpt);
                        goFastTrial = iGoTargTrial(goFastTrial);
                        goSlowTrial = iGoTargTrial(goSlowTrial);
                        
                        % Go fast and slow rts
                        if strcmp(mEpochName, 'checkerOn')
                            Data(kDataInd, jTarg).checkerOn.colorCoh(iColor).goFast.ssd(jSSDIndex).rt = trialData.rt(goFastTrial);
                            Data(kDataInd, jTarg).checkerOn.colorCoh(iColor).goSlow.ssd(jSSDIndex).rt = trialData.rt(goSlowTrial);
                            
                            % Use mean RTs during GO Slow trials to define
                            % a "responseOnset" alignement for stopStop
                            % trials
                            alignListStop = stopStopCheckerOnList + round(nanmean(trialData.rt(goSlowTrial)));
                            [alignedRasters, alignmentIndex] = spike_to_raster(trialData.(Opt.unitArray{kDataInd})(stopStopCheckerTrial), alignListStop);
                            Data(kDataInd, jTarg).responseOnset.colorCoh(iColor).stopStop.ssd(jSSDIndex).alignTime = alignmentIndex;
                            
                            sdf = spike_density_function(alignedRasters, Kernel);
                            Data(kDataInd, jTarg).responseOnset.colorCoh(iColor).stopStop.ssd(jSSDIndex).raster = alignedRasters;
                            Data(kDataInd, jTarg).responseOnset.colorCoh(iColor).stopStop.ssd(jSSDIndex).sdf = sdf;
                            Data(kDataInd, jTarg).responseOnset.colorCoh(iColor).stopStop.ssd(jSSDIndex).sdfMean = nanmean(sdf, 1);
                            clear sdf alignedRasters
                            
                        end
                        
                        alignListGoFast = trialData.(mEpochName)(goFastTrial);
                        Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goFast.ssd(jSSDIndex).alignTimeList = alignListGoFast;   % Keep track of trial-by-trial alignemnt
                        alignListGoSlow = trialData.(mEpochName)(goSlowTrial);
                        Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goSlow.ssd(jSSDIndex).alignTimeList = alignListGoSlow;   % Keep track of trial-by-trial alignemnt
                        
                        switch Opt.dataType
                            case 'neuron'
                                % Go Fast data
                                [alignedRasters, alignmentIndex] = spike_to_raster(trialData.(Opt.unitArray{kDataInd})(goFastTrial), alignListGoFast);
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goFast.ssd(jSSDIndex).alignTime = alignmentIndex;
                                
                                sdf = spike_density_function(alignedRasters, Kernel);
                                if ~isempty(sdf); yMax(mEpoch, iColor, jSSDIndex+1, 3) = nanmax(nanmean(sdf, 1)); end;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goFast.ssd(jSSDIndex).raster = alignedRasters;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goFast.ssd(jSSDIndex).sdf = sdf;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goFast.ssd(jSSDIndex).sdfMean = nanmean(sdf, 1);
                                clear sdf alignedRasters
                                
                                
                                % Go Slow data
                                [alignedRasters, alignmentIndex] = spike_to_raster(trialData.(Opt.unitArray{kDataInd})(goSlowTrial), alignListGoSlow);
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goSlow.ssd(jSSDIndex).alignTime = alignmentIndex;
                                
                                sdf = spike_density_function(alignedRasters, Kernel);
                                if ~isempty(sdf); yMax(mEpoch, iColor, jSSDIndex+1, 3) = nanmax(nanmean(sdf, 1)); end;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goSlow.ssd(jSSDIndex).raster = alignedRasters;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goSlow.ssd(jSSDIndex).sdf = sdf;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goSlow.ssd(jSSDIndex).sdfMean = nanmean(sdf, 1);
                                clear sdf alignedRasters
                                
                            case 'lfp'
                                % Go Fast data
                                [goLFP, alignmentIndex] 	= align_signals(trialData.lfpData(goFastTrial, kUnit), alignListGoFast, cropWindow);
                                satTrial                = signal_reject_saturate(goLFP, 'alignIndex', alignmentIndex);
                                iGoTrial(satTrial) = [];
                                goLFP(satTrial,:)     = [];
                                if Opt.baselineCorrect
                                    goLFP                 = signal_baseline_correct(goLFP, baseWindow, alignmentIndex);
                                end
                                if Opt.filterData
                                    goLFPMean = lowpass(nanmean(goLFP, 1)', Opt.stopHz);
                                else
                                    goLFPMean = nanmean(goLFP, 1);
                                end
                                
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goFast.ssd(jSSDIndex).alignTime = alignmentIndex;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goFast.ssd(jSSDIndex).lfp = goLFP;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goFast.ssd(jSSDIndex).lfpMean = goLFPMean;
                                clear goLFP
                                % Go Slow data
                                [goLFP, alignmentIndex] 	= align_signals(trialData.lfpData(goSlowTrial, kUnit), alignListGoSlow, cropWindow);
                                satTrial                = signal_reject_saturate(goLFP, 'alignIndex', alignmentIndex);
                                iGoTrial(satTrial) = [];
                                goLFP(satTrial,:)     = [];
                                if Opt.baselineCorrect
                                    goLFP                 = signal_baseline_correct(goLFP, baseWindow, alignmentIndex);
                                end
                                if Opt.filterData
                                    goLFPMean = lowpass(nanmean(goLFP, 1)', Opt.stopHz);
                                else
                                    goLFPMean = nanmean(goLFP, 1);
                                end
                                
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goSlow.ssd(jSSDIndex).alignTime = alignmentIndex;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goSlow.ssd(jSSDIndex).lfp = goLFP;
                                Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).goSlow.ssd(jSSDIndex).lfpMean = goLFPMean;
                                
                                clear goLFP
                            case 'erp'
                                disp('Need to implement latency matching in erp data')
                        end
                        
                    end % jSSD
                end % if Opt.doStops
                
                
                
                
            end % Epochs loop
        end %iColor
        
        
        
        
        
        
        
        
        a = yMax(:,:,1,:);   % the go trials
        b = yMax(:,:,2:end,:);  % all the stop trials
        c = yMin(:,:,1,:);   % go trials
        d = yMin(:,:,2:end,:);  % stop trials
        
        Data(kDataInd, jTarg).yMax = max([a(:); nanmean(b(:))]);
        Data(kDataInd, jTarg).yMin = min([c(:); nanmean(d(:))]);
        
        
        
        
        
        
        % print the figure if we're stepping through
        if Opt.plotFlag
            howProcess = Opt.howProcess;
            switch howProcess
                case {'step','print','each'}
                    ccm_session_data_plot(Data(kDataInd, jTarg), Opt)
                    if strcmp(howProcess, 'step')
                        pause
                    end
                    if nUnit > 1
                        clear Data
                    end
            end
        end
    end % jTarg
end % kDataInd

% if Opt.plotFlag
%     switch Opt.unitArray
%         case {'step','print'}
%             close(gcf)
%     end
% end








% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                   Opt.normalize
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% % If we want the signal to be normalized to the maximum value for each
% % unit/channel, loop back through and Opt.normalize here...
% %
% %  Haven't tested this in a while, probably needs work
%
% if Opt.normalize
%     for kDataInd = 1 : nUnit
%         for iColor = 1 : nSignal;
%             for mEpoch = 1 : length(epochArray)
%                 mEpochName = epochArray{mEpoch};
%
%                 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % Go trials
%                 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 if ~strcmp(mEpochName, 'stopSignalOn')  % No stop signals on go trials
%
%                     switch Opt.dataType
%                         case 'neuron'
%                             % Go to Target trials
%                             Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).sdf = ...
%                                 Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).sdf ./ Data(kDataInd, jTarg).yMax;
%                             Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).sdfMean = ...
%                                 nanmean(Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).sdf, 1);
%
%
%                         case 'lfp'
%                             % Go to Target trials
%                             Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).lfp = ...
%                                 Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).lfp ./ Data(kDataInd, jTarg).yMax;
%                             Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).lfpMean = ...
%                                 nanmean(Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).lfp, 1);
%
%
%                         case 'erp'
%                             % Go to Target trials
%                             Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).eeg = ...
%                                 Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).eeg ./ Data(kDataInd, jTarg).yMax;
%                             Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).lfpMean = ...
%                                 nanmean(Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(goOutcomeArray{g}).eeg, 1);
%
%                     end % switch Opt.dataType
%                 end % ~strcmp(mEpochName, 'stopSignalOn')
%
%
%
%
%
%                 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % Stop trials
%                 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 if Opt.doStops
%                     for jSSDIndex = 1 : nSSD
%
%                         switch Opt.dataType
%                             case 'neuron'
%                                 if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
%
%                                     % Stop Correct trials
%                                     Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).sdf = ...
%                                         Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).sdf ./ Data(kDataInd, jTarg).yMax;
%                                     Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).sdfMean = ...
%                                         nanmean(Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).sdf, 1);
%
%                                 end % ~strcmp(mEpochName, 'responseOnset')
%
%
%                             case 'lfp'
%
%                                     % Stop Correct trials
%                                     Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).lfp = ...
%                                         Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).lfp ./ Data(kDataInd, jTarg).yMax;
%                                     Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).lfpMean = ...
%                                         nanmean(Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).lfp, 1);
%
%                                 end % ~strcmp(mEpochName, 'responseOnset')
%
%
%                             case 'erp'
%
%                                 if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
%
%                                     % Stop Correct trials
%                                     Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).eeg = ...
%                                         Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).eeg ./ Data(kDataInd, jTarg).yMax;
%                                     Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).erp = ...
%                                         nanmean(Data(kDataInd, jTarg).(mEpochName).colorCoh(iColor).(stopOutcomeArray{s}).ssd(jSSDIndex).eeg, 1);
%
%                                 end % ~strcmp(mEpochName, 'responseOnset')
%
%                         end % switch Opt.dataType
%                     end % jSSD
%                 end % if Opt.doStops
%             end % mEpoch
%         end % for mEpoch = 1 : length(epochArray)
%         Data(kDataInd, jTarg).yMax = 1.1;
%         switch Opt.dataType
%             case 'neuron'
%                 Data(kDataInd, jTarg).yMin = 0;
%             case {'lfp', 'erp'}
%                 Data(kDataInd, jTarg).yMin = -1.1;
%         end
%
%
%
%     end % for kDataInd = 1 : nUnit
%
%
% end % if Opt.normalize







%     disp('completed data collection')

Data(1).unitArray       = Opt.unitArray;
Data(1).howProcess       = Opt.howProcess;
Data(1).pSignalArray    = pSignalArray;
Data(1).targAngleArray  = targAngleArray;
Data(1).ssdArray        = ssdArray;
Data(1).sessionID       = sessionID;
Data(1).subjectID       = subjectID;
% Data(1).hemisphere       = SessionData.hemisphere;
Data(1).spikeUnitArray       = SessionData.spikeUnitArray;
Data(1).Opt             = Opt;


if Opt.plotFlag &&...
        ~strcmp(Opt.howProcess, 'step') &&...
        ~strcmp(Opt.howProcess, 'each') &&...
        ~strcmp(Opt.howProcess, 'pop') &&...
        ~strcmp(Opt.howProcess, 'print')
    ccm_session_data_plot(Data, Opt)
elseif strcmp(Opt.howProcess, 'pop')
    ccm_session_data_plot_pop(Data, Opt)
end

% ################################################################################
% ############################   SUB FUNCTIONS    ################################
% ################################################################################

%     function [goFastTrial, goSlowTrial] = latency_match(goRT, stopRT)
%
%         goFastTrial = [];
%         goSlowTrial = [];
%         if isempty(stopRT) || isempty(goRT)
%             return
%         end
%
%         [rt, ind] = sort(goRT);
%
%         lastInd = 1;
%         while nanmean(rt(lastInd:end)) < nanmean(stopRT)
%             goFastTrial = [goFastTrial; ind(lastInd)];
%             lastInd = lastInd + 1;
%
%
%         end
%         goSlowTrial = ind(lastInd:end);
%         %         trialList = 1 : length(goRT);
%         %         remaining = 1 : length(goRT);
%         %         remove = [];
%         %         goFastTrial = remaining;
%         %         while nanmean(goRT(remaining)) > nanmean(stopRT)
%         %             [y,i] = max(goRT(remaining));
%         %             remove = [remove; i]
%         %             goSlowTrial = [goSlowTrial; i];
%         %
%         %
%         %         end
%         %         goFastTrial = remaining;
%
%
%     end

end

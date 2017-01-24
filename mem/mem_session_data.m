function [Data, options] = mem_session_data(subjectID, sessionID, options)

%
% function Unit = mem_session_data(subjectID, sessionID, dataType, varargin)
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
%    options.collapseTarg         = false, true;
%    options.filterData 	= false, true;
%    options.normalize      = false, true;
%    options.unitArray      = {'spikeUnit17a'},'each', units want to analyze
%
%
% Returns Unit structure with fields:
%
%   Unit.signalStrength(x).(condition).ssd(x).(epoch name)
%
%   condition can be:  rightTarg, goDist, stopTarg, stopDist, stopCorrect
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
clear Data

if nargin < 3
    options.dataType = 'neuron';
    
    options.figureHandle     = 1000;
    options.printPlot        = true;
    options.plotFlag         = true;
    options.collapseSignal   = false;
    options.collapseTarg      = false;
    options.doStops          = true;
    options.filterData       = false;
    options.stopHz           = 50;
    options.normalize        = false;
    options.unitArray        = 'print';
    options.baselineCorrect  = true;
    
    if nargin == 0
        Data = options;
        return
    end
end
collapseSignal  = options.collapseSignal;
doStops         = options.doStops;
normalize       = options.normalize;
filterData      = options.filterData;
unitArray       = options.unitArray;
baselineCorrect = options.baselineCorrect;



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

if ~(strcmp(SessionData.taskID, 'mem') || strcmp(SessionData.taskID, 'del'))
    fprintf('Not a memory guided saccade session, try again\n')
    return
end



nTarg = length(unique(trialData.targAngle));
nTrial = size(trialData.trialOutcome, 1);
% if nTarg ~= 2
%     disp('Don"t have code yet for memory guided saccade task with more than 2 targets')
%     Unit = struct();
%     return
% end


% CONSTANTS
Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;

tickWidth = 5;
rasYlim = 100;

cropWindow  = -499 : 500;  % used to extract a semi-small portion of signal for each epoch/alignemnt
baseWindow 	= -149 : 0;   % To baseline-shift the eeg signals, relative to event alignment index;
stopHz      = 50;

baseEpoch = -99 : 0;
visEpoch = mem_epoch_range('targOn', 'analyze');
movEpoch = mem_epoch_range('responseOnset', 'analyze');


% Set defaults
dataType = options.dataType;
switch dataType
    case 'neuron'
        dataArray     = SessionData.spikeUnitArray;
    case 'lfp'
        chNum = SessionData.lfpChannel;
        dataArray 	= num2cell(SessionData.lfpChannel);
        dataArray   = cellfun(@(x) sprintf('lfp_%s', num2str(x, '%02d')), dataArray, 'uniformoutput', false);
    case 'erp'
        dataArray     = eeg_electrode_map(subjectID);
end
if strcmp(unitArray, 'each') || strcmp(unitArray, 'step') || strcmp(unitArray, 'print')
    unitArray     = dataArray;
end

% Make sure user input a dataType that was recorded during the session
dataTypePossible = {'neuron', 'lfp', 'erp'};
if ~sum(strcmp(dataType, dataTypePossible))
    fprintf('%s Is not a valid data type \n', dataType)
    return
end
if isempty(unitArray)
    fprintf('Session %s apparently does not contain %s data \n', sessionID, dataType)
    return
end


% Make sure user input a dataType that was recorded during the session
dataTypePossible = {'neuron', 'lfp', 'erp'};
if ~sum(strcmp(dataType, dataTypePossible))
    fprintf('%s Is not a valid data type \n', dataType)
    return
end










% epochArray = {'fixWindowEntered', 'targOn', 'checkerOn', 'stopSignalOn', 'responseOnset', 'rewardOn'};
epochArray = {'fixWindowEntered', 'targOn', 'responseCueOn', 'responseOnset', 'rewardOn'};
nEpoch = length(epochArray);


% How many units were recorded?
% nUnit = size(dataArray, 2);
nUnit = length(unitArray);

% nUnit = 1;










[minTrialDuration,i] = min(trialData.trialDuration(~strcmp(trialData.trialOutcome, 'fixationAbort')));
[maxTrialDuration,i] = max(trialData.trialDuration(~strcmp(trialData.trialOutcome, 'fixationAbort')));


yMax = zeros(nUnit, nEpoch, nTarg);  % Keep track of maximum sdf values, for setting y-axis limits in plots
yMin = zeros(nUnit, nEpoch, nTarg);  % Keep track of maximum sdf values, for setting y-axis limits in plots

for kDataIndex = 1 : nUnit
    switch dataType
        case 'neuron'
            [a, kUnit] = ismember(unitArray{kDataIndex}, SessionData.spikeUnitArray);
        case 'lfp'
            [a, kUnit] = ismember(unitArray{kDataIndex}, SessionData.lfpChannel);
        case 'erp'
            [a, kUnit] = ismember(unitArray{kDataIndex}, eeg_electrode_map(subjectID));
    end
    
    Data(kDataIndex).subjectID = subjectID;
    Data(kDataIndex).sessionID = sessionID;
    Data(kDataIndex).name = unitArray{kDataIndex};
    
    rightTargTrial = mem_trial_selection(trialData, {'saccToTarget'}, 'right');
    leftTargTrial = mem_trial_selection(trialData, {'saccToTarget'}, 'left');
    
    for mEpoch = 1 : length(epochArray)
        mEpochName = epochArray{mEpoch};
        
        
        alignListR = trialData.(mEpochName)(rightTargTrial);
        Data(kDataIndex).rightTarg.(mEpochName).alignTimeList = alignListR;   % Keep track of trial-by-trial alignemnt
        alignListL = trialData.(mEpochName)(leftTargTrial);
        Data(kDataIndex).leftTarg.(mEpochName).alignTimeList = alignListL;   % Keep track of trial-by-trial alignemnt
        
        
        
        switch dataType
            case 'neuron'
                % Right Target trials
                [alignedRasters, alignIndex] = spike_to_raster(trialData.spikeData(rightTargTrial, kUnit), alignListR);
                % add NaN pad because often we want to display the sdf longer, and low-firing cells will otherwise get cut off
                alignedRasters = [alignedRasters, nan(length(alignListR), 6000)];
                Data(kDataIndex).rightTarg.(mEpochName).alignTime = alignIndex;
                sdf = spike_density_function(alignedRasters, Kernel);
                if ~isempty(sdf); yMax(kDataIndex, mEpoch, 1) = nanmax(nanmean(sdf, 1)); end;
                Data(kDataIndex).rightTarg.(mEpochName).raster = alignedRasters;
                Data(kDataIndex).rightTarg.(mEpochName).signalFn = sdf;
                Data(kDataIndex).rightTarg.(mEpochName).signalMean = nanmean(sdf, 1);
                
                % Left Target trials
                [alignedRasters, alignIndex] = spike_to_raster(trialData.spikeData(leftTargTrial, kUnit), alignListL);
                alignedRasters = [alignedRasters, nan(length(alignListL), 6000)];
                Data(kDataIndex).leftTarg.(mEpochName).alignTime = alignIndex;
                sdf = spike_density_function(alignedRasters, Kernel);
                if ~isempty(sdf); yMax(kDataIndex, mEpoch, 2) = nanmax(nanmean(sdf, 1)); end;
                Data(kDataIndex).leftTarg.(mEpochName).raster = alignedRasters;
                Data(kDataIndex).leftTarg.(mEpochName).signalFn = sdf;
                Data(kDataIndex).leftTarg.(mEpochName).signalMean = nanmean(sdf, 1);
                
                
                
            case 'lfp'
                % Right Target trials
                [targLFP, alignIndex] 	= align_signals(trialData.lfpData(rightTargTrial, kUnit), alignListR, cropWindow);
                satTrial                = signal_reject_saturate(targLFP, 'alignIndex', alignIndex);
                targLFP(satTrial,:)     = [];
                targLFP                 = signal_baseline_correct(targLFP, baseWindow, alignIndex);
                if filterData
                    targLFPMean = lowpass(nanmean(targLFP, 1)', stopHz);
                else
                    targLFPMean = nanmean(targLFP, 1);
                end
                Data(kDataIndex).rightTarg.(mEpochName).alignTime = alignIndex;
                
                Data(kDataIndex).rightTarg.(mEpochName).signalFn = targLFP;
                Data(kDataIndex).rightTarg.(mEpochName).signalMean = targLFPMean;
                if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                    if ~isempty(targLFP);
                        yMax(kDataIndex, mEpoch, 1) = nanmax(nanmean(targLFP, 1));
                        yMin(kDataIndex, mEpoch, 1) = nanmin(nanmean(targLFP, 1));
                    end
                end
                
                % Left Target trials
                [targLFP, alignIndex] 	= align_signals(trialData.lfpData(leftTargTrial, kUnit), alignListR, cropWindow);
                satTrial                = signal_reject_saturate(targLFP, 'alignIndex', alignIndex);
                targLFP(satTrial,:)     = [];
                targLFP                 = signal_baseline_correct(targLFP, baseWindow, alignIndex);
                if filterData
                    targLFPMean = lowpass(nanmean(targLFP, 1)', stopHz);
                else
                    targLFPMean = nanmean(targLFP, 1);
                end
                Data(kDataIndex).leftTarg.(mEpochName).alignTime = alignIndex;
                
                Data(kDataIndex).leftTarg.(mEpochName).signalFn = targLFP;
                Data(kDataIndex).leftTarg.(mEpochName).signalMean = targLFPMean;
                if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                    if ~isempty(targLFP);
                        yMax(kDataIndex, mEpoch, 1) = nanmax(nanmean(targLFP, 1));
                        yMin(kDataIndex, mEpoch, 1) = nanmin(nanmean(targLFP, 1));
                    end
                end
                
                
                
                
            case 'erp'
                % Right Target trials
                [targEEG, alignIndex] 	= align_signals(trialData.eegData(rightTargTrial, kUnit), alignListR, cropWindow);
                satTrial                = signal_reject_saturate(targEEG, 'alignIndex', alignIndex);
                targEEG(satTrial,:)     = [];
                targEEG                 = signal_baseline_correct(targEEG, baseWindow, alignIndex);
                if filterData
                    targERP = lowpass(nanmean(targEEG, 1)', stopHz);
                else
                    targERP = nanmean(targEEG, 1);
                end
                Data(kDataIndex).rightTarg.(mEpochName).alignTime = alignIndex;
                
                Data(kDataIndex).rightTarg.(mEpochName).signalFn = targEEG;
                Data(kDataIndex).rightTarg.(mEpochName).signalMean = targERP;
                if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                    if ~isempty(targEEG);
                        yMax(kDataIndex, mEpoch, 1) = nanmax(nanmean(targEEG, 1));
                        yMin(kDataIndex, mEpoch, 1) = nanmin(nanmean(targEEG, 1));
                    end
                end
                
                
                % Left Target trials
                [targEEG, alignIndex] 	= align_signals(trialData.eegData(leftTargTrial, kUnit), alignListR, cropWindow);
                satTrial                = signal_reject_saturate(targEEG, 'alignIndex', alignIndex);
                targEEG(satTrial,:)     = [];
                targEEG                 = signal_baseline_correct(targEEG, baseWindow, alignIndex);
                if filterData
                    targERP = lowpass(nanmean(targEEG, 1)', stopHz);
                else
                    targERP = nanmean(targEEG, 1);
                end
                Data(kDataIndex).leftTarg.(mEpochName).alignTime = alignIndex;
                
                Data(kDataIndex).leftTarg.(mEpochName).signalFn = targEEG;
                Data(kDataIndex).leftTarg.(mEpochName).signalMean = targERP;
                if ~strcmp(mEpochName, 'rewardOn') && ~strcmp(mEpochName, 'responseOnset')
                    if ~isempty(targEEG);
                        yMax(kDataIndex, mEpoch, 1) = nanmax(nanmean(targEEG, 1));
                        yMin(kDataIndex, mEpoch, 1) = nanmin(nanmean(targEEG, 1));
                    end
                end
        end
        
        
    end % mEpoch
    
    yLimMax = max(yMax(kDataIndex, :)) * 1.1;
    %          yLimMax = 55;
    switch dataType
        case 'neuron'
            yLimMin = 0;
        case {'lfp', 'erp'}
            yLimMin = min(yMin(kDataIndex, :)) * 1.1;
    end
    Data(kDataIndex).yMax = yLimMax;
    Data(kDataIndex).yMin = yLimMin;
    
    
    
    
    if strcmp(dataType, 'neuron')
        classify_cell_type
    end
    
    % print the figure if we're stepping through
    if options.plotFlag
        switch options.unitArray
            case {'step','print'}
                mem_session_data_plot(Data(kDataIndex), options)
                if strcmp(options.unitArray, 'step')
                    pause
                end
                clear Data
        end
    end
end % kUnitIndex





if options.plotFlag && ~strcmp(options.unitArray, 'step') && ~strcmp(options.unitArray, 'each') && ~strcmp(options.unitArray, 'print')
    mem_session_data_plot(Data, options)
end



Data(1).unitArray       = unitArray;
Data(1).dataArray       = dataArray;
Data(1).sessionID       = sessionID;
Data(1).subjectID       = subjectID;
Data(1).options         = options;








%**********************   SUBFUNCTIONS    %*************************
%*******************************************************************


%**********************   CLASSIFY CELL RESPONSE TYPE    %*************************
    function classify_cell_type
        if strcmp(dataType, 'neuron')
            %     for kDataIndex = 1 : nUnit
            
            
            vis = false;
            mov = false;
            alpha = .05;
            baseEpoch = -99 : 0;
            visEpoch = mem_epoch_range('targOn', 'analyze');
            movEpoch = mem_epoch_range('responseOnset', 'analyze');
            
            
            baseRasR = Data(kDataIndex).rightTarg.targOn.raster(:, baseEpoch + Data(kDataIndex).rightTarg.targOn.alignTime);
            baseRasL = Data(kDataIndex).leftTarg.targOn.raster(:, baseEpoch + Data(kDataIndex).leftTarg.targOn.alignTime);
            nBaseSpikeR = nansum(baseRasR, 2);
            nBaseSpikeL = nansum(baseRasL, 2);
            
            
            % Visual cell?
            rightRas = Data(kDataIndex).rightTarg.targOn.raster(:, visEpoch + Data(kDataIndex).rightTarg.targOn.alignTime);
            nRightSpike = nansum(rightRas, 2);
            [h, p, ci, stats] = ttest2(nBaseSpikeR, nRightSpike);
            if p < alpha && nanmean(nBaseSpikeR) < nanmean(nRightSpike)
                vis = true;
            end
            
            leftRas = Data(kDataIndex).leftTarg.targOn.raster(:, visEpoch + Data(kDataIndex).leftTarg.targOn.alignTime);
            nLeftSpike = nansum(leftRas, 2);
            [h, p, ci, stats] = ttest2(nBaseSpikeL, nLeftSpike);
            if p < alpha && nanmean(nBaseSpikeL) < nanmean(nLeftSpike)
                vis = true;
            end
            
            
            % Movement cell?
            rightRas = Data(kDataIndex).rightTarg.responseOnset.raster(:, movEpoch + Data(kDataIndex).rightTarg.responseOnset.alignTime);
            nRightSpike = nansum(rightRas, 2);
            [h, p, ci, stats] = ttest2(nBaseSpikeR, nRightSpike);
            if p < alpha && nanmean(nBaseSpikeR) < nanmean(nRightSpike)
                mov = true;
            end
            
            leftRas = Data(kDataIndex).leftTarg.responseOnset.raster(:, movEpoch + Data(kDataIndex).leftTarg.responseOnset.alignTime);
            nLeftSpike = nansum(leftRas, 2);
            [h, p, ci, stats] = ttest2(nBaseSpikeL, nLeftSpike);
            if p < alpha && nanmean(nBaseSpikeL) < nanmean(nLeftSpike)
                mov = true;
            end
            
            
            
            if ~vis && ~mov
                Data(kDataIndex).cellType = nan;
            elseif vis && ~mov
                Data(kDataIndex).cellType = 'visual';
            elseif ~vis && mov
                Data(kDataIndex).cellType = 'movement';
            elseif vis && mov
                Data(kDataIndex).cellType = 'visuomovement';
            end
            %     end % kUnitIndex
        end
    end


















%**********************   PLOTTING    %*************************
    function mem_session_data_plot(PlotData, options)
        
        if options.plotFlag
            cMap = ccm_colormap([0 1]);
            kernelMethod = 'gaussian';
            SIGMA = 20;
            kernelMethod = 'postsynaptic potential';
            GROWTH = 1;
            DECAY = 20;
            % kernelIn = [SIGMA];
            
            
            [nUnitPlot, nTargPair] = size(PlotData);
            figureHandle    = options.figureHandle;
            
            
            targLineW = 2;
            
            
            for kDataPlot = 1 : nUnitPlot
                %       colormap([1 1 1; cMap(1,:); cMap(2,:)])
                nRow = 2;
                nEpoch = length(epochArray);
                nColumn = nEpoch;
                figureHandle = figureHandle + 1;
                if options.printPlot
                    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
                else
                    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
                end
                clf
                for mEpochPlot = 1 : nEpoch
                    mEpochName = epochArray{mEpochPlot};
                    epochRangeDisplay = mem_epoch_range(mEpochName, 'plot');
                    epochRangeAnalysis = mem_epoch_range(mEpochName, 'analyze');
                    
                    
                    % _______  Set up axes  ___________
                    % axes names
                    axSig = 1;
                    axRas = 2;
                    
                    % Set up plot axes
                    % SDFs
                    ax(axSig, mEpochPlot) = axes('units', 'centimeters', 'position', [xAxesPosition(axSig, mEpochPlot) yAxesPosition(axSig, mEpochPlot) axisWidth axisHeight]);
                    set(ax(axSig, mEpochPlot), 'ylim', [yLimMin yLimMax], 'xlim', [epochRangeDisplay(1) epochRangeDisplay(end)])
                    cla
                    hold(ax(axSig, mEpochPlot), 'on')
                    title(epochArray{mEpochPlot})
                    
                    
                    % Rasters
                    ax(axRas, mEpochPlot) = axes('units', 'centimeters', 'position', [xAxesPosition(axRas, mEpochPlot) yAxesPosition(axRas, mEpochPlot) axisWidth axisHeight]);
                    set(ax(axRas, mEpochPlot), 'ylim', [0 rasYlim], 'xlim', [epochRangeDisplay(1) epochRangeDisplay(end)])
                    cla
                    hold(ax(axRas, mEpochPlot), 'on')
                    
                    
                    if mEpochPlot > 1
                        set(ax(axSig, mEpochPlot), 'yticklabel', [])
                        set(ax(axRas, mEpochPlot), 'yticklabel', [])
                        
                        set(ax(axSig, mEpochPlot), 'ycolor', [1 1 1])
                        set(ax(axRas, mEpochPlot), 'ycolor', [1 1 1])
                    end
                    
                    
                    
                    alignRightTarg = PlotData(kDataPlot).rightTarg.(mEpochName).alignTime;
                    alignLeftTarg = PlotData(kDataPlot).leftTarg.(mEpochName).alignTime;
                    
                    axes(ax(axSig, mEpochPlot))
                    fillX = [epochRangeAnalysis(1), epochRangeAnalysis(end), epochRangeAnalysis(end), epochRangeAnalysis(1)];
                    fillY = [yLimMin yLimMin yLimMax yLimMax];
                    fillColor = [1 1 .5];
                    h = fill(fillX, fillY, fillColor);
                    set(h, 'edgecolor', 'none');
                    
                    if strcmp(mEpochName, 'targOn')
                        fillX = [baseEpoch(1), baseEpoch(end), baseEpoch(end), baseEpoch(1)];
                        fillY = [yLimMin yLimMin yLimMax yLimMax];
                        fillColor = [.8 .8 .8];
                        b = fill(fillX, fillY, fillColor);
                        set(b, 'edgecolor', 'none');
                    end
                    
                    % Alignemnt line
                    plot(ax(axSig, mEpochPlot), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
                    
                    
                    
                    switch dataType
                        case 'neuron'
                            % Right SDFs
                            if ~isempty(alignRightTarg)
                                sdfRightTarg = PlotData(kDataPlot).rightTarg.(mEpochName).signalMean;
                                plot(ax(axSig, mEpochPlot), epochRangeDisplay, sdfRightTarg(alignRightTarg + epochRangeDisplay), 'color', cMap(2,:), 'linewidth', targLineW)
                            end
                            % Left SDFs
                            if ~isempty(alignLeftTarg)
                                sdfLeftTarg = PlotData(kDataPlot).leftTarg.(mEpochName).signalMean;
                                plot(ax(axSig, mEpochPlot), epochRangeDisplay, sdfLeftTarg(alignLeftTarg + epochRangeDisplay), 'color', cMap(1,:), 'linewidth', targLineW)
                            end
                            
                            axes(ax(axRas, mEpochPlot))
                            colormap([1 1 1; cMap(1,:); cMap(2,:)])
                            rasRightTarg =  PlotData(kDataPlot).rightTarg.(mEpochName).raster;
                            rightTargRas = fat_raster(rasRightTarg, tickWidth);
                            rightTargRas = rightTargRas .* 2;
                            %                  colormap([1 1 1; cMap(2,:)])
                            imagesc(epochRangeDisplay, 1 : size(rasRightTarg, 1), rightTargRas(:, alignRightTarg + epochRangeDisplay))
                            
                            rasLeftTarg =  PlotData(kDataPlot).leftTarg.(mEpochName).raster;
                            leftTargRas = fat_raster(rasLeftTarg, tickWidth);
                            %                 colormap([1 1 1; cMap(1,:)])
                            imagesc(epochRangeDisplay, 1+size(rasRightTarg, 1) : size(rasRightTarg, 1) + size(rasLeftTarg, 1), leftTargRas(:, alignLeftTarg + epochRangeDisplay))
                            
                            plot(ax(axRas, mEpochPlot), [1 1], [0 rasYlim * .9], '-k', 'linewidth', 2)
                            
                            
                        case 'lfp'
                            % Right LFPs
                            if ~isempty(alignRightTarg)
                                lfpRightTarg = PlotData(kDataPlot).rightTarg.(mEpochName).signalMean;
                                plot(ax(axSig, mEpochPlot), epochRangeDisplay, lfpRightTarg(alignRightTarg + epochRangeDisplay), 'color', cMap(1,:), 'linewidth', targLineW)
                            end
                            % Left LFPs
                            if ~isempty(alignLeftTarg)
                                lfpLeftTarg = PlotData(kDataPlot).leftTarg.(mEpochName).signalMean;
                                plot(ax(axSig, mEpochPlot), epochRangeDisplay, lfpLeftTarg(alignLeftTarg + epochRangeDisplay), 'color', cMap(2,:), 'linewidth', targLineW)
                            end
                        case 'erp'
                            % Right LFPs
                            if ~isempty(alignRightTarg)
                                erpRightTarg = PlotData(kDataPlot).rightTarg.(mEpochName).signalMean;
                                plot(ax(axSig, mEpochPlot), epochRangeDisplay, erpRightTarg(alignRightTarg + epochRangeDisplay), 'color', cMap(1,:), 'linewidth', targLineW)
                            end
                            % Left LFPs
                            if ~isempty(alignLeftTarg)
                                erpLeftTarg = PlotData(kDataPlot).leftTarg.(mEpochName).signalMean;
                                plot(ax(axSig, mEpochPlot), epochRangeDisplay, erpLeftTarg(alignLeftTarg + epochRangeDisplay), 'color', cMap(2,:), 'linewidth', targLineW)
                            end
                            
                    end
                    
                end % mEpochPlot
                
                %                             legend(ax(axGo, 1), {num2cell(pSignalArray'), num2str(pSignalArray')})
                
                %         colorbar('peer', ax(axGo, 1), 'location', 'west')
                %         colorbar('peer', ax(axStopGo, 1), 'location', 'west')
                h=axes('Position', [0 0 1 1], 'Visible', 'Off');
                titleString = sprintf('%s \t %s', sessionID, PlotData(kDataPlot).name);
                text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top', 'color', 'k')
                if options.printPlot
                    localFigurePath = local_figure_path;
                    %          print(figureHandle,[localFigurePath, sessionID, '_', dataArray{kDataPlot}, '_mem_session_' dataType],'-dpdf', '-r300')
                    %             print(figureHandle,[localFigurePath, sessionID, '_', dataArray{kDataPlot}, '_mem_session_' dataType],'-djpeg')
%                     print(figureHandle,[local_figure_path, sessionID, '_mem_', PlotData(kDataPlot).name, '_',dataType,'.pdf'],'-dpdf', '-r300')
 micalaFolder = '/Volumes/SchallLab/Users/Paul/micala/mem/';
print(figureHandle,[micalaFolder, sessionID, '_mem_', PlotData(kDataPlot).name, '_',dataType,'.pdf'],'-dpdf', '-r300')
               end
            end % kUnitIndex
            
            
        end % plotFlag
        
    end
end
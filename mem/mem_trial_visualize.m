function mem_trial_visualize(subjectID, sessionID, options)

%
% function Unit = mem_session_data(subjectID, sessionID, dataType, varargin)
%
% Step through trials and visualize the data
%
% input:
%   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%
%   options: A structure with various ways to select/organize data: If
%   ccm_session_data.m is called without input arguments, the default
%   options structure is returned. options has the following fields with
%   possible values (default listed first):


%%
subjectID   = 'broca';
sessionID   = 'bp211n01';
jUnit       = 1;

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);

nTrial = size(trialData.trialOutcome, 1);

% CONSTANTS
TICK_WIDTH = 5;
sdfColor = [0 0 1];
% Set up plot
printPlot = true;
figureHandle = 87;
nRow = 7;
nColumn = 1;
if printPlot
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
else
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
end
axisHeight = axisHeight*1.4
clf

axEye = 1;
axRas = 2;
axH10 = 3;
axH20 = 4;
axG10 = 5;
axG20 = 6;
axPSP = 7;

% Eye trace axes
ax(axEye) = axes('units', 'centimeters', 'position', [xAxesPosition(axEye) yAxesPosition(axEye) axisWidth axisHeight]);

ax(axRas) = axes('units', 'centimeters', 'position', [xAxesPosition(axRas) yAxesPosition(axRas) axisWidth axisHeight]);
ax(axH10) = axes('units', 'centimeters', 'position', [xAxesPosition(axH10) yAxesPosition(axH10) axisWidth axisHeight]);
ax(axH20) = axes('units', 'centimeters', 'position', [xAxesPosition(axH20) yAxesPosition(axH20) axisWidth axisHeight]);
ax(axG10) = axes('units', 'centimeters', 'position', [xAxesPosition(axG10) yAxesPosition(axG10) axisWidth axisHeight]);
ax(axG20) = axes('units', 'centimeters', 'position', [xAxesPosition(axG20) yAxesPosition(axG20) axisWidth axisHeight]);
ax(axPSP) = axes('units', 'centimeters', 'position', [xAxesPosition(axPSP) yAxesPosition(axPSP) axisWidth axisHeight]);

for i = 1 : nTrial
    
    
    %   Eye Traces
    % ************************
    eyeX            = trialData.eyeX{i};
    eyeY            = trialData.eyeY{i};
    
    axes(ax(axEye))
    cla
    hold(ax(axEye), 'on')
    set(ax(axEye), 'xlim', [0 length(eyeX)])
    plot(ax(axEye), eyeX)
    plot(ax(axEye), eyeY)
    yS = ylim;
    plot(ax(axEye), [trialData.fixWindowEntered(i) trialData.fixWindowEntered(i)], [yS(1) yS(2)], 'k')
    plot(ax(axEye), [trialData.targOn(i) trialData.targOn(i)], [yS(1) yS(2)], 'k')
    plot(ax(axEye), [trialData.responseCueOn(i) trialData.responseCueOn(i)], [yS(1) yS(2)], 'k')
    plot(ax(axEye), [trialData.responseOnset(i) trialData.responseOnset(i)], [yS(1) yS(2)], 'k')
    plot(ax(axEye), [trialData.rewardOn(i) trialData.rewardOn(i)], [yS(1) yS(2)], 'k')
    
    
    %   Raster
    % ************************
    iRas            = spike_to_raster(trialData.spikeData{i,jUnit});
    iTick           = fat_raster(iRas, TICK_WIDTH);
    %    return
    axes(ax(axRas))
    cla
    hold(ax(axRas), 'on')
    set(ax(axRas), 'xlim', [0 length(eyeX)])
    colormap([1 1 1; sdfColor])
    imagesc(iTick)
    
    
    %   PSTH 10 ms
    % ************************
    binWidth        = 10;
    [PSTH, binCenters] = peristimulus_time_histogram(iRas, binWidth, 0);
    
    axes(ax(axH10))
    cla
    hold(ax(axH10), 'on')
    set(ax(axH10), 'xlim', [0 length(eyeX)])
    bar(binCenters, PSTH, 'hist');
    %     bar(binCenters, PSTH, 'hist', 'faceColor', sdfColor);
    
    
    %   PSTH 20 ms
    % ************************
    binWidth        = 20;
    [PSTH, binCenters] = peristimulus_time_histogram(iRas, binWidth, 0);
    
    axes(ax(axH20))
    cla
    hold(ax(axH20), 'on')
    set(ax(axH20), 'xlim', [0 length(eyeX)])
    bar(binCenters, PSTH, 'hist');
    %     bar(binCenters, PSTH, 'hist', 'color', sdfColor);
    
    
    %   SDF Gaussian 10 ms
    % ************************
    Kernel.method = 'gaussian';
    Kernel.sigma = 10;
    [sdf, kShape] = spike_density_function(iRas, Kernel);
    
    axes(ax(axG10))
    cla
    hold(ax(axG10), 'on')
    set(ax(axG10), 'xlim', [0 length(eyeX)])
    plot(ax(axG10), sdf, 'color', sdfColor)
    
    %   SDF Gaussian 20 ms
    % ************************
    Kernel.method = 'gaussian';
    Kernel.sigma = 20;
    [sdf, kShape] = spike_density_function(iRas, Kernel);
    
    axes(ax(axG20))
    cla
    hold(ax(axG20), 'on')
    set(ax(axG20), 'xlim', [0 length(eyeX)])
    plot(ax(axG20), sdf, 'color', sdfColor)
    
    %   PSP growth = 1ms;  decay = 20ms
    % ************************
    Kernel.method = 'postsynaptic potential';
    Kernel.growth = 1;
    Kernel.decay = 20;
    [sdf, kShape] = spike_density_function(iRas, Kernel);
    
    axes(ax(axPSP))
    cla
    hold(ax(axPSP), 'on')
    set(ax(axPSP), 'xlim', [0 length(eyeX)])
    plot(ax(axPSP), sdf, 'color', sdfColor)
    
    
    
savePlot = input('save?', 's');
if strcmp(savePlot, 'y')
         localFigurePath = local_figure_path;
         print(figureHandle,[localFigurePath, sprintf('%s_trial%.2d', sessionID, i)],'-dpdf', '-r300')
end
end

%%

clear Data

if nargin < 3
    options.dataType = 'neuron';
    
    options.figureHandle     = 2000;
    options.printPlot        = false;
    options.plotFlag         = true;
    options.collapseTarg      = false;
    options.filterData       = false;
    options.normalize        = false;
    options.unitArray        = 'each';
    
    if nargin == 0
        Data = options;
        return
    end
end
normalize       = options.normalize;
filterData      = options.filterData;
unitArray       = options.unitArray;
plotFlag        = options.plotFlag;
printPlot       = options.printPlot;
dataType        = options.dataType;
figureHandle    = options.figureHandle;



% Make sure user input a dataType that was recorded during the session
dataTypePossible = {'neuron', 'lfp', 'erp'};
if ~sum(strcmp(dataType, dataTypePossible))
    fprintf('%s Is not a valid data type \n', dataType)
    return
end

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);

if ~strcmp(SessionData.taskID, 'mem')
    fprintf('Not a memory guided saccade session, try again\n')
    return
end

Data.dataType = dataType;


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
        if strcmp(unitArray, 'each')
            unitArray     = dataArray;
        end
    case 'lfp'
        dataArray 	= num2cell(SessionData.lfpChannel);
        if strcmp(unitArray, 'each')
            unitArray 	= dataArray;
        end
    case 'erp'
        dataArray     = eeg_electrode_map(subjectID);
        if strcmp(unitArray, 'each')
            unitArray     = dataArray;
        end
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
                Data(kDataIndex).rightTarg.(mEpochName).alignTime = alignIndex;
                sdf = spike_density_function(alignedRasters, Kernel);
                if ~isempty(sdf); yMax(kDataIndex, mEpoch, 1) = nanmax(nanmean(sdf, 1)); end;
                Data(kDataIndex).rightTarg.(mEpochName).raster = alignedRasters;
                Data(kDataIndex).rightTarg.(mEpochName).signalFn = sdf;
                Data(kDataIndex).rightTarg.(mEpochName).signalMean = nanmean(sdf, 1);
                
                % Left Target trials
                [alignedRasters, alignIndex] = spike_to_raster(trialData.spikeData(leftTargTrial, kUnit), alignListL);
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
    
    
end % kUnitIndex









%%
%**********************   CLASSIFY CELL RESPONSE TYPE    %*************************
if strcmp(dataType, 'neuron')
    for kDataIndex = 1 : nUnit
        
        
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
    end % kUnitIndex
end








%%
%**********************   PLOTTING    %*************************


if plotFlag
    cMap = ccm_colormap([0 1]);
    kernelMethod = 'gaussian';
    SIGMA = 20;
    kernelMethod = 'postsynaptic potential';
    GROWTH = 1;
    DECAY = 20;
    % kernelIn = [SIGMA];
    
    
    
    
    targLineW = 2;
    
    
    for kDataIndex = 1 : nUnit
        %       colormap([1 1 1; cMap(1,:); cMap(2,:)])
        nRow = 2;
        nEpoch = length(epochArray);
        nColumn = nEpoch;
        figureHandle = figureHandle + 1;
        if printPlot
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
        else
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
        end
        clf
        for mEpoch = 1 : nEpoch
            mEpochName = epochArray{mEpoch};
            epochRangeDisplay = mem_epoch_range(mEpochName, 'plot');
            epochRangeAnalysis = mem_epoch_range(mEpochName, 'analyze');
            
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
            
            % _______  Set up axes  ___________
            % axes names
            axSig = 1;
            axRas = 2;
            
            % Set up plot axes
            % SDFs
            ax(axSig, mEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(axSig, mEpoch) yAxesPosition(axSig, mEpoch) axisWidth axisHeight]);
            set(ax(axSig, mEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [epochRangeDisplay(1) epochRangeDisplay(end)])
            cla
            hold(ax(axSig, mEpoch), 'on')
            title(epochArray{mEpoch})
            
            
            % Rasters
            ax(axRas, mEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(axRas, mEpoch) yAxesPosition(axRas, mEpoch) axisWidth axisHeight]);
            set(ax(axRas, mEpoch), 'ylim', [0 rasYlim], 'xlim', [epochRangeDisplay(1) epochRangeDisplay(end)])
            cla
            hold(ax(axRas, mEpoch), 'on')
            
            
            if mEpoch > 1
                set(ax(axSig, mEpoch), 'yticklabel', [])
                set(ax(axRas, mEpoch), 'yticklabel', [])
                
                set(ax(axSig, mEpoch), 'ycolor', [1 1 1])
                set(ax(axRas, mEpoch), 'ycolor', [1 1 1])
            end
            
            
            
            alignRightTarg = Data(kDataIndex).rightTarg.(mEpochName).alignTime;
            alignLeftTarg = Data(kDataIndex).leftTarg.(mEpochName).alignTime;
            
            axes(ax(axSig, mEpoch))
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
            plot(ax(axSig, mEpoch), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
            
            
            
            switch dataType
                case 'neuron'
                    % Right SDFs
                    if ~isempty(alignRightTarg)
                        sdfRightTarg = Data(kDataIndex).rightTarg.(mEpochName).signalMean;
                        plot(ax(axSig, mEpoch), epochRangeDisplay, sdfRightTarg(alignRightTarg + epochRangeDisplay), 'color', cMap(2,:), 'linewidth', targLineW)
                    end
                    % Left SDFs
                    if ~isempty(alignLeftTarg)
                        sdfLeftTarg = Data(kDataIndex).leftTarg.(mEpochName).signalMean;
                        plot(ax(axSig, mEpoch), epochRangeDisplay, sdfLeftTarg(alignLeftTarg + epochRangeDisplay), 'color', cMap(1,:), 'linewidth', targLineW)
                    end
                    
                    axes(ax(axRas, mEpoch))
                    colormap([1 1 1; cMap(1,:); cMap(2,:)])
                    rasRightTarg =  Data(kDataIndex).rightTarg.(mEpochName).raster;
                    rightTargRas = fat_raster(rasRightTarg, tickWidth);
                    rightTargRas = rightTargRas .* 2;
                    %                  colormap([1 1 1; cMap(2,:)])
                    imagesc(epochRangeDisplay, 1 : size(rasRightTarg, 1), rightTargRas(:, alignRightTarg + epochRangeDisplay))
                    
                    rasLeftTarg =  Data(kDataIndex).leftTarg.(mEpochName).raster;
                    leftTargRas = fat_raster(rasLeftTarg, tickWidth);
                    %                 colormap([1 1 1; cMap(1,:)])
                    imagesc(epochRangeDisplay, 1+size(rasRightTarg, 1) : size(rasRightTarg, 1) + size(rasLeftTarg, 1), leftTargRas(:, alignLeftTarg + epochRangeDisplay))
                    
                    plot(ax(axRas, mEpoch), [1 1], [0 rasYlim * .9], '-k', 'linewidth', 2)
                    
                    
                case 'lfp'
                    % Right LFPs
                    if ~isempty(alignRightTarg)
                        lfpRightTarg = Data(kDataIndex).rightTarg.(mEpochName).signalMean;
                        plot(ax(axSig, mEpoch), epochRangeDisplay, lfpRightTarg(alignRightTarg + epochRangeDisplay), 'color', cMap(1,:), 'linewidth', targLineW)
                    end
                    % Left LFPs
                    if ~isempty(alignLeftTarg)
                        lfpLeftTarg = Data(kDataIndex).leftTarg.(mEpochName).signalMean;
                        plot(ax(axSig, mEpoch), epochRangeDisplay, lfpLeftTarg(alignLeftTarg + epochRangeDisplay), 'color', cMap(2,:), 'linewidth', targLineW)
                    end
                case 'erp'
                    % Right LFPs
                    if ~isempty(alignRightTarg)
                        erpRightTarg = Data(kDataIndex).rightTarg.(mEpochName).signalMean;
                        plot(ax(axSig, mEpoch), epochRangeDisplay, erpRightTarg(alignRightTarg + epochRangeDisplay), 'color', cMap(1,:), 'linewidth', targLineW)
                    end
                    % Left LFPs
                    if ~isempty(alignLeftTarg)
                        erpLeftTarg = Data(kDataIndex).leftTarg.(mEpochName).signalMean;
                        plot(ax(axSig, mEpoch), epochRangeDisplay, erpLeftTarg(alignLeftTarg + epochRangeDisplay), 'color', cMap(2,:), 'linewidth', targLineW)
                    end
                    
            end
            
        end % mEpoch
        
        %                             legend(ax(axGo, 1), {num2cell(pSignalArray'), num2str(pSignalArray')})
        
        %         colorbar('peer', ax(axGo, 1), 'location', 'west')
        %         colorbar('peer', ax(axStopGo, 1), 'location', 'west')
        h=axes('Position', [0 0 1 1], 'Visible', 'Off');
        titleString = sprintf('%s \t %s', sessionID, Data(kDataIndex).name);
        text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top', 'color', 'k')
        if printPlot
            localFigurePath = local_figure_path;
            print(figureHandle,[localFigurePath, sessionID, '_', dataArray{kDataIndex}, '_mem_session_' dataType],'-dpdf', '-r300')
        end
    end % kUnitIndex
    
    
end % plotFlag


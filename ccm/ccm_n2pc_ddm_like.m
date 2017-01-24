function data = ccm_n2pc_ddm_like(subjectID, sessionID, varargin)
% DESCRIPTION
% Analyzes eeg data from an electrode (or from a list of electrodes) to
% determine whether the signal has "ddm-like" properties.
% See also ccm_ddm_like.m for original, neuronal spike- based analyses
%
%
% SYNTAX
% ddmLike = ccm_n2pc_ddm_like(subjectID, sessionID, varargin)
%
% % % % % % % INPUT:
% % % % % % %   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
% % % % % % %   sessionID: e.g. 'bp111n01', 'Allsaccade'
% % % % % % %
% % % % % % % alignEvent  - 	e.g. 'checkerOn', 'targOn'
% % % % % % %
% % % % % % % electrodeArray   - e.g. {'o1', o2'}
% % % % % % %
% % % % % % % OUTPUT:
% % % % % % % ddmLike       - a structure containg the following variables:
% % % % % % %

%
% EXAMPLE
% [ddmLike] = ccm_n2pc_ddm_like('Broca', 'bp063n01', 'alignEvent', 'targOn')
%
% .........................................................................

% CONTENTS
% 1.SET LB, UB, X0 VALUES
% 1.1.Starting value (z0)
% 1.1.1.GO units
% 1.1.2.STOP unit
% 1.2.Threshold (zc)
% 1.2.1.GO units
% 1.2.2.STOP unit
% 1.3.Accumulation rate correct (vCor)
% 1.3.1.GO units
% 1.3.2.STOP unit
% 1.4.Accumulation rate incorrect (vIncor)





% ccm_n2pc_ddm_like(subjectID, sessionID, varargin)
%
% Analyzes
% input:
%   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%
%   varargin: property names and their values:
%           'plotFlag': 0 or 1
%           'printPlot': 0 or 1: If set to 1, this prints the figure in the local_figures folder
%           'unitArray': a single unit, like 'spikeUnit17a', or an array of units, like {'spikeUnit17a', 'spikeUnit17b'}
%           'latencyMatchMethod': 'ssrt' or 'rt': do latency matching based on 'ssrt' or using the rt distributions to latcency match

%%
% subjectID = 'Broca';
% sessionID = 'bp063n01';
%%
% Set defaults
plotFlag    = 1;
printPlot   = 0;
figureHandle = 6000;
alignEvent = 'checkerOn';
electrodeArray  = eeg_electrode_map(subjectID);
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
        case 'electrodeArray'
            electrodeArray = varargin{i+1};
        otherwise
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELCARE CONSTNANTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epochOffset = 100;
preSaccadeBuffer = 40;

preSaccRejMS    = 20; % Used to reject EEG data 20 ms before a saccade begins
baseDuration    = -99;   % relative to start of subject fixating
n2pcEpoch       = [125 250];
rejectWindow    = baseDuration : n2pcEpoch(2);
satThreshold    = 50; % How many ms will we allow a signal to be saturated?
baseDuration    = -99;   % relative to start of subject fixating
displayRange    = [-100, 500];  % relative to
MIN_RT          = 150;
MAX_RT          = 1200;
STD_MULTIPLE    = 3;

% electrodeArray  = eeg_electrode_map(subjectID);
nElectrode      = length(electrodeArray);


channelArray = [];
for i = 1 : nElectrode
    channelArray = [channelArray, find(strcmp(electrodeArray{i}, electrodeArray))];
end





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND PROCESS DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID, [ccm_min_vars, 'eegData']);
pSignalArray    = ExtraVar.pSignalArray;
ssdArray        = ExtraVar.ssdArray;
pSignalArray    = pSignalArray(pSignalArray ~= .5);

if ~ismember('eegData', fieldnames(trialData))
    fprintf('Session %s does not contain eeg data \n', sessionID)
    ddmLike = [];
    return
end

% Make sure the task is correct
if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end




% make sure there are as many electrodes as we expect- otherwise, alert
% user and return with empty data
if size(trialData.trialOutcome.eegData, 2) < nElectrode - 1  % Broca had a severed electode at one point
    fprintf('%s %s should have %d electrodes but seems to have %d for the session', subjectID, sessionID, nElectrode, size(trialData.trialOutcome.eegData, 2));
    data = [];
    return
end




% Get rid of trials with outlying RTs
[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTD);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);
trialData(isnan(trialData.rt), :) = [];

signalLeftP = 100 .* pSignalArray(pSignalArray < .5);
signalRightP = 100 .* pSignalArray(pSignalArray > .5);

% Get default trial selection options
selectOpt       = ccm_trial_selection;
selectOpt.ssd       = 'none';
selectOpt.outcome     = {'goCorrectTarget'};


selectOpt.targDir   	= 'left';
selectOpt.rightCheckerPct = signalLeftP;
leftTrial = ccm_trial_selection(trialData, selectOpt);
selectOpt.targDir   	= 'right';
selectOpt.rightCheckerPct = signalRightP;
rightTrial = ccm_trial_selection(trialData, selectOpt);

% alignTime = trialData.(alignEvent);
% leftTrial(isnan(alignTime(leftTrial)) | isnan(trialData.rt(leftTrial))) = [];
% rightTrial(isnan(alignTime(rightTrial)) | isnan(trialData.rt(rightTrial))) = [];

% Rearrange trialData dataset with only the relevant trials
trialDataSub   = [trialData(leftTrial,:); trialData(rightTrial,:)];
nTrial      = size(trialData.trialOutcomeSub, 1);
leftTrial   = 1 : length(leftTrial);
rightTrial  = 1 + length(leftTrial) : nTrial;
alignTimeList = trialDataSub.(alignEvent);

%%




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP THROUGH THE CHANNELS, ANALYZING EACH
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iChannelInd = 1 : nElectrode;
    
    
    iChannel = channelArray(iChannelInd);
    
    leftTrialDelete = 0;
    rightTrialDelete = 0;
    
    % Go to Target trials
    [targEEG, alignIndex]    = align_signals(trialDataSub.eegData(:, iChannel), alignTimeList);
    
    
    % Loop through trials to process, delete, etc
    % ----------------------------------------
    satTrial = zeros(nTrial, 1);
    for iTrial = 1 : nTrial
        % Truncate signal just before saccade
        targEEG(iTrial, alignIndex + trialDataSub.rt(iTrial) - preSaccRejMS : end) = nan;
        
        % Reject trials with saturated signals
        if sum(diff(targEEG(iTrial, alignIndex + rejectWindow)) == 0) > satThreshold
            satTrial(iTrial) = 1;
            
            % Keep track of trials with saturated signals to re-adjust the
            % trial number (in the rearranged dataset "trialData".
            if sum(leftTrial == iTrial)
                leftTrialDelete = leftTrialDelete + 1;
            elseif sum(rightTrial == iTrial)
                rightTrialDelete = rightTrialDelete + 1;
            end
        end
        
        % Baseline-shift the signals
        iBase = mean(targEEG(iTrial, alignIndex + baseDuration : alignIndex));
        targEEG(iTrial, :)  = targEEG(iTrial, :) - iBase;
        
    end %  iTrial = 1 : nRightTrial
    
    % Separate eeg data into Left and Right target trials
    eegLeft = nanmean(targEEG(leftTrial,:), 1);
    eegRight = nanmean(targEEG(rightTrial,:), 1);
    
    % Delete rejected trials
    %    trialData.rt(find(satTrial))          = nan;
    targEEG(find(satTrial),:) = [];
    trialDataSub(find(satTrial),:) = [];
    alignTimeList(find(satTrial),:) = [];
    nTrial = size(trialData.trialOutcomeSub, 1);
    
    % Adjust the "new" leftTiral and rightTrial here
    leftTrial = 1 : length(leftTrial) - leftTrialDelete;
    rightTrial  = 1 + length(leftTrial) : length(leftTrial) + length(rightTrial) - rightTrialDelete;
    
    
    
    
    
    medianLeftRT     = round(nanmedian(trialDataSub.rt(leftTrial)));
    medianRightRT    = round(nanmedian(trialDataSub.rt(rightTrial)));
    epochEnd         = [medianLeftRT * ones(length(leftTrial), 1); medianRightRT * ones(length(rightTrial), 1)];
    epochEnd         = epochEnd - preSaccadeBuffer;
    % replace epoch-cutoffs for trials with rts shorter than the median RT
    epochEnd(leftTrial)  = alignIndex + min(trialDataSub.rt(leftTrial) - preSaccadeBuffer, epochEnd(leftTrial));
    epochEnd(rightTrial) = alignIndex + min(trialDataSub.rt(rightTrial) - preSaccadeBuffer, epochEnd(rightTrial));
    epochBegin           = alignIndex + epochOffset * ones(nTrial, 1);
    epochDuration        = epochEnd - epochBegin;
    t = epochEnd < epochBegin;
    
    
    % For ddm-like consideration, thake the mean of the signal within the
    % epoch that would've been used for spikes (see ccm_n2pc_like.m)
    eegMeanEpoch = cellfun(@(x,y,z) nanmean(x(y:z)), num2cell(targEEG,2), num2cell(epochBegin), num2cell(epochEnd), 'uniformoutput', false);
    eegMeanEpoch = cell2mat(eegMeanEpoch);
    
    %%
    
    %     data(iChannel).spikeRate      = spikeRate;
    data(iChannel).leftTrial      = leftTrial;
    data(iChannel).rightTrial     = rightTrial;
    data(iChannel).signalP        = trialDataSub.targ1CheckerProp;
    %     data(iChannel).alignedRasters = alignedRasters;
    data(iChannel).alignIndex     = alignIndex;
    data(iChannel).epochOffset    = epochOffset;
    data(iChannel).eegMeanEpoch    = eegMeanEpoch;
    data(iChannel).alignedSignal   = targEEG;
    data(iChannel).rt   = trialDataSub.rt;
    
    
    [ddmData]                    = ding_gold_ddm_like(data(iChannel), 'dataType', 'eeg');
    choiceDependent(iChannel)    = ddmData.choiceDependent;
    coherenceDependent(iChannel) = ddmData.coherenceDependent;
    ddmLike(iChannel)            = ddmData.ddmLike;
    
    % ________________________________________________________________
    % PLOT THE DATA
    
    if plotFlag
        
        % SET UP PLOT
        lineW = 2;
        plotEpochRange = [-200 : 300];
        cMap = ccm_colormap(pSignalArray);
        leftColor = cMap(1,:) .* .8;
        rightColor = cMap(end,:) .* .8;
        nRow = 2;
        nColumn = 3;
        figureHandle = figureHandle + 1;
        if printPlot
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
        else
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
        end
        clf
        axChoice = 1;
        axCoh = 2;
        axCohL = 3;
        axCohR = 4;
        
        
        
        ax(axChoice) = axes('units', 'centimeters', 'position', [xAxesPosition(axChoice, 2) yAxesPosition(axChoice, 2) axisWidth axisHeight]);
        cla
        hold(ax(axChoice), 'on')
        switch choiceDependent(iChannel)
            case true
                choiceStr = 'YES';
            otherwise
                choiceStr = 'NO';
        end
        tt = sprintf('Choice dependence: %s', choiceStr);
        title(tt)
        
        ax(axCohL) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 1) yAxesPosition(2, 1) axisWidth axisHeight]);
        cla
        hold(ax(axCohL), 'on')
        title('Coherence dependence')
        
        ax(axCohR) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 3) yAxesPosition(2, 3) axisWidth axisHeight]);
        cla
        hold(ax(axCohR), 'on')
        title('Coherence dependence')
        
        ax(axCoh) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2) yAxesPosition(2, 2) axisWidth axisHeight]);
        cla
        hold(ax(axCoh), 'on')
        switch coherenceDependent(iChannel)
            case true
                cohStr = 'YES';
            otherwise
                cohStr = 'NO';
        end
        tt = sprintf('Coherence dependence: %s', cohStr);
        title(tt)
        
        
        signalMax = max(max(nanmean(eegLeft(:, alignIndex + plotEpochRange), 1)), max(nanmean(eegRight(:, alignIndex + plotEpochRange), 1)));
        yMax = 1.1 * signalMax;
        signalMin = min(min(nanmean(eegLeft(:, alignIndex + plotEpochRange), 1)), min(nanmean(eegRight(:, alignIndex + plotEpochRange), 1)));
        yMin = 1.1 * signalMin;
        fillX = [epochOffset, nanmean(epochEnd)-alignIndex, nanmean(epochEnd)-alignIndex, epochOffset];
        fillY = [yMin yMin yMax yMax];
        fillColor = [1 1 .5];
        
        
        
        
        
        
        
        
        % CHOICE DEPENDENCE PLOTTING(LEFT VS. RIGHT CHOICE FOR CORRECT TRIALS)
        axes(ax(axChoice))
        h = fill(fillX, fillY, fillColor);
        set(h, 'edgecolor', 'none');
        plot(ax(axChoice), plotEpochRange, nanmean(eegLeft(:, alignIndex + plotEpochRange), 1), 'color', leftColor, 'linewidth', lineW)
        plot(ax(axChoice), plotEpochRange, nanmean(eegRight(:, alignIndex + plotEpochRange), 1), 'color', rightColor, 'linewidth', lineW)
        plot(ax(axChoice), [1 1], [yMin yMax], '-k', 'linewidth', 2);
        set(ax(axChoice), 'ylim', [yMin yMax])
        
        
        
        
        
        
        
        % COHERENCE DEPENDENCE PLOTTING
        
        % Leftward trials
        
        plot(ax(axCohL), [0 0], [yMin yMax], '-k', 'linewidth', 2);
        axes(ax(axCohL))
        h = fill(fillX, fillY, fillColor);
        set(h, 'edgecolor', 'none');
        
        for i = 1 : length(signalLeftP)
            iProp = pSignalArray(i);
            
            % Determine color to use for plot based on which checkerboard color
            % proportion being used. Normalize the available color spectrum to do
            % it
            sigColor = cMap(i,:);
            
            leftTrialData = trialDataSub(leftTrial,:);
            signalTrial = leftTrial(leftTrialData.targ1CheckerProp == iProp);
            % Go to Target trials
            iEegLeft = nanmean(targEEG(signalTrial, :), 1);
            
            plot(ax(axCohL), plotEpochRange, iEegLeft(alignIndex + plotEpochRange), 'color', sigColor, 'linewidth', lineW)
            set(ax(axCohL), 'ylim', [yMin yMax])
            
            scatter(ax(axCoh), trialDataSub.targ1CheckerProp(signalTrial), eegMeanEpoch(signalTrial), 'o', 'markeredgecolor', sigColor, 'markerfacecolor', sigColor, 'sizedata', 20)
        end % for i = 1 : length(signalLeftP)
        
        
        
        % Rightward trials
        
        plot(ax(axCohR), [0 0], [yMin yMax], '-k', 'linewidth', 2);
        axes(ax(axCohR))
        h = fill(fillX, fillY, fillColor);
        set(h, 'edgecolor', 'none');
        for i = (i+1) : (length(signalLeftP) + length(signalRightP))
            iProp = pSignalArray(i);
            
            % Determine color to use for plot based on which checkerboard color
            % proportion being used. Normalize the available color spectrum to do
            % it
            sigColor = cMap(i,:);
            
            rightTrialData = trialDataSub(rightTrial,:);
            signalTrial = rightTrial(rightTrialData.targ1CheckerProp == iProp);
            % Go to Target trials
            iEegRight = nanmean(targEEG(signalTrial, :), 1);
            
            plot(ax(axCohR), plotEpochRange, iEegRight(alignIndex + plotEpochRange), 'color', sigColor, 'linewidth', lineW)
            set(ax(axCohR), 'ylim', [yMin yMax])
            
            scatter(ax(axCoh), trialDataSub.targ1CheckerProp(signalTrial), eegMeanEpoch(signalTrial), 'o', 'markeredgecolor', sigColor, 'markerfacecolor', sigColor, 'sizedata', 30)
        end % for i = 1 : length(signalRightP)
        
        
        
        
        % regressions on trial-by-trial signal metric in the epoch
        xLeft = (signalLeftP(1) : .01 : signalLeftP(end));
        xRight = (signalRightP(1) : .01 : signalRightP(end));
        switch ddmData.leftIsIn
            case true
                yLeft = ddmData.pIn(1) .* xLeft + ddmData.pIn(2);
                yRight = ddmData.pOut(1) .* xRight + ddmData.pOut(2);
            case false
                yRight = ddmData.pIn(1) .* xRight + ddmData.pIn(2);
                yLeft = ddmData.pOut(1) .* xLeft + ddmData.pOut(2);
        end
        plot(ax(axCoh), xLeft, yLeft, '-k', 'lineWidth', lineW)
        plot(ax(axCoh), xRight, yRight, '-k', 'lineWidth', lineW)
        set(ax(axCoh), 'Xlim', [.9 * signalLeftP(1) 1.1 * signalRightP(end)])
        set(ax(axCoh), 'xtick', pSignalArray)
        set(ax(axCoh), 'xtickLabel', pSignalArray*100)
        
        
        h=axes('Position', [0 0 1 1], 'Visible', 'Off');
        if choiceDependent(iChannel) && coherenceDependent(iChannel)
            ddmStr = 'YES';
        else
            ddmStr = 'NO';
        end
        titleString = sprintf('%s \t %s \t DDM-Like: %s', sessionID, electrodeArray{iChannel}, ddmStr);
        text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')
        
        if printPlot
            localFigurePath = local_figure_path;
            print(figureHandle,[localFigurePath, sessionID, '_ccm_ddm_like_', electrodeArray{iChannel}],'-dpdf', '-r300')
        end
    end % if plotFlag
    
end






end
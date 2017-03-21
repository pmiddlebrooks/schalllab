function unitInfo = ccm_neuron_choice(subjectID, sessionID, unitArray, options)
%
%
% function unitInfo = ccm_neuron_choice(subjectID, sessionID, unitArray, options)
%
% Classify neurons according to modified criteria of Ding & Gold 2012 Cereb Ctx.
% They query "choice" and "coherence" dependence within 4 defined epochs:
% Stim, Sacc, Post, Reward.


%%
% Set default options or return a default options structure
if nargin < 4
    
    options.multiUnit            = false;
    options.plotFlag            = true;
    options.printPlot           = true;
    options.figureHandle           = 54;
    
    options.ms2Std = 75;
    % Return just the default options struct if no input
    if nargin == 0
        unitInfo           = options;
        return
    end
end


% ________________________________________________________________
% DELCARE CONSTNANTS AND PREPARE DATA
% Constants:

% These will determine the trial-to-trial epochs used for analyses:
epochOffset = 150;  % When to begin spike rate analysis after stimulus (checkerboard) onset
preSaccadeBuffer = 10; % When to cut off spike rate analysis before saccade onset
minEpochDuration = 10; % Only include trials for which the determined epoch is this long

epochRangeChecker = 1 : 1000;  % Make this big enough to catch really late cancel times.

% For now use a convolved spike density function.. but see Ding and Gold,
% who do calculations on an "unfiltered spike density function".
Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;

% Initialize choice and coherence dependence to null assumption (i.e. they
% aren't)
choiceDependent     = false;
coherenceDependent  = false;
ddmLike             = false;

% some constants
alphaChoice     = .05;   % alpha criteria for choice dependence
alphaCoherence  = .05;   % alpha criteria for coherence dependence





% ________________________________________________________________
% Load the data


    
    
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID, [ccm_min_vars, 'spikeData'], options.multiUnit);
pSignalArray = ExtraVar.pSignalArray;
pSignalArray = pSignalArray(pSignalArray ~= .5);

% Get rid of all trials except GoCorrect Target trials for this analysis,
% with targetCheckerProp above or below .5
keepTrial = strcmp(trialData.trialOutcome, 'goCorrectTarget');
trialData = structfun(@(x) x(keepTrial,:), trialData, 'uni', false);
keepCoh = ~(trialData.targ1CheckerProp == .5);
trialData = structfun(@(x) x(keepCoh,:), trialData, 'uni', false);
[a, spikeUnit] = ismember(unitArray, SessionData.spikeUnitArray);


% RT must be after the sum of the following to be included in analyses:
MIN_RT = 120;
MAX_RT = 1200;
nSTD = 3;
% Get rid of trials with outlying RTs
[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTD);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);
trialData = structfun(@(x) x(~isnan(trialData.rt),:), trialData, 'uni', false);
nTrial = length(trialData.rt);

signalLeftP = pSignalArray(pSignalArray < .5);
signalRightP = pSignalArray(pSignalArray > .5);


leftTrial = trialData.targ1CheckerProp < .5;
rightTrial = trialData.targ1CheckerProp > .5;
easyLeftTrial = trialData.targ1CheckerProp == min(pSignalArray);
easyRightTrial = trialData.targ1CheckerProp == max(pSignalArray);

medianLeftRT = nanmedian(trialData.rt(easyLeftTrial));
medianRightRT = nanmedian(trialData.rt(easyRightTrial));


for iUnit = 1 : length(spikeUnit)
    fprintf('%s: %s\n',sessionID, (unitArray{iUnit}))
    
    %   Get the saccade receptive field
    %   Get neural data from the session/unit:
    optCollapse             = ccm_options;
    optCollapse.plotFlag    = false;
    optCollapse.printFlag    = false;
    optCollapse.collapseTarg = true;
    optCollapse.collapseSignal = true;
    optCollapse.unitArray   = unitArray(iUnit);
    optCollapse.doStops   = false;
    optCollapse.multiUnit   = options.multiUnit;
    UnitCollapse                = ccm_session_data(subjectID, sessionID, optCollapse);
    rf = ccm_find_saccade_rf(UnitCollapse);

    
    
    % Go to Target trials
    alignmentTimeList = trialData.checkerOn;
    [alignedRasters, alignmentIndex] = spike_to_raster(trialData.spikeData(:, spikeUnit(iUnit)), alignmentTimeList);
    sdf = spike_density_function(alignedRasters, Kernel);
    alignedRasters = num2cell(alignedRasters, 2);
    
    
    
    %    CHOICE DEPENDENCE
    % =================================================================
    
    
    % Choice Selection Time -- via Hanes et al 1998 (p.822) differential sdf test
    % ------------------------------
    % Figure out when signals diverge w.r.t. choice direction. Use the
    % easiest left and right coherence, so the epoch begins as early as
    % possbile (to capture as much coherence divergence during
    % coherence-dependence analysis
    choiceSelectionTime = nan; % Initialize to NaN
    
    goTargEasyRightFn      = nanmean(sdf(easyRightTrial,:), 1);
    goTargEasyLeftFn      = nanmean(sdf(easyLeftTrial,:), 1);
    
    rightLeftDiffFn = goTargEasyRightFn - goTargEasyLeftFn;
    
    sdfDiffCheckerOn = 500;
    % Differential sdf go - stopStop
    
    % Baseline sdf differential to caclulate SDF
    %             stdDiffFn = nanmean(sdf(rightTrial,:), 1) - nanmean(sdf(leftTrial,:), 1);
    sdfDiffBase = rightLeftDiffFn(alignmentIndex + (-sdfDiffCheckerOn+1 : epochRangeChecker(end)))';
    
    % Get the standard deviation of the differential sdf during the
    % 500 ms before checkerboard onset
    stdDiff = std(sdfDiffBase(1:sdfDiffCheckerOn));
    
    % GoTarg - stopStop SDF from checker Onset to end of checker epoch.
    sdfDiff = abs(rightLeftDiffFn(alignmentIndex : alignmentIndex + epochRangeChecker(end))');
    
    
    % are there times at which the difference between sdfs is
    % greater than 2 standard deviations of the difference 500
    % ms before checkerboard onset? Check from checkerboard onset
    % to end of the checkerboard epoch.
    % So std2Ind starts at checkerboard onset.
    std2Ind = sdfDiff > 2*stdDiff;
    
    % Look for a sequence of options.ms2Std ms for which the go sdf is 2
    % std greater than the stop sdf. Determein whether there was a time
    % after the checkerboard onset that the differential
    % sdf was > 2*Std for at least options.ms2Std ms.
    riseAbove2Std = find([0; diff(std2Ind)] == 1);
    sinkBelow2Std = find([0; diff(std2Ind)] == -1);
    if ~isempty(riseAbove2Std)
        % Get rid of occasions for which the signals differ
        % going into the epoch (and therefore they will
        % cease to differ before they begin again to
        % differ)
        if ~isempty(sinkBelow2Std)
            sinkBelow2Std(sinkBelow2Std < riseAbove2Std(1)) = [];
        end
        
        % If there's one more riseAbove2Std than sinkBelow2Std, the last riseAbove2Std
        % will last until the end of the sdf: Add to
        % singkBelowStd the end of the epoch
        if length(riseAbove2Std) > length(sinkBelow2Std)
            sinkBelow2Std = [sinkBelow2Std; epochRangeChecker(end)];
        end
        
        % Now riseAbove2Std length should be equal. See if
        % any of the riseAbove2Std streaks go longer than
        % 50ms
        ind = find(sinkBelow2Std - riseAbove2Std >= options.ms2Std, 1);
        if ~isempty(ind)
            choiceSelectionTime = riseAbove2Std(ind);
        end
    end
    
    
    
    
    % ________________________________________________________________
    % DETERMINE THE ONSET OF THE CHOICE DEPENDENCE a la Ding & Gold
    tChoice = nan;  % Initialize to nan; in the case of eeg signals there may not be a tChoice (noise in signal)
    nanThreshold = .3; % Establish a threshold of non-signal fraction over which we cease looking for tChoice
    slideWindowWidth = 50; % ms, Ding and Gold used a 100ms sliding window
    slideWindowStep = 10; % ms, D&G used 10 ms steps
    
    
    choiceDependenceFound = false;
    iStepInd = 0;
    while ~choiceDependenceFound
        
        
        iEpochBegin = (alignmentIndex + (iStepInd * slideWindowStep) + epochOffset) * ones(nTrial, 1);
        iEpochEnd = iEpochBegin + slideWindowWidth;
        epochDuration = iEpochEnd - iEpochBegin;
        
        
        if iEpochEnd(1) > length(alignedRasters{1})
            break
        end
%         switch dataType
%             case 'spikes'
                nSpike = cellfun(@(x,y,z) sum(x(y:z)), alignedRasters, num2cell(iEpochBegin), num2cell(iEpochEnd), 'uniformoutput', false);
                iSpikeRate = cell2mat(nSpike) .* 1000 ./ epochDuration;
                % Choice dependence
                leftMetric = iSpikeRate(leftTrial);
                rightMetric = iSpikeRate(rightTrial);
%             case 'eeg'
%                 eegMeanEpoch = cellfun(@(x,y,z) nanmean(x(y:z)), num2cell(alignedSignal,2), num2cell(iEpochBegin), num2cell(iEpochEnd), 'uniformoutput', false);
%                 eegMeanEpoch = cell2mat(eegMeanEpoch);
%                 % Choice dependence
%                 leftMetric = eegMeanEpoch(leftTrial);
%                 rightMetric = eegMeanEpoch(rightTrial);
%             otherwise
%         end
        
        % Break out of the loop if we're so far out that tChoice is
        % meaningless
        if (sum(isnan(leftMetric)) / length(leftMetric) > nanThreshold || ...
                sum(isnan(rightMetric)) / length(rightMetric) > nanThreshold)
            break
        end
        
        
        [p, h, stats] = ranksum(leftMetric, rightMetric);
        
        if p < alphaChoice
            choiceDependenceFound = true;
            tChoice = (iStepInd * slideWindowStep) + epochOffset;
        end
        iStepInd = iStepInd + 1;
    end
    
    
    
    
    
    
    
    % Define the beginning and end of the epoch to analyze
    % Initialize end of epoch as median RTs from easiect choice
    % conditions
    epochEnd = ceil([medianLeftRT * ones(sum(leftTrial), 1); medianRightRT * ones(sum(rightTrial), 1)]);
    % Use Choice selection time if possible, to define
    % beginning of epoch. Otherwise, use the RTs
    if ~isnan(choiceSelectionTime)
        epochBegin = choiceSelectionTime * ones(nTrial, 1);
    else
        % epochBegin = epochEnd - epochDuration;
        epochBegin = ceil(.5 * epochEnd);
    end
    
    % Replace epoch-cutoffs for trials with rts shorter than the median RT
    earlyRTTrial = trialData.rt < epochEnd + preSaccadeBuffer;
    epochEnd(earlyRTTrial) = trialData.rt(earlyRTTrial) - preSaccadeBuffer;
    
    % Adjust for the alignment index
    epochEnd = epochEnd + alignmentIndex;
    epochBegin = epochBegin + alignmentIndex;
    
    % If there are trials with negative epochs because of the
    negativeEpochTrial = epochEnd < epochBegin + minEpochDuration;
    epochBegin(negativeEpochTrial) = epochEnd(negativeEpochTrial) - minEpochDuration;
    
    epochDuration = epochEnd - epochBegin;
    %    excludeTrial = find(epochDuration < minEpochDuration);
    %    leftTrial(ismember(leftTrial, excludeTrial))= [];
    %    rightTrial(ismember(rightTrial, excludeTrial))= [];
    
    nSpike = cellfun(@(x,y,z) sum(x(y:z)), alignedRasters, num2cell(epochBegin), num2cell(epochEnd), 'uniformoutput', false);
    spikeRate = cell2mat(nSpike) .* 1000 ./ epochDuration;
    
    
    
    
    
    % Choice dependence
    
    [p, h, stats]   = ranksum(spikeRate(leftTrial), spikeRate(rightTrial));
    
    if p < alphaChoice
        choiceDependent = true;
    end
    
    if nanmedian(spikeRate(leftTrial)) >= nanmedian(spikeRate(rightTrial))
        leftIsIn    = true;
        inTrial     = leftTrial;
        outTrial    = rightTrial;
    else
        leftIsIn    = false;
        inTrial     = rightTrial;
        outTrial    = leftTrial;
    end
    
    
    
    
    %    COHERENCE DEPENDENCE
    % =================================================================
    
    % For IN trials
    
    % Regress spikeRate vs signalStrength into RF
    [coeffIn, sIn]          = polyfit(trialData.targ1CheckerProp(inTrial), spikeRate(inTrial), 1);
    [yPredIn, deltaIn]  = polyval(coeffIn, trialData.targ1CheckerProp(inTrial), sIn);
    statsIn             = regstats(trialData.targ1CheckerProp(inTrial), spikeRate(inTrial));
    rIn                 = corr(trialData.targ1CheckerProp(inTrial), spikeRate(inTrial));
    
    slopeIn     = coeffIn(1);
    signSlopeIn = sign(slopeIn);
    fTestIn     = statsIn.fstat.f;
    pValIn      = statsIn.fstat.pval;
    
    
    % Regress spikeRate vs signalStrength out of RF
    [coeffOut, sOut]        = polyfit(trialData.targ1CheckerProp(outTrial), spikeRate(outTrial), 1);
    [yPredOut, deltaOut] = polyval(coeffOut, trialData.targ1CheckerProp(outTrial), sOut);
    statsOut            = regstats(trialData.targ1CheckerProp(outTrial), spikeRate(outTrial));
    rOut                = corr(trialData.targ1CheckerProp(outTrial), spikeRate(outTrial));
    
    slopeOut    = coeffOut(1);
    signSlopeOut = -sign(slopeOut);  % NOTE THIS IS NEGATED BECAUSE IF IN CONDITION GOES HARD TO EASY, OUT GOES EASY TO HARD, VICE VERSA
    fTestOut    = statsOut.fstat.f;
    pValOut     = statsOut.fstat.pval;
    
    
    % Decision tree to determine whether the neuron/signal was "coherence dependent"
    if pValIn < alphaCoherence
        if pValOut > alphaCoherence
            coherenceDependent = true;
        elseif pValOut < alphaCoherence
            % slopeOut must have opposite sign than slopeIn
            if signSlopeIn ~= signSlopeOut
                coherenceDependent = true;
            end
        end
    elseif pValIn > alphaCoherence
        if pValOut < alphaCoherence
            coherenceDependent = true;
        elseif pValOut > alphaCoherence
            % slopeOut must have opposite sign than slopeIn
            if signSlopeIn ~= signSlopeOut
                coherenceDependent = true;
            end
        end
    end
    
    if choiceDependent && coherenceDependent
        ddmLike = true;
    end
    
    
    
  
    
    
    
    
    
        %%
    
    % ________________________________________________________________
    % PLOT THE DATA
    
    if options.plotFlag
        
        sdfLeft = nanmean(sdf(leftTrial, :), 1);
        sdfRight = nanmean(sdf(rightTrial, :), 1);
        
        % SET UP PLOT
        lineW = 2;
        plotEpochRange = [-200 : 300];
        plotEpochRange = [-49 : 450];
        cMap = ccm_colormap(pSignalArray);
        leftColor = cMap(1,:) .* .8;
        rightColor = cMap(end,:) .* .8;
        nRow = 2;
        nColumn = 3;
        options.figureHandle = options.figureHandle + 1;
        if options.printPlot
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, options.figureHandle);
        else
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, options.figureHandle);
        end
        clf
        axChoice = 1;
        axCoh = 2;
        axCohL = 3;
        axCohR = 4;
        
        
        
        ax(axChoice) = axes('units', 'centimeters', 'position', [xAxesPosition(axChoice, 2) yAxesPosition(axChoice, 2) axisWidth axisHeight]);
        cla
        hold(ax(axChoice), 'on')
        switch choiceDependent(iUnit)
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
        
        ax(axCoh) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2)*1.3 yAxesPosition(2, 2)*.8 axisWidth/2 axisHeight*.8]);
        cla
        hold(ax(axCoh), 'all')
        switch coherenceDependent(iUnit)
            case true
                cohStr = 'YES';
            otherwise
                cohStr = 'NO';
        end
        tt = sprintf('Coherence dependence: %s', cohStr);
        title(tt)
        
        
        sdfMax = max(max(sdfLeft(alignmentIndex + plotEpochRange)), max(sdfRight(alignmentIndex + plotEpochRange)));
        yMax = 1.1 * sdfMax;
        %       fillX = [mean(epochBegin)-alignmentIndex, mean(epochEnd)-alignmentIndex, mean(epochEnd)-alignmentIndex, mean(epochBegin)-alignmentIndex];
        fillXLeft = [mean(epochBegin(leftTrial))-alignmentIndex, mean(epochEnd(leftTrial))-alignmentIndex, mean(epochEnd(leftTrial))-alignmentIndex, mean(epochBegin(leftTrial))-alignmentIndex];
        fillXRight = [mean(epochBegin(rightTrial))-alignmentIndex, mean(epochEnd(rightTrial))-alignmentIndex, mean(epochEnd(rightTrial))-alignmentIndex, mean(epochBegin(rightTrial))-alignmentIndex];
        fillY = [.1 .1 yMax yMax];
        fillColor = [1 1 .5];
        
        % CHOICE DEPENDENCE PLOTTING(LEFT VS. RIGHT CHOICE FOR CORRECT TRIALS)
        axes(ax(axChoice))
        %       h = fill(fillX, fillY, fillColor);
        %       set(h, 'edgecolor', 'none');
        plot(ax(axChoice), plotEpochRange, sdfLeft(alignmentIndex + plotEpochRange), 'color', leftColor, 'linewidth', lineW)
        plot(ax(axChoice), plotEpochRange, sdfRight(alignmentIndex + plotEpochRange), 'color', rightColor, 'linewidth', lineW)
        plot(ax(axChoice), [1 1], [0 yMax], '-k', 'linewidth', 2);
        set(ax(axChoice), 'xlim', [plotEpochRange(1) plotEpochRange(end)], 'ylim', [0 yMax])
        
        
        
        
        
        
        
        % COHERENCE PLOTTING
        
        minColorGun = 0;
        maxColorGun = 1;
        epochRange = ccm_epoch_range('checkerOn', 'plot');
        % Leftward trials
        
        plot(ax(axCohL), [0 0], [0 yMax], '-k', 'linewidth', 2);
        axes(ax(axCohL))
        h = fill(fillXLeft, fillY, fillColor);
        set(h, 'edgecolor', 'none');
        
%         % Rightward trials

plot(ax(axCohR), [0 0], [0 yMax], '-k', 'linewidth', 2);
        axes(ax(axCohR))
        h = fill(fillXRight, fillY, fillColor);
        set(h, 'edgecolor', 'none');

        
         for i = 1 : length(pSignalArray)
             iProp = pSignalArray(i);
            
            % Determine color to use for plot based on which checkerboard color
            % proportion being used. Normalize the available color spectrum to do
            % it
            inhColor = cMap(i,:);
            
            iTrial = trialData.targ1CheckerProp == iProp;
            iSdf = nanmean(sdf(iTrial,:), 1);
            if iProp < .5
            plot(ax(axCohL), plotEpochRange, iSdf(alignmentIndex + plotEpochRange), 'color', inhColor, 'linewidth', lineW)
            set(ax(axCohL), 'xlim', [plotEpochRange(1) plotEpochRange(end)], 'ylim', [0 yMax])
            elseif iProp > .5
             plot(ax(axCohR), plotEpochRange, iSdf(alignmentIndex + plotEpochRange), 'color', inhColor, 'linewidth', lineW)
            set(ax(axCohR), 'xlim', [plotEpochRange(1) plotEpochRange(end)], 'ylim', [0 yMax])
            end
         end
             
             
       
        boxplot(ax(axCoh), spikeRate, trialData.targ1CheckerProp, 'position', pSignalArray, 'colors', cMap, 'plotstyle', 'compact')
        
        % regressions on trial-by-trial spike rates in the epoch
        xLeft = (signalLeftP(1) : .01 : signalLeftP(end));
        xRight = (signalRightP(1) : .01 : signalRightP(end));
        switch leftIsIn
            case true
                yLeft = coeffIn(1) .* xLeft + coeffIn(2);
                yRight = coeffOut(1) .* xRight + coeffOut(2);
            case false
                yRight = coeffIn(1) .* xRight + coeffIn(2);
                yLeft = coeffOut(1) .* xLeft + coeffOut(2);
        end
        plot(ax(axCoh), xLeft, yLeft, '-k', 'lineWidth', lineW)
        plot(ax(axCoh), xRight, yRight, '-k', 'lineWidth', lineW)
        set(ax(axCoh), 'Xlim', [signalLeftP(1)-.02 signalRightP(end)+.02])
        set(ax(axCoh), 'xtick', pSignalArray)
        set(ax(axCoh), 'xtickLabel', pSignalArray*100)
        
        
        h=axes('Position', [0 0 1 1], 'Visible', 'Off');
        if choiceDependent(iUnit) && coherenceDependent(iUnit)
            ddmStr = 'YES';
        else
            ddmStr = 'NO';
        end
        titleString = sprintf('%s\t %s\t DDM-Like: %s', sessionID, SessionData.spikeUnitArray{spikeUnit(iUnit)}, ddmStr);
        text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')
        
        if options.printPlot
            if ~isdir(fullfile(local_figure_path, subjectID, 'choice'))
                mkdir(fullfile(local_figure_path, subjectID, 'choice'))
            end
            
            print(options.figureHandle,fullfile(local_figure_path, subjectID, 'choice', [sessionID, '_', SessionData.spikeUnitArray{spikeUnit(iUnit)}, '_ccm_ddm_like', '.pdf']),'-dpdf', '-r300')
        end
    end % if options.plotFlag

    
    
    
end

unitInfo.rf                     = rf;
unitInfo.choiceDependent        = choiceDependent;
unitInfo.coherenceDependent     = coherenceDependent;
unitInfo.ddmLike                = ddmLike;
unitInfo.tChoice                = tChoice;
unitInfo.choiceSelectionTime    = choiceSelectionTime;
unitInfo.leftIsIn               = leftIsIn;
unitInfo.coeffIn                = coeffIn;
unitInfo.coeffOut               = coeffOut;




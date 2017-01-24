function Data = ccm_inhibition_vs_accuracy(subjectID, sessionID, options)

% IN THE WORKS.... for calculaing separate ssrts, should I use the entire
% go RT distribution both for noncanceled stop correct and error choices,
% or use correct and error go RT distributions for each??


%
% function data = ccm_inhibition(subjectID, sessionID, plotFlag, figureHandle)
%
% Response inhibition analyses for choice countermanding task.
%
%
% If called without any arguments, returns a default options structure.
% If options are input but one is not specified, it assumes default.
%
% input:
%   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%
% Possible options are (default listed first):
%     options.collapseSignal    = Collapse across signal strength (difficulty conditions)?
%            false, true
%     options.collapseTarg 	= collapse angle/directions of the CORRECT
%     TARGET WITHIN a signal strength (so for signal strengths with correct
%     targets on the left, all left targets will be treated as one if set
%     to true
%           false, true
%     options.include50 	= if there is a 50% signal condition, do you
%           want to include it in analyses?
%           false, true
%
%     options.plotFlag       = true, false;
%     options.printPlot       = false, true;
%     options.figureHandle  = optional way to assign the figure to a handle
%
%
% Returns data structure with fields:
%
%   nGo
%   nGoRight
%   nStopIncorrect
%   nStopIncorrectRight
%   goRightLogical
%   goRightSignalStrength
%   stopRightLogical
%   stopRightSignalStrength


% Set default options or return a default options structure
if nargin < 3
    options.collapseSignal   	= false;
    options.collapseTarg        = false;
    options.include50           = false;
    
    options.plotFlag            = true;
    options.printPlot           = false;
    options.figureHandle      	= 401;
    
    % Return just the default options struct if no input
    if nargin == 0
        Data           = options;
        return
    end
end

include50       = options.include50;
plotFlag        = options.plotFlag;
printPlot       = options.printPlot;
figureHandle    = options.figureHandle;

usePreSSD = true;
useCorrectOrAll = 'all';
% ***********************************************************************
% Inhibition Function:
%       &
% SSD vs. Proportion of Response trials
%       &
% SSD vs. Proportion(Correct Choice)
% ***********************************************************************

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
ssdArray = ExtraVar.ssdArray;
pSignalArray = ExtraVar.pSignalArray;
targAngleArray = ExtraVar.targAngleArray;
nTarg = length(targAngleArray);
nTrial = size(trialData.trialOutcome, 1);



% Truncate RTs
MIN_RT = 120;
MAX_RT = 1200;
nSTD   = 3;
[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTD);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);
allRT(outlierTrial) = [];




ssdArrayRaw = trialData.stopSignalOn - trialData.responseCueOn;

% If there weren't stop trials, skip all stop-related analyses
if isempty(ssdArray)
    Data = [];
    disp('ccm_inhibition.m: No stop trials or stop trial analyses not requested');
    return
end


if ~include50
    pSignalArray(pSignalArray == .5) = [];
end

if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Constnats, Conditions setup
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Which Signal Strength levels to analyze
switch options.collapseSignal
    case true
        nSignal = 2;
    case false
        nSignal = length(pSignalArray);
end


% If collapsing into all left and all right need to note here that there are "2" angles to deal with
% (important for calling ccm_trial_selection.m)
leftTargInd = (targAngleArray < -89) & (targAngleArray > -270) | ...
    (targAngleArray > 90) & (targAngleArray < 269);
rightTargInd = ~leftTargInd;
if options.collapseTarg
    nTargPair = 1;
else
    nTargPair = sum(rightTargInd);
    % do nothing, all target angles will be considered separately
end















% Loop through all right targets (or collapse them if desired) and
% account for all target pairs if the session had more than one target
% pair
for kTarg = 1 : nTargPair
    
    
    
    ssd                 = cell(nSignal, 1);
    stopStopTrial       = cell(nSignal, length(ssdArray));
    stopTargRespondTrial       = cell(nSignal, length(ssdArray));
    stopDistRespondTrial       = cell(nSignal, length(ssdArray));
    goTrialTotal        = cell(nSignal, 1);
    goTotalRT           = cell(nSignal, 1);
    stopRespondProb     = nan(nSignal, length(ssdArray));
    stopTargetProb      = nan(nSignal, length(ssdArray));
    inhibitionFn        = cell(nSignal, 1);
    inhFnGoMinusSSD     = cell(nSignal, 1);
    goRTMinusSSD        = cell(nSignal, 1);
    goRTMean            = nan(nSignal, 1);
    nStop               = nan(nSignal, length(ssdArray));
    nStopStop           = nan(nSignal, length(ssdArray));
    nStopTarg           = nan(nSignal, length(ssdArray));
    nStopDist           = nan(nSignal, length(ssdArray));
    stopAccuracy     	= nan(nSignal, length(ssdArray));
    ssrtGrandCor           = nan(nSignal, 1);
    ssrtMeanCor            = nan(nSignal, 1);
    ssrtIntegrationCor     = cell(nSignal, 1);
    ssrtIntegrationWeightedCor = nan(nSignal, 1);
    ssrtIntegrationSimpleCor = nan(nSignal, 1);
    ssrtGrandErr           = nan(nSignal, 1);
    ssrtMeanErr            = nan(nSignal, 1);
    ssrtIntegrationErr     = cell(nSignal, 1);
    ssrtIntegrationWeightedErr = nan(nSignal, 1);
    ssrtIntegrationSimpleErr = nan(nSignal, 1);
    stopTargRT          = cell(nSignal, length(ssdArray));
    stopDistRT          = cell(nSignal, length(ssdArray));
    stopRespondRT       = nan(nSignal, length(ssdArray));
    stopRespondRTPredict = nan(nSignal, length(ssdArray));
    
    conditionSSD        = cell(nSignal, 1);
    
    
    % Get default trial selection options
    optSelect = ccm_trial_selection;
    
    
    if plotFlag
        figureHandle = figureHandle + 1;
        minColorGun = .25;
        maxColorGun = 1;
        nRow = 3;
        nColumn = 3;
        screenOrSave = 'save';
        if printPlot
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
        else
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
        end
        choicePlotXMargin = .03;
        ssdMargin = 20;
        ylimArray = [];
        
        % axes names
        axInhGrand = 1;
        axInhEachCorr = 2;
        axInhEachErr = 7;
        axInhRTSSD = 6;
        SSDvPCorrect = 3;
        ssrtPRight = 4;
        SSDvSigStrength = 5;
        
        
        % Set up plot axes
        % inhibition function: grand (collapsed over all signal strengths)
        ax(axInhGrand) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
        cla
        hold(ax(axInhGrand), 'on')
        
        % inhibition function for each signal strength
        ax(axInhEachCorr) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
        cla
        hold(ax(axInhEachCorr), 'on')

        % inhibition function for each signal strength
        ax(axInhEachErr) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
        cla
        hold(ax(axInhEachErr), 'on')

        % goRTMean - SSD inhibition function for each signal strength
        ax(axInhRTSSD) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 3) yAxesPosition(1, 3) axisWidth axisHeight]);
        cla
        hold(ax(axInhRTSSD), 'on')
        
        % p(right) vs ssrt
        ax(ssrtPRight) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 3) yAxesPosition(2, 3) axisWidth axisHeight]);
        cla
        hold(ax(ssrtPRight), 'on')
        
        % SSD vs p(correct)
        ax(SSDvPCorrect) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 1) yAxesPosition(3, 1) axisWidth axisHeight]);
        cla
        hold(ax(SSDvPCorrect), 'on')
        
        % SSD vs p(correct)
        ax(SSDvSigStrength) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 3) yAxesPosition(3, 3) axisWidth axisHeight]);
        cla
        hold(ax(SSDvSigStrength), 'on')
    end
    
    
    
    
    
    
    
    
    
    for iPropIndex = 1 : nSignal
        % Which Signal Strength levels to analyze
        switch options.collapseSignal
            case true
                if iPropIndex == 1
                    iPct = pSignalArray(pSignalArray < .5) * 100;
                elseif iPropIndex == 2
                    iPct = pSignalArray(pSignalArray > .5) * 100;
                end
            case false
                iPct = pSignalArray(iPropIndex) * 100;
        end
        optSelect.rightCheckerPct = iPct;
        
        
        
        % If collapsing into all left and all right or all up/all down,
        % need to note here that there are "2" angles to deal with
        % (important for calling ccm_trial_selection.m)
        if options.collapseTarg && iPct(1) > 50
            kAngle = 'right';
        elseif options.collapseTarg && iPct(1) < 50
            kAngle = 'left';
        else
            if iPct(1) > 50
                rightTargArray = targAngleArray(rightTargInd);
                kAngle = rightTargArray(kTarg);
            elseif iPct(1) < 50
                leftTargArray = targAngleArray(leftTargInd);
                kAngle = leftTargArray(kTarg);
            end
        end
        optSelect.targDir = kAngle;
        
        
        
        % Get go RT distribution, to be used below for calculating SSRT and for
        % predicting stop Incorrect RTs (to compare with observed, as per the
        % race model).
        optSelect.ssd = 'none';
                optSelect.outcome     = {'goCorrectTarget','targetHoldAbort'};
        goTargTrial = ccm_trial_selection(trialData, optSelect);
        iGoTargInd = zeros(nTrial, 1);
        iGoTargInd(goTargTrial) = 1;

        optSelect.outcome     = {'goCorrectDistractor','distractorHoldAbort'};
        goDistrial = ccm_trial_selection(trialData, optSelect);
        iGoDistInd = zeros(nTrial, 1);
        iGoDistInd(goDistrial) = 1;
       
        
        
        if sum(iGoTargInd)
            iGoTargTrial = find(iGoTargInd);
            goTrialTotal{iPropIndex} = iGoTargTrial;  % Keep track of totals for grand inhibition fnct
            
            iGoTargTrialRT = allRT(iGoTargTrial);
            goTotalRT{iPropIndex} = iGoTargTrialRT;
            goRTMean(iPropIndex) = nanmean(iGoTargTrialRT);
        end
        iGoTargRTCum = 1/length(iGoTargTrial):1/length(iGoTargTrial):1; %y-axis of a cumulative prob dist
        iGoTargRTSort = sort(iGoTargTrialRT);

                
        if sum(iGoDistInd)
            iGoDistTrial = find(iGoDistInd);
            goTrialTotal{iPropIndex} = iGoDistTrial;  % Keep track of totals for grand inhibition fnct
            
            iGoDistTrialRT = allRT(iGoDistTrial);
            goTotalRT{iPropIndex} = iGoDistTrialRT;
            goRTMean(iPropIndex) = nanmean(iGoDistTrialRT);
        end
        iGoDistRTCum = 1/length(iGoDistTrial):1/length(iGoDistTrial):1; %y-axis of a cumulative prob dist
        iGoDistRTSort = sort(iGoDistTrialRT);

        
        % Get stopping data within each SSD
        for jSSDIndex = 1 : length(ssdArray)
            %             tic
            jSSD = ssdArray(jSSDIndex);
            optSelect.ssd       = jSSD;
            
            % stop correct trials
            optSelect.outcome       = {'stopCorrect'};
            stopStopTrial = ccm_trial_selection(trialData, optSelect);
            stopStopTrial{iPropIndex, jSSDIndex} = stopStopTrial;  % Keep track of totals for grand inhibition fnct
            nStopStop(iPropIndex, jSSDIndex) = length(stopStopTrial);
            
            % stop incorrect target trials
            optSelect.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
            %             optSelect.outcome       = {'stopIncorrectTarget', 'targetHoldAbort'};
            stopTargTrial = ccm_trial_selection(trialData, optSelect);
            if ~usePreSSD
                preSSDTrial = trialData.rt(stopTargTrial) < trialData.ssd(stopTargTrial);
                stopTargTrial(preSSDTrial) = [];
            end
            nStopTarg(iPropIndex, jSSDIndex) = length(stopTargTrial);
            
            % stop incorrect distractor trials
            optSelect.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
            %             optSelect.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort'};
            stopDistTrial = ccm_trial_selection(trialData, optSelect);
            if ~usePreSSD
                preSSDTrial = trialData.rt(stopDistTrial) < trialData.ssd(stopDistTrial);
                stopDistTrial(preSSDTrial) = [];
            end
            nStopDist(iPropIndex, jSSDIndex) = length(stopDistTrial);
            
            
                                stopIncorrectTrial = union(stopTargTrial, stopDistTrial);

            if nStopTarg(iPropIndex, jSSDIndex) > 0 || nStopDist(iPropIndex, jSSDIndex) > 0
            stopAccuracy(iPropIndex, jSSDIndex) = length(stopTargTrial) / (length(stopTargTrial) + length(stopDistTrial));
            elseif ~isempty(stopStopTrial) && isempty(stopIncorrectTrial)
            stopAccuracy(iPropIndex, jSSDIndex) = .5;  % Fake number 
            end
            
            stopTargRT{iPropIndex, jSSDIndex} = trialData.responseOnset(stopTargTrial) - trialData.responseCueOn(stopTargTrial);
            stopDistRT{iPropIndex, jSSDIndex} = trialData.responseOnset(stopDistTrial) - trialData.responseCueOn(stopDistTrial);
            

            stopTargRespondTrial{iPropIndex, jSSDIndex} = stopTargTrial;  % Keep track of totals for grand inhibition fnct
            stopDistRespondTrial{iPropIndex, jSSDIndex} = stopDistTrial;  % Keep track of totals for grand inhibition fnct
            
            % Inhibition function data points:
            stopRespondProb(iPropIndex, jSSDIndex) = length(stopIncorrectTrial) / (length(stopStopTrial) + length(stopIncorrectTrial));
            nStop(iPropIndex, jSSDIndex) = length(stopStopTrial) + length(stopIncorrectTrial);
            
            % p(Correct choice) vs. SSD data points:
            stopTargetProb(iPropIndex, jSSDIndex) = length(stopTargTrial) / (length(stopTargTrial) + length(stopDistTrial));
            
            
            % Predicted mean stopIncorrectRT (based on goRTs and p(noncanceled):
            jRTIndex = find(iGoTargRTCum >= stopRespondProb(iPropIndex, jSSDIndex), 1);   %match estimated p(noncan|SSD) to p(RT)
            stopRespondRTPredict(iPropIndex, jSSDIndex) = nanmean(iGoTargRTSort(1 : jRTIndex-1));
            % Observed mean stop IncorrectRT:
            if length(allRT(stopIncorrectTrial)) > 1
                stopRespondRT(iPropIndex, jSSDIndex) = nanmean(allRT(stopIncorrectTrial));
            end
        end % iSSDIndex
        
        
        % Inhibition function data points
        iStopProbRespond = stopRespondProb(iPropIndex, :);
        keepSSD = ~isnan(iStopProbRespond);
        iStopProbRespond = iStopProbRespond(keepSSD);
        
        % p(Correct choice) vs. SSD data points:
        iStopProbTarget = stopTargetProb(iPropIndex, :);
        keepSSDChoice = ~isnan(iStopProbTarget);
        iStopProbTarget = iStopProbTarget(keepSSDChoice);
        iSSDArrayChoice = ssdArray(keepSSDChoice)';
        [p, s] = polyfit(iSSDArrayChoice, iStopProbTarget, 1);
        xVal = min(iSSDArrayChoice) : max(iSSDArrayChoice);
        
        
        % Inhibition function calculation
        iNStop = nStop(iPropIndex, :);
        iNStop = iNStop(keepSSD);
        iSSDArray = ssdArray(keepSSD);
        ssd{iPropIndex} = iSSDArray;
        
        [fitParameters, lowestSSE] = Weibull(iSSDArray, iStopProbRespond, iNStop);
        ssdTimePoints = ssdArray(1) : ssdArray(end);
        inhibitionFn{iPropIndex} = weibull_curve(fitParameters, ssdTimePoints);
        
        
        
        
        % SSRT: get go RTs and number of stop trials (already have other
        % necessary variables)

        % Correct Choices:
        iNStopT = iNStop .* stopAccuracy(iPropIndex, keepSSD);
%         iNStopT(iNStopT == 0) = [];
        ssrt                            = get_ssrt(iSSDArray, iStopProbRespond, iNStopT, iGoTargTrialRT, fitParameters);
        ssrtGrandCor(iPropIndex)        = ssrt.grand;
        ssrtIntegrationWeightedCor(iPropIndex) = ssrt.integrationWeighted;
        ssrtIntegrationSimpleCor(iPropIndex) = ssrt.integrationSimple;
        ssrtIntegrationCor{iPropIndex}  = ssrt.integration;
        ssrtMeanCor(iPropIndex)         = ssrt.mean;
        
        % Error Choices:
        iNStopD = iNStop .* (1 - stopAccuracy(iPropIndex, keepSSD));
%         iNStopD(iNStopD == 0) = [];       
        ssrt                            = get_ssrt(iSSDArray, iStopProbRespond, iNStopD, iGoDistTrialRT, fitParameters);
        ssrtGrandErr(iPropIndex)        = ssrt.grand;
        ssrtIntegrationWeightedErr(iPropIndex) = ssrt.integrationWeighted;
        ssrtIntegrationSimpleErr(iPropIndex) = ssrt.integrationSimple;
        ssrtIntegrationErr{iPropIndex}  = ssrt.integration;
        ssrtMeanErr(iPropIndex)         = ssrt.mean;
        
        
        
        
        
        % Mean SSD in each signal strength
        optSelect.outcome = {'targetHoldAbort', 'stopIncorrectTarget', 'stopIncorrectPreSSDTarget', ...
            'distractorHoldAbort', 'stopIncorrectDistractor', 'stopIncorrectPreSSDDistractor'};
        optSelect.ssd = 'collapse';
        allStopTrial = ccm_trial_selection(trialData, optSelect);
        
        conditionSSD{iPropIndex} = ssdArrayRaw(allStopTrial);
        
        
        % Inhibition functions using goRT - SSD in each signal strength:
        goRTMinusSSD{iPropIndex} = goRTMean(iPropIndex) - iSSDArray;
        offsetVal = 1000;
        [fitParameters, lowestSSE] = Weibull(-goRTMinusSSD{iPropIndex} + offsetVal, iStopProbRespond', iNStop);
        goSSDTimepoints = offsetVal + (min(-goRTMinusSSD{iPropIndex}) : max(-goRTMinusSSD{iPropIndex}));
        inhFnGoMinusSSD{iPropIndex} = weibull_curve(fitParameters, goSSDTimepoints);
        
        
        
        
        
        
        
        if plotFlag
            % Determine color to use for plot based on which checkerboard color
            % proportion being used.
            if options.collapseSignal
                cMap = ccm_colormap([0 1]);
            else
                cMap = ccm_colormap(pSignalArray);
            end
            inhColor = cMap(iPropIndex,:);
            
            
%             plot(ax(axInhEach), ssdTimePoints, inhibitionFn{iPropIndex}, 'color', inhColor, 'linewidth', 2)
            plot(ax(axInhEachCorr), ssd{iPropIndex}, iStopProbRespond, 'color', inhColor, 'linewidth', 2)
            %         if iPropIndex > 1 && iPropIndex < 4
            plot(ax(axInhRTSSD), (goSSDTimepoints - offsetVal), inhFnGoMinusSSD{iPropIndex}, '-', 'color', inhColor, 'linewidth', 2);%'linewidth', 2)
            plot(ax(axInhRTSSD), -goRTMinusSSD{iPropIndex}, iStopProbRespond, '.', 'color', inhColor, 'markersize', 15);%'linewidth', 2)
            %         end
            plot(ax(ssrtPRight), pSignalArray(iPropIndex), ssrtGrandCor(iPropIndex), '.', 'markersize', 30, 'color', inhColor)
            plot(ax(ssrtPRight), pSignalArray(iPropIndex), ssrtGrandErr(iPropIndex), '.', 'markersize', 30, 'color', 'r')
            plot(ax(SSDvPCorrect), xVal, p(1) * xVal + p(2), 'color', inhColor, 'linewidth', 2)
            plot(ax(SSDvSigStrength), pSignalArray(iPropIndex), mean(conditionSSD{iPropIndex}), '.',  'color', inhColor, 'markersize', 30)
        end
        
        
        
        
        
    end % iPropIndex
    
    
    
    if ~options.collapseSignal
        flipPropArray = pSignalArray;
        flipPropArray(flipPropArray > .5) = abs(1 - flipPropArray(flipPropArray > .5));
        [p, s] = polyfit(flipPropArray, ssrtGrandCor, 1);
        [y, delta] = polyval(p, flipPropArray, s);
        stats = regstats(flipPropArray, ssrtGrandCor);
        fprintf('p-value for regression: %.4f\n', stats.tstat.pval(2))
        R = corrcoef(flipPropArray, ssrtGrandCor);
        Rsqrd = R(1, 2)^2;
        cov(flipPropArray, ssrtGrandCor);
        xVal = min(flipPropArray) : .001 : .5;
        yVal = p(1) * xVal + p(2);
        regressionLineY = [yVal(1:end-1), fliplr(yVal)];
        regressionLineX = min(flipPropArray) : .001 : abs(1 - min(flipPropArray));
    end
    
    
    if plotFlag
        % inhibition function
        %     xlim(ax(axInhEach), [ssdTimePoints(1) 400])
        xlim(ax(axInhEachCorr), [ssdTimePoints(1) ssdTimePoints(end)])
        %     set(ax(axInhEach), 'xtick', ssdArray)
        %     set(ax(axInhEach), 'xtickLabel', ssdArray)
        ylim(ax(axInhEachCorr), [0 1]);
        set(get(ax(axInhEachCorr), 'ylabel'), 'String', 'p(Respond | stop)')
        
        
        xlim(ax(axInhRTSSD), [goSSDTimepoints(1) - offsetVal goSSDTimepoints(1) - offsetVal+400])
        
        %         ssrt vs signal strength
        set(ax(ssrtPRight), 'Ylim', [0 round(2*mean(ssrtGrandCor))])
        set(ax(ssrtPRight), 'Ylim', [0 300])
        set(get(ax(ssrtPRight), 'ylabel'), 'String', 'SSRT')
        plot(ax(ssrtPRight), [.5 .5], ylim, '--k')
        if ~options.collapseSignal
            xlim(ax(ssrtPRight), [pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
            set(ax(ssrtPRight), 'xtick', pSignalArray)
            set(ax(ssrtPRight), 'xtickLabel', pSignalArray*100)
            plot(ax(ssrtPRight), regressionLineX, regressionLineY, 'b')
        else
            xlim(ax(ssrtPRight), [pSignalArray(1) - choicePlotXMargin pSignalArray(2) + choicePlotXMargin])
            set(ax(ssrtPRight), 'xtick', [pSignalArray(1) pSignalArray(2)])
            set(ax(ssrtPRight), 'xtickLabel', {'left','right'})
        end
        
        
        % ssd vs p(correct)
        xlim(ax(SSDvPCorrect), [ssdArray(1)-ssdMargin ssdArray(end)+ssdMargin])
        set(ax(SSDvPCorrect), 'xtick', ssdArray)
        set(ax(SSDvPCorrect), 'xtickLabel', ssdArray)
        ylim(ax(SSDvPCorrect), [-.05 1.05]);
        set(get(ax(SSDvPCorrect), 'ylabel'), 'String', 'p(Correct)')
        
        if options.collapseSignal
            xlim(ax(SSDvSigStrength), [pSignalArray(1) - choicePlotXMargin pSignalArray(2) + choicePlotXMargin])
            set(ax(SSDvSigStrength), 'xtick', [pSignalArray(1) pSignalArray(2)])
            set(ax(SSDvSigStrength), 'xticklabel', {'left','right'})
        else
            xlim(ax(SSDvSigStrength), [pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
            set(ax(SSDvSigStrength), 'xtick', pSignalArray)
        end
        set(get(ax(SSDvSigStrength), 'ylabel'), 'String', 'Mean SSD')
        
        
        
        figure(64)
        clf
        %      surf(ax(axPred), stopRespondRT)
        %      surf(ax(axPred), stopRespondRTPredict)
        colormap(hsv)
        %      surf(ssdArray, signalStrength, stopRespondRT)
        surface(stopRespondRT)
        hold on
    end % if plotflag
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % Also calculate a grand inhibition function, SSRT, and predicted vs.
    % observed non-canceled RTs across discriminability levels
    goGrandRT = [];
    for iPropIndex = 1 : nSignal
        goGrandRT = [goGrandRT; goTotalRT{iPropIndex}];
    end
    goGrandRT = sort(goGrandRT);
    rtCum = 1/length(goGrandRT):1/length(goGrandRT):1; %y-axis of a cumulative prob dist
    
    % Grand Inhibition function
    nStopTargGrand = zeros(length(ssdArray), 1);
    nStopDistGrand = zeros(length(ssdArray), 1);
    nStopCorrectGrand = zeros(length(ssdArray), 1);
    stopRespondProbGrand = nan(length(ssdArray), 1);
    nStopGrand = zeros(length(ssdArray), 1);
    for iSSDIndex = 1 : length(ssdArray)
        for jPropIndex = 1 : nSignal;
            nStopTargGrand(iSSDIndex) = nStopTargGrand(iSSDIndex) + length(stopTargRespondTrial{jPropIndex, iSSDIndex});
            nStopDistGrand(iSSDIndex) = nStopDistGrand(iSSDIndex) + length(stopDistRespondTrial{jPropIndex, iSSDIndex});
            nStopCorrectGrand(iSSDIndex) = nStopCorrectGrand(iSSDIndex) + length(stopStopTrial{jPropIndex, iSSDIndex});
        end
        nStopGrand(iSSDIndex) = nStopTargGrand(iSSDIndex) + nStopCorrectGrand(iSSDIndex);
        % Inhibition function data points
        stopRespondProbGrand(iSSDIndex) = nStopTargGrand(iSSDIndex) / nStopGrand(iSSDIndex);
        
        
        % Predicted mean stopIncorrectRT (based on goRTs and p(noncanceled):
        iRTIndex = find(rtCum >= stopRespondProbGrand(iSSDIndex), 1);   %match estimated p(noncan|SSD) to p(RT)
        stopRespondRTPredictGrand(iSSDIndex) = nanmean(goGrandRT(1 : iRTIndex-1));
        % Observed mean stop IncorrectRT:
        emptyData = cellfun(@isempty, stopTargRespondTrial(:, iSSDIndex));
        stopInocrrectTrial = cell2mat(stopTargRespondTrial(~emptyData, iSSDIndex));
        %         stopInocrrectTrial = cell2mat(cellfun(@(x) x', stopRespondTrial(:, iSSDIndex), 'uniformoutput', false));
        stopTargRespondRTObserveGrand(iSSDIndex) = nanmean(allRT(stopInocrrectTrial));
        
        emptyData = cellfun(@isempty, stopDistRespondTrial(:, iSSDIndex));
        stopInocrrectTrial = cell2mat(stopDistRespondTrial(~emptyData, iSSDIndex));
        %         stopInocrrectTrial = cell2mat(cellfun(@(x) x', stopRespondTrial(:, iSSDIndex), 'uniformoutput', false));
        stopDistRespondRTObserveGrand(iSSDIndex) = nanmean(allRT(stopInocrrectTrial));
        
    end
    keepSSD = ~isnan(stopRespondProbGrand);
    stopRespondProbGrand = stopRespondProbGrand(keepSSD);
    nStopGrand = nStopGrand(keepSSD);
    ssdArray = ssdArray(keepSSD);
    %     [fitParameters, lowestSSE] = Weibull(ssdArray, stopRespondProbGrand);
    [fitParameters, lowestSSE] = Weibull(ssdArray, stopRespondProbGrand, nStopGrand);
    %     [fitParameters, lowestSSE] = Weibull_fast(ssdArray, stopRespondProbGrand);
    ssdTimePoints = ssdArray(1) : ssdArray(end);
    inhibitionFnGrand = weibull_curve(fitParameters, ssdTimePoints);
    
    if plotFlag
        plot(ax(axInhGrand), ssdTimePoints, inhibitionFnGrand, 'color', 'g', 'linewidth', 2)
        plot(ax(axInhGrand), ssdArray, stopRespondProbGrand, '.k', 'markersize', 25)
        set(ax(axInhGrand),'YLim',[0 1])
        set(get(ax(axInhGrand), 'ylabel'), 'String', 'p(Respond | stop)')
        %     set(ax(axInhGrand),'XLim',[0 800])
        set(ax(axInhGrand),'XLim',[ssdTimePoints(1) ssdTimePoints(end)])
        %     set(ax(axInhGrand),'XLim',[ssdTimePoints(1) 400])
    end
    
    % Grand SSRT
    % [SSRT both_SSRTs ssrtI ssrtM meanSSD] = get_ssrt(ssdArray, stopRespondProbGrand, nStopGrand, goGrandRT, fitParameters);
    % ssrtGrand = SSRT;
    % ssrtMeanGrand = ssrtM;
    % ssrtIntegrationGrand = ssrtI;
    
    ssrtCollapse        = get_ssrt(ssdArray, stopRespondProbGrand, nStopGrand, goGrandRT, fitParameters);
    ssrtCollapseGrand   = ssrtCollapse.grand;
    ssrtCollapseIntegrationWeighted = ssrtCollapse.integrationWeighted;
    ssrtCollapseIntegration = ssrtCollapse.integration;
    ssrtCollapseMean    = ssrtCollapse.mean;
    
    if plotFlag && ~options.collapseSignal
        plot(ax(ssrtPRight), xlim, [ssrtCollapseGrand ssrtCollapseGrand], '--k')
        
    end
    
    
    
    %         if printPlot
    %             localFigurePath = local_figure_path;
    %             print(figureHandle,[localFigurePath, sessionID, '_ccm_inhibition'],'-dpdf', '-r300')
    %         end
    
    
    
    
    
    Data(kTarg).ssrtGrand               = ssrtGrandCor;
    Data(kTarg).ssrtMean               = ssrtMeanCor;
    Data(kTarg).ssrtIntegration               = ssrtIntegrationCor;
    Data(kTarg).ssrtIntegrationWeighted               = ssrtIntegrationWeightedCor;
    Data(kTarg).ssrtIntegrationSimple               = ssrtIntegrationSimpleCor;
    Data(kTarg).ssrtCollapseGrand     	= ssrtCollapseGrand;
    Data(kTarg).ssrtCollapseIntegrationWeighted     	= ssrtCollapseIntegrationWeighted;
    Data(kTarg).ssrtCollapseIntegration     	= ssrtCollapseIntegration;
    Data(kTarg).ssrtCollapseMean     	= ssrtCollapseMean;
    Data(kTarg).ssd                = ssd;
    Data(kTarg).pSignalArray       = pSignalArray;
    Data(kTarg).goRTMinusSSD       = goRTMinusSSD;
    Data(kTarg).ssdArray           = ssdArray;
    Data(kTarg).stopRespondProb    = stopRespondProb;
    Data(kTarg).nStop              = nStop;
    Data(kTarg).nStopStop          = nStopStop;
    Data(kTarg).nStopTarg          = nStopTarg;
    Data(kTarg).nStopDist          = nStopDist;
    Data(kTarg).stopTargetProb     = stopTargetProb;
    Data(kTarg).inhibitionFn       = inhibitionFn;
    % Data(kTarg).ssrt               = ssrt;
    % Data(kTarg).ssrtMean               = ssrtMean;
    % Data(kTarg).ssrtIntegration               = ssrtIntegration;
    Data(kTarg).stopRespondProbGrand = stopRespondProbGrand;
    Data(kTarg).inhibitionFnGrand  = inhibitionFnGrand;
    % Data(kTarg).ssrtGrand          = ssrtGrand;
    % Data(kTarg).ssrtMeanGrand          = ssrtMeanGrand;
    % Data(kTarg).ssrtIntegrationGrand          = ssrtIntegrationGrand;
    
    Data(kTarg).stopRespondRTPredictGrand = stopRespondRTPredictGrand;
    Data(kTarg).stopTargRespondRTObserveGrand = stopTargRespondRTObserveGrand;
    Data(kTarg).stopDistRespondRTObserveGrand = stopDistRespondRTObserveGrand;
    Data(kTarg).stopRespondRTPredict = stopRespondRTPredict;
    Data(kTarg).stopRespondRT = stopRespondRT;
    
end % for kTarg = 1 : nTargPair

end % function













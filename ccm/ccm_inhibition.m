function Data = ccm_inhibition(subjectID, sessionID, options)

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

%%
% sessionID = 'bp093n02';
% Set default options or return a default options structure
if nargin < 3
    options.collapseSignal   	= false;
    options.collapseTarg        = true;
    options.include50           = false;
    options.USE_PRE_SSD         = true;
    options.USE_TWO_COLORS         = false;
    options.INCLUDE_GO_OMISSION = true;
    
    options.plotFlag            = true;
    options.printPlot           = false;
    options.figureHandle      	= 403;
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
INCLUDE_GO_OMISSION = options.INCLUDE_GO_OMISSION;


useCorrectOrAll = 'correct';
plotSurface = false;


% ***********************************************************************
% Inhibition Function:
%       &
% SSD vs. Proportion of Response trials
%       &
% SSD vs. Proportion(Correct Choice)
% ***********************************************************************

% Load the behavior-only version of the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID, ccm_min_vars);
ssdArray = ExtraVar.ssdArray;
pSignalArray = unique(trialData.targ1CheckerProp);
if options.USE_TWO_COLORS
    if length(pSignalArray) == 6
        pSignalArray([2 5]) = [];
    elseif length(pSignalArray) == 7
        pSignalArray([2 4 6]) = [];
    end
end
targAngleArray = unique(trialData.targAngle);
nTrial = size(trialData.trialOutcome, 1);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Constnats, Conditions setup
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Truncate RTs
MIN_RT = 120;
MAX_RT = 1200;
nSTDTruncate   = 4;

% To include Go Omissions for SSRT calculations, set all goIncorrect RTs to
% lowest goIndorrect RT
goIncorrectRT = min(MAX_RT-1, nanmedian(trialData.rt(strcmp(trialData.trialOutcome, 'goIncorrect'))));
trialData.rt(strcmp(trialData.trialOutcome, 'goIncorrect')) = goIncorrectRT;

[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTDTruncate);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);
allRT(outlierTrial) = [];




ssdArrayRaw = trialData.stopSignalOn - trialData.responseCueOn;

% If there weren't stop trials, skip all stop-related analyses
if isempty(ssdArray)
    Data = [];
    %     disp('ccm_inhibition.m: No stop trials or stop trial analyses not requested');
    return
end


if ~include50
    pSignalArray(pSignalArray == .5) = [];
end

if include50
    middleDiff = ceil(length(pSignalArray)/1);
else
    middleDiff = [];
end

choiceLevels = [1:floor(length(pSignalArray)/2) middleDiff floor(length(pSignalArray)/2):-1:1];


if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end



% Which Signal Strength levels to analyze
switch options.collapseSignal
    case true
        nSignal = 2;
    case false
        nSignal = length(pSignalArray);
end
nSSD = length(ssdArray);

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

goColor = [0 0 0];
stopColor = [1 0 0];














% Loop through all right targets (or collapse them if desired) and
% account for all target pairs if the session had more than one target
% pair
for kTarg = 1 : nTargPair
    
    
    
    ssd                 = cell(nSignal, 1);
    stopStopTrial       = cell(nSignal, nSSD);
    stopRespondTrial  	= cell(nSignal, nSSD);
    stopTargTrial       = cell(nSignal, nSSD);
    stopDistTrial       = cell(nSignal, nSSD);
    
    goTotalRT           = cell(nSignal, nSSD);
    goTargRT           = cell(nSignal, nSSD);
    goDistRT           = cell(nSignal, nSSD);
    goRTMean            = nan(nSignal, nSSD);
    
    stopTargRT          = cell(nSignal, nSSD);
    stopDistRT          = cell(nSignal, nSSD);
    stopRespondRT       = nan(nSignal, nSSD);
    stopRespondRTPredict = nan(nSignal, nSSD);
    stopRespondProb     = nan(nSignal, nSSD);
    stopTargetProb      = nan(nSignal, nSSD);
    inhibitionFn        = cell(nSignal, 1);
    inhFnGoMinusSSD     = cell(nSignal, 1);
    goRTMinusSSD        = cell(nSignal, 1);
    nStop               = nan(nSignal, nSSD);
    nStopStop           = nan(nSignal, nSSD);
    nStopTarg           = nan(nSignal, nSSD);
    nStopDist           = nan(nSignal, nSSD);
    ssrtGrand           = nan(nSignal, 1);
    ssrtMean            = nan(nSignal, 1);
    ssrtIntegration     = cell(nSignal, 1);
    ssrtIntegrationWeighted = nan(nSignal, 1);
    ssrtIntegrationSimple = nan(nSignal, 1);
    
    conditionSSD        = cell(nSignal, 1);
    weibullParams        = cell(nSignal, 1);
    
    
    % Get default ccm options
    optSelect = ccm_options;
    
    
    if plotFlag
        %         figureHandle = figureHandle + 1;
        minColorGun = .25;
        maxColorGun = 1;
        nRow = 3;
        nColumn = 4;
        screenOrSave = 'save';
        if printPlot
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
        else
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
        end
        clf
        choicePlotXMargin = .03;
        ssdMargin = 20;
        ylimArray = [];
        
        % axes names
        axInhGrand = 1;
        axInhEach = 2;
        SSDvPCorrect = 3;
        ssrtPRight = 4;
        SSDvSigStrength = 5;
        axInhRTSSD = 6;
        nTrialPerCoh = 7;
        axWeibullParam = 8;
        
        
        % Set up plot axes
        % inhibition function: grand (collapsed over all signal strengths)
        ax(axInhGrand) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
        cla
        hold(ax(axInhGrand), 'on')
        
        % inhibition function for each signal strength
        ax(axInhEach) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
        cla
        hold(ax(axInhEach), 'on')
        
        % goRTMean - SSD inhibition function for each signal strength
        ax(axInhRTSSD) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 3) yAxesPosition(1, 3) axisWidth axisHeight]);
        cla
        hold(ax(axInhRTSSD), 'on')
        
        % p(right) vs ssrt
        ax(ssrtPRight) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 3) yAxesPosition(2, 3) axisWidth axisHeight]);
        cla
        hold(ax(ssrtPRight), 'on')
        
        %         % SSD vs p(correct)
        %         ax(SSDvPCorrect) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 1) yAxesPosition(3, 1) axisWidth axisHeight]);
        %         cla
        %         hold(ax(SSDvPCorrect), 'on')
        
        % SSD vs p(correct)
        ax(nTrialPerCoh) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 2) yAxesPosition(3, 2) axisWidth axisHeight]);
        cla
        hold(ax(nTrialPerCoh), 'on')
        
        % SSD vs p(correct)
        ax(SSDvSigStrength) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 3) yAxesPosition(3, 3) axisWidth axisHeight]);
        cla
        hold(ax(SSDvSigStrength), 'on')
        
        % Weibull params w.r.t. Color Coherence
        ax(axWeibullParam) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 4) yAxesPosition(1, 4) axisWidth axisHeight]);
        cla
        hold(ax(axWeibullParam), 'on')
        
        
        
        
        
        
        
        
    end
    
    
    
    
    
    
    
    % Get stopping data within each SSD
    for jSSDInd = 1 : nSSD
        
        
        
        
        for iPropInd = 1 : nSignal
            % Which Signal Strength levels to analyze
            switch options.collapseSignal
                case true
                    if iPropInd == 1
                        iPct = pSignalArray(pSignalArray < .5) * 100;
                    elseif iPropInd == 2
                        iPct = pSignalArray(pSignalArray > .5) * 100;
                    end
                case false
                    iPct = pSignalArray(iPropInd) * 100;
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
            
            
            
            
            
            
            %    GATHER STOP TRIALS FROM DIFFERENT OUTCOMES
            %  ============================================================
            jSSD = ssdArray(jSSDInd);
            optSelect.ssd       = jSSD;
            
            
            % Stop correct trials
            optSelect.outcome       = {'stopCorrect'};
            jStopCorrectTrial = ccm_trial_selection(trialData, optSelect);
            stopStopTrial{iPropInd, jSSDInd} = jStopCorrectTrial;  % Keep track of totals for grand inhibition fnct
            nStopStop(iPropInd, jSSDInd) = length(jStopCorrectTrial);
            
            
            
            % Stop incorrect target trials
            optSelect.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
            %             optSelect.outcome       = {'stopIncorrectTarget', 'targetHoldAbort'};
            jStopTargTrial = ccm_trial_selection(trialData, optSelect);
            
            
            
            % Stop incorrect distractor trials
            optSelect.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
            %             optSelect.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort'};
            jStopDistTrial = ccm_trial_selection(trialData, optSelect);
            
            
            
            
            
            %    GATHER GO TRIALS FROM DIFFERENT OUTCOMES
            %  ============================================================
            
            % Get go RT distribution, to be used below for calculating SSRT and for
            % predicting stop Incorrect RTs (to compare with observed, as per the
            % race model).
            optSelect.ssd = 'none';
            
            % Go trial correct target
            optSelect.outcome           = {'goCorrectTarget', 'targetHoldAbort'};
            jGoTargTrial                 = ccm_trial_selection(trialData, optSelect);
            
            
            % Go trial errors
            optSelect.outcome           = {'goCorrectDistractor', 'distractorHoldAbort'};
            jGoDistTrial                 = ccm_trial_selection(trialData, optSelect);
            
            % Go omission trials
            optSelect.outcome           = {'goIncorrect'};
            jGoOmissionTrial                 = ccm_trial_selection(trialData, optSelect);
            
            
            
            
            % ADJUST VALID (STOP AND GO) TRIALS BASED ON WHETHER RT WAS MADE ON STOP TRIALS BEFORE SSD
            %  ============================================================
            if ~options.USE_PRE_SSD
                stopTargPreSSDTrial = trialData.rt(jStopTargTrial) < jSSD;
                jStopTargTrial(stopTargPreSSDTrial) = [];
                
                goTargPreSSDTrial = trialData.rt(jGoTargTrial) < jSSD;
                jGoTargTrial(goTargPreSSDTrial) = [];
                
                stopDistPreSSDTrial = trialData.rt(jStopDistTrial) < jSSD;
                jStopDistTrial(stopDistPreSSDTrial) = [];
                
                goDistPreSSDTrial = trialData.rt(jGoDistTrial) < jSSD;
                jGoDistTrial(goDistPreSSDTrial) = [];
            end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            jGoTrial = sort([jGoDistTrial; jGoTargTrial]);
            % Include Go Omissions (see consensus stop signal task paper)
            if INCLUDE_GO_OMISSION
            jGoTrial = sort([jGoTrial; jGoOmissionTrial]);
            end
            
            jGoTrialRT = allRT(jGoTrial);
            goTotalRT{iPropInd, jSSDInd} = jGoTrialRT;
            goTargRT{iPropInd, jSSDInd} = allRT(jGoTargTrial);
            goDistRT{iPropInd, jSSDInd} = allRT(jGoDistTrial);
            goRTMean(iPropInd, jSSDInd) = nanmean(jGoTrialRT);
            iRTCum = 1/length(jGoTrial):1/length(jGoTrial):1; %y-axis of a cumulative prob dist
            iRTSort = sort(jGoTrialRT);
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            nStopTarg(iPropInd, jSSDInd) = length(jStopTargTrial);
            nStopDist(iPropInd, jSSDInd) = length(jStopDistTrial);
            
            
            stopTargRT{iPropInd, jSSDInd} = trialData.rt(jStopTargTrial);
            stopDistRT{iPropInd, jSSDInd} = trialData.rt(jStopDistTrial);
            
            % stop incorrect trials for inhibition: do we want stop incorrect to the target or to
            % target/distractor?
            stopIncorrectTrial = union(jStopTargTrial, jStopDistTrial);
            stopIncorrectTrial = stopIncorrectTrial(:);
            stopRespondTrial{iPropInd, jSSDInd} = stopIncorrectTrial;
            
            stopIncorrectTrialTarg = jStopTargTrial;
            if size(stopIncorrectTrialTarg, 2) > 1
                stopIncorrectTrialTarg = stopIncorrectTrialTarg';
            end
            stopTargTrial{iPropInd, jSSDInd} = stopIncorrectTrialTarg;  % Keep track of totals for grand inhibition fnct
            
            stopIncorrectTrialDist = jStopDistTrial;
            if size(stopIncorrectTrialDist, 2) > 1
                stopIncorrectTrialDist = stopIncorrectTrialDist';
            end
            stopDistTrial{iPropInd, jSSDInd} = stopIncorrectTrialDist;  % Keep track of totals for grand inhibition fnct
            
            
            % Inhibition function data points:
            stopRespondProb(iPropInd, jSSDInd) = length(stopIncorrectTrial) / (length(jStopCorrectTrial) + length(stopIncorrectTrial));
            nStop(iPropInd, jSSDInd) = length(jStopCorrectTrial) + length(stopIncorrectTrial);
            
            % p(Correct choice) vs. SSD data points:
            stopTargetProb(iPropInd, jSSDInd) = length(jStopTargTrial) / (length(jStopTargTrial) + length(jStopDistTrial));
            
            
            % Predict noncanceled stop RTs:
            % ----------------------------------------------------------
            % Predicted mean stopIncorrectRT (based on goRTs and p(noncanceled):
            jRTIndex = find(iRTCum >= stopRespondProb(iPropInd, jSSDInd), 1);   %match estimated p(noncan|SSD) to p(RT)
            stopRespondRTPredict(iPropInd, jSSDInd) = nanmean(iRTSort(1 : jRTIndex-1));
            % Observed mean stop IncorrectRT:
            if length(allRT(stopIncorrectTrial)) > 1
                stopRespondRT(iPropInd, jSSDInd) = nanmean(allRT(stopIncorrectTrial));
            end
            
            
            
            
            
        end % iSSDIndex
        
    end % iPropIndex
    
    
    
    
    
    
    
    
    
    
    
    for iPropInd = 1 : nSignal
        % Which Signal Strength levels to analyze
        switch options.collapseSignal
            case true
                if iPropInd == 1
                    iPct = pSignalArray(pSignalArray < .5) * 100;
                elseif iPropInd == 2
                    iPct = pSignalArray(pSignalArray > .5) * 100;
                end
            case false
                iPct = pSignalArray(iPropInd) * 100;
        end
        
        
        
        % Inhibition function data points
        iStopProbRespond    = stopRespondProb(iPropInd, :);
        keepSSD             = ~isnan(iStopProbRespond);
        iStopProbRespond    = iStopProbRespond(keepSSD);
        
        % p(Correct choice) vs. SSD data points:
        iStopProbTarget     = stopTargetProb(iPropInd, :);
        keepSSDChoice       = ~isnan(iStopProbTarget);
        iStopProbTarget     = iStopProbTarget(keepSSDChoice);
        iSSDArrayChoice     = ssdArray(keepSSDChoice)';
        [p, s]              = polyfit(iSSDArrayChoice, iStopProbTarget, 1);
        xVal                = min(iSSDArrayChoice) : max(iSSDArrayChoice);
        
        
        % Inhibition function calculation
        iNStop = nStop(iPropInd, :);
        iNStop = iNStop(keepSSD);
        iSSDArray = ssdArray(keepSSD);
        ssd{iPropInd} = iSSDArray;
        
        [fitParameters, lowestSSE] = Weibull(iSSDArray, iStopProbRespond, iNStop);
        ssdTimePoints = ssdArray(1) : ssdArray(end);
        inhibitionFn{iPropInd} = weibull_curve(fitParameters, ssdTimePoints);
        weibullParams{iPropInd} = fitParameters;
        
        % SSRT: get go RTs and number of stop trials (already have other
        % necessary variables)
        iNStop = nStop(iPropInd, :);
        iNStop(iNStop == 0) = [];
        
        if options.USE_PRE_SSD
            ssrt = get_ssrt(iSSDArray, iStopProbRespond, iNStop, goTotalRT{iPropInd, 1}, fitParameters);
            ssrtGrand(iPropInd) = ssrt.grand;
            ssrtIntegrationWeighted(iPropInd) = ssrt.integrationWeighted;
            ssrtIntegrationSimple(iPropInd) = ssrt.integrationSimple;
            ssrtIntegration{iPropInd} = ssrt.integration;
            ssrtMean(iPropInd) = ssrt.mean;
        else
            ssrt = get_ssrt_post_ssd(iSSDArray, iStopProbRespond, iNStop, goTotalRT(iPropInd, keepSSD));
            ssrtIntegrationWeighted(iPropInd) = ssrt.integrationWeighted;
            ssrtIntegration{iPropInd} = ssrt.integration;
        end
        
        
        % Mean SSD in each signal strength
        optSelect = ccm_options;
        optSelect.outcome = {'targetHoldAbort', 'stopIncorrectTarget', 'stopIncorrectPreSSDTarget', ...
            'distractorHoldAbort', 'stopIncorrectDistractor', 'stopIncorrectPreSSDDistractor'};
        optSelect.rightCheckerPct = iPct;
        optSelect.ssd = 'collapse';
        allStopTrial = ccm_trial_selection(trialData, optSelect);
        
        conditionSSD{iPropInd} = ssdArrayRaw(allStopTrial);
        
        
        % Inhibition functions using goRT - SSD in each signal strength:
        goRTMinusSSD{iPropInd} = goRTMean(iPropInd, keepSSD)' - iSSDArray;
        offsetVal = 1000;
        [fitParameters, lowestSSE] = Weibull(-goRTMinusSSD{iPropInd} + offsetVal, iStopProbRespond', iNStop);
        goSSDTimepoints = offsetVal + (min(-goRTMinusSSD{iPropInd}) : max(-goRTMinusSSD{iPropInd}));
        inhFnGoMinusSSD{iPropInd} = weibull_curve(fitParameters, goSSDTimepoints);
        
        
        
        
        
        
        
        if plotFlag
            % Determine color to use for plot based on which checkerboard color
            % proportion being used.
            if options.collapseSignal
                cMap = ccm_colormap([0 1]);
            else
                cMap = ccm_colormap(pSignalArray);
            end
            inhColor = cMap(iPropInd,:);
            
            inhEachMap = ccm_colormap([0 1]);
            if pSignalArray(iPropInd) < .5
                iInhColor = inhEachMap(1,:);
            elseif pSignalArray(iPropInd) > .5
                iInhColor = inhEachMap(2,:);
            else
                iInhColor = [.1 .1 .1];
            end
            inhLineW = 2*max(choiceLevels)/choiceLevels(iPropInd);
            
            %             plot(ax(axInhEach), ssd{iPropInd}, iStopProbRespond, 'color', inhColor, 'linewidth', 2)
            plot(ax(axInhEach), ssdTimePoints, inhibitionFn{iPropInd}, 'color', iInhColor, 'linewidth', inhLineW)
            plot(ax(ssrtPRight), pSignalArray(iPropInd), ssrtIntegrationWeighted(iPropInd), '.', 'markersize', 30, 'color', inhColor)
            %             plot(ax(axInhRTSSD), -goRTMinusSSD{iPropInd}, iStopProbRespond, '.', 'color', inhColor, 'markersize', 15);%'linewidth', 2)
            %             plot(ax(axInhRTSSD), (goSSDTimepoints - offsetVal), inhFnGoMinusSSD{iPropInd}, '-', 'color', inhColor, 'linewidth', 2);%'linewidth', 2)
            %             plot(ax(SSDvPCorrect), xVal, p(1) * xVal + p(2), 'color', inhColor, 'linewidth', 2)
            plot(ax(SSDvSigStrength), pSignalArray(iPropInd), mean(conditionSSD{iPropInd}), '.',  'color', inhColor, 'markersize', 30)
            
            plot(ax(nTrialPerCoh), pSignalArray(iPropInd), length(goTotalRT{iPropInd,1}), '.', 'markersize', 30, 'color', goColor);
            plot(ax(nTrialPerCoh), pSignalArray(iPropInd), sum(nStopTarg(iPropInd,:))+sum(nStopDist(iPropInd,:)), '.', 'markersize', 30, 'color', stopColor);
            plot(ax(nTrialPerCoh), pSignalArray(iPropInd), sum(nStopStop(iPropInd,:)), '*', 'markersize', 25, 'color', stopColor);
            if options.USE_PRE_SSD
                
                %                 plot(ax(ssrtPRight), pSignalArray(iPropInd), ssrtGrand(iPropInd), '.', 'markersize', 30, 'color', 'r')
            end
            yyaxis left
            plot(ax(axWeibullParam), pSignalArray(iPropInd), weibullParams{iPropInd}(1), '.', 'markersize', 30, 'color', inhColor)
            ylim([0 600])
            yyaxis right
            plot(ax(axWeibullParam), pSignalArray(iPropInd), weibullParams{iPropInd}(2), '*', 'markersize', 20, 'color', 'b')
            ylim([0 12])
            
        end
        
        
        
        
        
        
    end
    
    
    
    
    
    
    
    
    
    
    % Regress the SSRT w.r.t. color coherence
    %     if options.USE_PRE_SSD
    %         regressSSRT = ssrtGrand;
    %     else
    regressSSRT = ssrtIntegrationWeighted;
    %     end
    
    if ~options.collapseSignal
        flipPropArray = pSignalArray;
        flipPropArray(flipPropArray > .5) = abs(1 - flipPropArray(flipPropArray > .5));
        [p, s] = polyfit(flipPropArray, regressSSRT, 1);
        [y, delta] = polyval(p, flipPropArray, s);
        stats = regstats(flipPropArray, regressSSRT);
        %         fprintf('p-value for regression: %.4f\n', stats.tstat.pval(2))
        R = corrcoef(flipPropArray, regressSSRT);
        Rsqrd = R(1, 2)^2;
        cov(flipPropArray, regressSSRT);
        xVal = min(flipPropArray) : .001 : .5;
        yVal = p(1) * xVal + p(2);
        regressionLineY = [yVal(1:end-1), fliplr(yVal)];
        regressionLineX = min(flipPropArray) : .001 : abs(1 - min(flipPropArray));
    end
    
    
    if plotFlag
        % inhibition function
        xlim(ax(axInhEach), [ssdTimePoints(1) ssdTimePoints(end)])
        ylim(ax(axInhEach), [0 1]);
        set(get(ax(axInhEach), 'ylabel'), 'String', 'p(Respond | stop)')
        
        
        xlim(ax(axInhRTSSD), [goSSDTimepoints(1) - offsetVal goSSDTimepoints(1) - offsetVal+400])
        
        % SSRT vs color coherence
        set(ax(ssrtPRight), 'Ylim', [min([0; 1.2*min(regressSSRT(:)); 1.2*min(ssrtIntegrationWeighted(:))]) max([0; 1.2*max(regressSSRT(:)); 1.2*max(ssrtIntegrationWeighted(:))])])
        set(get(ax(ssrtPRight), 'ylabel'), 'String', 'SSRT')
        plot(ax(ssrtPRight), [.5 .5], ylim, '--k')
        if ~options.collapseSignal
            % SSRT vs color coherence
            xlim(ax(ssrtPRight), [pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
            set(ax(ssrtPRight), 'xtick', pSignalArray)
            set(ax(ssrtPRight), 'xtickLabel', pSignalArray*100)
            plot(ax(ssrtPRight), regressionLineX, regressionLineY, 'b')
            
            % Trial Number vs color coherence
            xlim(ax(nTrialPerCoh), [pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
            set(ax(nTrialPerCoh), 'xtick', pSignalArray)
            set(ax(nTrialPerCoh), 'xtickLabel', pSignalArray*100)
            
            % Weibull parameters
            xlim(ax(axWeibullParam), [pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
            set(ax(axWeibullParam), 'xtick', pSignalArray)
            set(ax(axWeibullParam), 'xtickLabel', pSignalArray*100)
        else
            % SSRT vs color coherence
            xlim(ax(ssrtPRight), [pSignalArray(1) - choicePlotXMargin pSignalArray(2) + choicePlotXMargin])
            set(ax(ssrtPRight), 'xtick', [pSignalArray(1) pSignalArray(2)])
            set(ax(ssrtPRight), 'xtickLabel', {'left','right'})
            
            % Trial Number vs color coherence
            xlim(ax(nTrialPerCoh), [pSignalArray(1) - choicePlotXMargin pSignalArray(2) + choicePlotXMargin])
            set(ax(nTrialPerCoh), 'xtick', [pSignalArray(1) pSignalArray(2)])
            set(ax(nTrialPerCoh), 'xtickLabel', {'left','right'})
            
            % Weibull parameters
            xlim(ax(axWeibullParam), [pSignalArray(1) - choicePlotXMargin pSignalArray(2) + choicePlotXMargin])
            set(ax(axWeibullParam), 'xtick', [pSignalArray(1) pSignalArray(2)])
            set(ax(nTrialaxWeibullParamPerCoh), 'xtickLabel', {'left','right'})
        end
        
        
        
        %         % ssd vs p(correct)
        %         xlim(ax(SSDvPCorrect), [ssdArray(1)-ssdMargin ssdArray(end)+ssdMargin])
        %         set(ax(SSDvPCorrect), 'xtick', ssdArray)
        %         set(ax(SSDvPCorrect), 'xtickLabel', ssdArray)
        %         ylim(ax(SSDvPCorrect), [-.05 1.05]);
        %         set(get(ax(SSDvPCorrect), 'ylabel'), 'String', 'p(Correct)')
        %
        %         if options.collapseSignal
        %             xlim(ax(SSDvSigStrength), [pSignalArray(1) - choicePlotXMargin pSignalArray(2) + choicePlotXMargin])
        %             set(ax(SSDvSigStrength), 'xtick', [pSignalArray(1) pSignalArray(2)])
        %             set(ax(SSDvSigStrength), 'xticklabel', {'left','right'})
        %         else
        %             xlim(ax(SSDvSigStrength), [pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
        %             set(ax(SSDvSigStrength), 'xtick', pSignalArray)
        %         end
        
        % Mean SSD per signal strength
        set(get(ax(SSDvSigStrength), 'ylabel'), 'String', 'Mean SSD')
        set(ax(SSDvSigStrength), 'Ylim', [200 700])
        if ~options.collapseSignal
            xlim(ax(SSDvSigStrength), [pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
            set(ax(SSDvSigStrength), 'xtick', pSignalArray)
            set(ax(SSDvSigStrength), 'xtickLabel', pSignalArray*100)
        else
            xlim(ax(SSDvSigStrength), [pSignalArray(1) - choicePlotXMargin pSignalArray(2) + choicePlotXMargin])
            set(ax(SSDvSigStrength), 'xtick', [pSignalArray(1) pSignalArray(2)])
            set(ax(SSDvSigStrength), 'xtickLabel', {'left','right'})
        end
        
        
        
        if plotSurface
            figure(64)
            clf
            %      surf(ax(axPred), stopRespondRT)
            %      surf(ax(axPred), stopRespondRTPredict)
            colormap(hsv)
            %      surf(ssdArray, signalStrength, stopRespondRT)
            surface(stopRespondRT)
            hold on
        end
    end % if plotflag
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % Also calculate a grand inhibition function, SSRT, and predicted vs.
    % observed non-canceled RTs across discriminability levels
    goGrandRT = [];
    for iPropInd = 1 : nSignal
        goGrandRT = [goGrandRT; goTotalRT{iPropInd}];
    end
    goGrandRT = sort(goGrandRT);
    rtCum = 1/length(goGrandRT):1/length(goGrandRT):1; %y-axis of a cumulative prob dist
    
    % Grand Inhibition function
    nStopIncorrectGrand = zeros(nSSD, 1);
    nStopCorrectGrand = zeros(nSSD, 1);
    stopRespondProbGrand = nan(nSSD, 1);
    nStopGrand = zeros(nSSD, 1);
    for iSSDIndex = 1 : nSSD
        for jPropIndex = 1 : nSignal;
            nStopIncorrectGrand(iSSDIndex) = nStopIncorrectGrand(iSSDIndex) + length(stopRespondTrial{jPropIndex, iSSDIndex});
            nStopCorrectGrand(iSSDIndex) = nStopCorrectGrand(iSSDIndex) + length(stopStopTrial{jPropIndex, iSSDIndex});
        end
        nStopGrand(iSSDIndex) = nStopIncorrectGrand(iSSDIndex) + nStopCorrectGrand(iSSDIndex);
        % Inhibition function data points
        stopRespondProbGrand(iSSDIndex) = nStopIncorrectGrand(iSSDIndex) / nStopGrand(iSSDIndex);
        
        
        % Predicted mean stopIncorrectRT (based on goRTs and p(noncanceled):
        iRTIndex = find(rtCum >= stopRespondProbGrand(iSSDIndex), 1);   %match estimated p(noncan|SSD) to p(RT)
        stopRespondRTPredictGrand(iSSDIndex) = nanmean(goGrandRT(1 : iRTIndex-1));
        % Observed mean stop IncorrectRT:
        emptyData = cellfun(@isempty, stopRespondTrial(:, iSSDIndex));
        stopInocrrectTrial = cell2mat(stopRespondTrial(~emptyData, iSSDIndex));
        %         stopInocrrectTrial = cell2mat(cellfun(@(x) x', stopRespondTrial(:, iSSDIndex), 'uniformoutput', false));
        stopRespondRTObserveGrand(iSSDIndex) = nanmean(allRT(stopInocrrectTrial));
        
    end
    %     keepSSD = ~isnan(stopRespondProbGrand);
    %     stopRespondProbGrand = stopRespondProbGrand(keepSSD);
    %     nStopGrand = nStopGrand(keepSSD);
    %     ssdArrayGrand = ssdArray(keepSSD);
    ssdArrayGrand = ssdArray;
    %     [fitParameters, lowestSSE] = Weibull(ssdArrayGrand, stopRespondProbGrand, nStopGrand);
    [fitParameters, lowestSSE] = Weibull(ssdArray, stopRespondProbGrand, nStopGrand);
    ssdTimePoints = ssdArrayGrand(1) : ssdArrayGrand(end);
    inhibitionFnGrand = weibull_curve(fitParameters, ssdTimePoints);
    
    if plotFlag
        plot(ax(axInhGrand), ssdTimePoints, inhibitionFnGrand, 'color', 'g', 'linewidth', 2)
        plot(ax(axInhGrand), ssdArrayGrand, stopRespondProbGrand, '.k', 'markersize', 25)
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
    
    ssrtCollapse        = get_ssrt(ssdArrayGrand, stopRespondProbGrand, nStopGrand, goGrandRT, fitParameters);
    ssrtCollapseGrand   = ssrtCollapse.grand;
    ssrtCollapseIntegrationWeighted = ssrtCollapse.integrationWeighted;
    ssrtCollapseIntegration = ssrtCollapse.integration;
    ssrtCollapseMean    = ssrtCollapse.mean;
    
    if plotFlag && ~options.collapseSignal
        plot(ax(ssrtPRight), xlim, [ssrtCollapseIntegrationWeighted ssrtCollapseIntegrationWeighted], '--k')
        
    end
    
    
    
    %         if printPlot
    %             localFigurePath = local_figure_path;
    %             print(figureHandle,[localFigurePath, sessionID, '_ccm_inhibition'],'-dpdf', '-r300')
    %         end
    
    
    
    Data(kTarg).allRT                   = allRT;
    
    
    Data(kTarg).ssrtGrand                   = ssrtGrand;
    Data(kTarg).ssrtMean                    = ssrtMean;
    Data(kTarg).ssrtIntegration             = ssrtIntegration;
    Data(kTarg).ssrtIntegrationWeighted  	= ssrtIntegrationWeighted;
    Data(kTarg).ssrtIntegrationSimple     	= ssrtIntegrationSimple;
    Data(kTarg).ssrtCollapseGrand           = ssrtCollapseGrand;
    Data(kTarg).ssrtCollapseIntegrationWeighted 	= ssrtCollapseIntegrationWeighted;
    Data(kTarg).ssrtCollapseIntegration             = ssrtCollapseIntegration;
    Data(kTarg).ssrtCollapseMean            = ssrtCollapseMean;
    Data(kTarg).ssd                = ssd;
    
    Data(kTarg).goTotalRT       = goTotalRT;
    Data(kTarg).goTargRT       = goTargRT;
    Data(kTarg).goDistRT       = goDistRT;
    Data(kTarg).stopTargRT       = stopTargRT;
    Data(kTarg).stopDistRT       = stopDistRT;
    
    
    
    
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
    Data(kTarg).stopRespondRTObserveGrand = stopRespondRTObserveGrand;
    Data(kTarg).stopRespondRTPredict = stopRespondRTPredict;
    Data(kTarg).stopRespondRT = stopRespondRT;
    
    Data(kTarg).weibullParams = weibullParams;
    
end % for kTarg = 1 : nTargPair

clear trialData
end % function













function Data = ccm_psychometric(subjectID, sessionID, options)

% function [goTargRT, goDistRT, stopTargRT, stopDistRT] = ccm_chronometric(subjectID, sessionID, plotFlag)
%
% Psychometric analyses for choice countermanding task.
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
% Returns Data structure with fields:
%
%   nGo
%   nGoRight
%   nStopIncorrect
%   nStopIncorrectRight
%   nGoMatch
%   nGoRightMatch
%   goRightLogical
%   goRightSignalStrength
%   stopRightLogical
%   stopRightSignalStrength
%   goRightLogicalMatch
%   goRightSignalStrengthMatch
%   goPsychFn;
%   stopPsychFn;
%   goPsychFnMatch;


% Set default options or return a default options structure
if nargin < 3
    options.collapseTarg        = false;
    options.include50           = false;
    options.doStops            = true;
    options.USE_TWO_COLORS    = false;
    
    options.plotFlag            = true;
    options.printPlot           = false;
    options.figureHandle      	= 400;
    
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


%
% % ***********************************************************************
% % Psychometric Function: Proportion(Red Checker) vs Probability(go Right)
% % ***********************************************************************
%
%%
% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID, ccm_min_vars);
% [trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
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
distAngleArray = targAngleArray;
nTarg = length(targAngleArray);
nTrial = size(trialData.trialOutcome, 1);

if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end


allRT = trialData.responseOnset - trialData.responseCueOn;

if ~include50
    pSignalArray(pSignalArray == .5) = [];
end

% Which Signal Strength levels to analyze
% switch options.collapseSignal
%     case true
%         nSignal = 2;
%     case false
nSignal = length(pSignalArray);
% end


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
    
    
    if plotFlag
%         figureHandle = figureHandle + 1;
        % axes names
        pRightvSignal = 1;
        nRow = 3;
        nColumn = 4;
        if printPlot
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
        else
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
        end
        
        choicePlotXMargin = .03;
        ssdMargin = 20;
        ylimArray = [];
        
        ax(pRightvSignal) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 1) yAxesPosition(2, 1) axisWidth axisHeight]);
        hold(ax(pRightvSignal), 'on')
        goColor = [0 0 0];
        stopColor = [1 0 0];
        goColorMatch = [.5 .5 .5];
        
    end
    
    
    
    
    
    
    nTrial                   = size(trialData.trialOutcome, 1);
    goRightLogical           = nan(nTrial, 1);
    goRightSignalStrength    = nan(nTrial, 1);
    stopRightLogical         = nan(nTrial, length(ssdArray));
    stopRightSignalStrength  = nan(nTrial, length(ssdArray));
    goRightLogicalMatch      = nan(nTrial, 1);
    goRightSignalStrengthMatch = nan(nTrial, 1);
    
    
    
    stopRightProb        = zeros(length(pSignalArray), 1);
    nStopIncorrect       = nan(length(ssdArray), length(pSignalArray));
    nStopIncorrectRight  = nan(length(ssdArray), length(pSignalArray));
    
    goRightProb          = zeros(length(pSignalArray), 1);
    nGo                  = nan(1, length(pSignalArray));
    nGoRight             = nan(1, length(pSignalArray));
    
    goRightProbMatch     = zeros(length(pSignalArray), 1);
    nGoMatch             = nan(1, length(pSignalArray));
    nGoRightMatch        = nan(1, length(pSignalArray));
    
    % Get default trial selection options
    optSelect       = ccm_trial_selection;
    
    for iPropIndex = 1 : nSignal
        iPct = pSignalArray(iPropIndex) * 100;
        optSelect.rightCheckerPct = iPct;
        
        
        
        % If collapsing into all left and all right or all up/all down,
        % need to note here that there are "2" angles to deal with
        % (important for calling ccm_trial_selection.m)
        rightTargArray = targAngleArray(rightTargInd);
        leftDistArray = distAngleArray(rightTargInd);
        leftTargArray = targAngleArray(leftTargInd);
        rightDistArray = distAngleArray(leftTargInd);
        if options.collapseTarg && iPct(1) > 50
            kTargAngle = rightTargArray;
            kDistAngle = leftDistArray;
        elseif options.collapseTarg && iPct(1) < 50
            kTargAngle = leftTargArray;
            kDistAngle = rightDistArray;
        else
            if iPct(1) > 50
                kTargAngle = rightTargArray(kTarg);
                kDistAngle = leftDistArray(kTarg);
            elseif iPct(1) < 50
                kTargAngle = leftTargArray(kTarg);
                kDistAngle = rightDistArray(kTarg);
            end
        end
        
        
        % Get correct go probabilities rightward
        optSelect.ssd       = 'none';
        
        % All completed go trials
        optSelect.targDir     = [kTargAngle, kDistAngle];
        optSelect.responseDir  = 'collapse';
        optSelect.outcome     = {'goCorrectTarget', 'targetHoldAbort'};
        goTarget              = ccm_trial_selection(trialData, optSelect);
        optSelect.outcome     = {'goCorrectDistractor', 'distractorHoldAbort'};
        goDistractor          = ccm_trial_selection(trialData, optSelect);
        goTrial               = union(goTarget, goDistractor);
        goRightLogical(goTrial) = 0;
        nGo(iPropIndex)       = length(goTrial);
        
        % Rightward go trials
        optSelect.responseDir  = 'right';
        optSelect.outcome     = {'goCorrectTarget', 'targetHoldAbort'};
        goRightTarget         = ccm_trial_selection(trialData, optSelect);
        optSelect.outcome     = {'goCorrectDistractor', 'distractorHoldAbort'};
        goRightDistractor     = ccm_trial_selection(trialData, optSelect);
        goRightTrial          = union(goRightTarget, goRightDistractor);
        goRightLogical(goRightTrial) = 1;
        nGoRight(iPropIndex)  = length(goRightTrial);
        
        %     nGoRight(iPropIndex)
        %     nGo(iPropIndex)
        
        goRightProb(iPropIndex) = nGoRight(iPropIndex) / nGo(iPropIndex) ;
        goRightSignalStrength(goTrial) = trialData.targ1CheckerProp(goTrial);
        
        
        
        % Get incorrect stop probabilities rightward
        iStopTrial = [];
        if ~isempty(ssdArray) && options.doStops
            for jSSDIndex = 1 : length(ssdArray)
                jSSD = ssdArray(jSSDIndex);
                optSelect.ssd       = jSSD;
                
                % All stop incorrect trials
                optSelect.targDir     = [kTargAngle, kDistAngle];
                optSelect.responseDir  = 'collapse';
                optSelect.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget', ...
                    'distractorHoldAbort', 'stopIncorrectDistractor', 'stopIncorrectPreSSDDistractor'};
                stopTrial = ccm_trial_selection(trialData, optSelect);
                iStopTrial = [iStopTrial; stopTrial];
                stopRightLogical(stopTrial, jSSDIndex) = 0;
                nStopIncorrect(jSSDIndex, iPropIndex) = length(stopTrial);
                
                % All stop incorrect rightward trials
                optSelect.responseDir  = 'right';
                optSelect.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
                stopTargetRight = ccm_trial_selection(trialData, optSelect);
                %             optSelect.targDir       = 'left';
                optSelect.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
                stopDistractorRight = ccm_trial_selection(trialData, optSelect);
                stopTrialRight = union(stopTargetRight, stopDistractorRight);
                stopRightLogical(stopTrialRight, jSSDIndex) = 1;
                nStopIncorrectRight(jSSDIndex, iPropIndex) = length(stopTrialRight);
                
                stopRightSignalStrength(stopTrial, jSSDIndex) = trialData.targ1CheckerProp(stopTrial);
            end  % SSD loop
            stopRightProb(iPropIndex) = nansum(nStopIncorrectRight(:, iPropIndex)) / nansum(nStopIncorrect(:, iPropIndex)) ;
        end
        
        
        % Get lateny-matched go probabilities rightward
        goTrialMatch = goTrial;
        goTrialMatchRT = allRT(goTrialMatch);
        goRightTrialMatch = goRightTrial;
        stopTrialRT = allRT(iStopTrial);
        deleteTrial = [];
        while nanmean(goTrialMatchRT) > nanmean(stopTrialRT)
            % While mean of go RT distribution > stop RTs, remove latest go RT
            % trials
            [y,i] = max(goTrialMatchRT);
            iDeleteTrial = goTrialMatch(i);
            deleteTrial = [deleteTrial; iDeleteTrial];
            goTrialMatchRT(i) = [];
            goTrialMatch(i) = [];
            % If that trial was also a go right trial, get rid of it too
            goRightTrialMatch(goRightTrialMatch == iDeleteTrial) = [];
        end
        %     nanmean(goTrialMatchRT)
        %     nanmean(stopTrialRT)
        %     deleteTrial
        
        
        goRightLogicalMatch(goTrialMatch) = 0;
        goRightSignalStrengthMatch(goTrialMatch) = trialData.targ1CheckerProp(goTrialMatch);
        nGoMatch(iPropIndex) = length(goTrialMatch);
        
        goRightLogicalMatch(goRightTrialMatch) = 1;
        goRightSignalStrengthMatch(goRightTrialMatch) = trialData.targ1CheckerProp(goRightTrialMatch);
        nGoRightMatch(iPropIndex) = length(goRightTrialMatch);
        
        goRightProbMatch(iPropIndex) = nGoRightMatch(iPropIndex) / nGoMatch(iPropIndex) ;
        
        
        
    end
    % goRightLogical(isnan(goRightLogical)) = [];
    
    
    
    
    
    
    
    
    xData = pSignalArray(1) : .001 : pSignalArray(end);
    
    % Need to change this to a psychometric function
    [bestFitParams, maxLikelihood] = psychometric_weibull_fit(goRightSignalStrength, goRightLogical);
    goParams.threshold = bestFitParams(1);
    goParams.slope = bestFitParams(2);
    goParams.max = bestFitParams(3);
    goParams.min = bestFitParams(4);
    goPsychFn = weibull_curve(bestFitParams, xData);
    
    
    if ~isempty(ssdArray) && options.doStops
        [bestFitParams, maxLikelihood] = psychometric_weibull_fit(nanmean(stopRightSignalStrength, 2), nanmean(stopRightLogical, 2));
        stopParams.threshold = bestFitParams(1);
        stopParams.slope = bestFitParams(2);
        stopParams.max = bestFitParams(3);
        stopParams.min = bestFitParams(4);
        stopPsychFn = weibull_curve(bestFitParams, xData);
    else
        stopPsychFn = [];
        stopParams = [];
    end
    
    
    
    [bestFitParams, maxLikelihood] = psychometric_weibull_fit(goRightSignalStrengthMatch, goRightLogicalMatch);
    goMatchParams.threshold = bestFitParams(1);
    goMatchParams.slope = bestFitParams(2);
    goMatchParams.max = bestFitParams(3);
    goMatchParams.min = bestFitParams(4);
    goPsychFnMatch = weibull_curve(bestFitParams, xData);
    
    
    
    
    if plotFlag
        
        cla
        %     plot(ax(pRightvSignal), propPoints/100, goPsychometricFn, '-', 'color', goColor, 'linewidth', 2)
        
        %
        plot(ax(pRightvSignal), xData, goPsychFn, '-', 'color', goColor, 'linewidth', 2)
        plot(ax(pRightvSignal), xData, goPsychFnMatch, '-', 'color', goColorMatch, 'linewidth', 2)
        %
        
        plot(ax(pRightvSignal), pSignalArray, goRightProb, 'o', 'color', goColor, 'linewidth', 2, 'markerfacecolor', goColor, 'markeredgecolor', goColor)
        plot(ax(pRightvSignal), pSignalArray, goRightProbMatch, 'o', 'color', goColorMatch, 'linewidth', 2, 'markerfacecolor', goColorMatch, 'markeredgecolor', goColorMatch)
        if ~isempty(ssdArray) && options.doStops
            %         plot(ax(pRightvSignal), propPoints/100, stopPsychometricFn, '-', 'color', stopColor, 'linewidth', 2)
            plot(ax(pRightvSignal), xData, stopPsychFn, '-', 'color', stopColor, 'linewidth', 2)
            plot(ax(pRightvSignal), pSignalArray, stopRightProb, 'o', 'color', stopColor, 'linewidth', 2, 'markerfacecolor', stopColor, 'markeredgecolor', stopColor)
        end
        set(ax(pRightvSignal), 'xtick', pSignalArray)
        set(ax(pRightvSignal), 'xtickLabel', pSignalArray*100)
        set(get(ax(pRightvSignal), 'ylabel'), 'String', 'p(Right)')
        set(ax(pRightvSignal),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
        set(ax(pRightvSignal),'YLim',[0 1])
        plot(ax(pRightvSignal), [.5 .5], ylim, '--k')
        
        
        
        % ANOVA calculations
        if ~isempty(ssdArray) && options.doStops
            anovaData = [];
            groupInh = {};
            groupSig = [];
            goTarg = goRightProb;
            stopTarg = stopRightProb;
            [h,p] = ttest2(goTarg, stopTarg);
            fprintf('Psychometric T-test:\nStop vs. Go: \tp = %.4f\n', p)
            
        end
        
        %         if printPlot
        %             localFigurePath = local_figure_path;
        %             print(figureHandle,[localFigurePath, sessionID, '_', Unit(kUnitIndex).name,'_ccm_psychometric'],'-dpdf', '-r300')
        %         end
        
    end % if plotFlag
    
    
    Data(kTarg).nGo            = nGo;
    Data(kTarg).nGoRight       = nGoRight;
    Data(kTarg).nStopIncorrect = nStopIncorrect;
    Data(kTarg).nStopIncorrectRight = nStopIncorrectRight;
    Data(kTarg).nGoMatch            = nGoMatch;
    Data(kTarg).nGoRightMatch       = nGoRightMatch;
    Data(kTarg).goRightLogical = goRightLogical;
    Data(kTarg).goRightSignalStrength = goRightSignalStrength;
    Data(kTarg).stopRightLogical = stopRightLogical;
    Data(kTarg).stopRightSignalStrength = stopRightSignalStrength;
    Data(kTarg).goRightLogicalMatch = goRightLogicalMatch;
    Data(kTarg).goRightSignalStrengthMatch = goRightSignalStrengthMatch;
    
    Data(kTarg).goPsychFn      = goPsychFn;
    Data(kTarg).stopPsychFn    = stopPsychFn;
    Data(kTarg).goPsychFnMatch      = goPsychFnMatch;
    Data(kTarg).goParams       = goParams;
    Data(kTarg).stopParams      = stopParams;
    Data(kTarg).goMatchParams    = goMatchParams;
    
    Data(kTarg).pSignalArray = pSignalArray;
end % kTarg

clear trialData


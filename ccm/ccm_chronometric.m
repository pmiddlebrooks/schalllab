function Data = ccm_chronometric(subjectID, sessionID, options)
% function [goTargRT, goDistRT, stopTargRT, stopDistRT] = ccm_chronometric(subjectID, sessionID, plotFlag)
%
% Chronometric analyses for choice countermanding task.
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
%   pSignalArrayLeft
%   pSignalArrayRight
%   goLeftToTarg
%   goRightToTarg
%   goLeftToDist
%   goRightToDist
%   stopLeftToTarg
%   stopRightToTarg
%   stopLeftToDist
%   stopRightToDist



% Set default options or return a default options structure
if nargin < 3
    options.collapseTarg        = false;
    options.include50           = false;
    options.doStops              = true;
    options.USE_TWO_COLORS         = false;
    
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

%%
% Load the data
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
distAngleArray = targAngleArray;


% Flag to determine whether we want to include stop trial analyses for the
% session
MIN_RT = 120;
MAX_RT = 1200;
nSTD   = 3;

if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end




% Truncate RTs
[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTD);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);
allRT(outlierTrial) = [];

% Get rid of trials that have mircorsaccades before the response (but after
% go cue)
% msData = ccm_microsaccade_before_RT(subjectID, sessionID);
% allRT(msData.msTrial) = nan;

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
    
    
    % axes names
    pRightvRT = 1;
    
    
    if plotFlag
%         figureHandle = figureHandle + 1;
        nRow = 3;
        nColumn = 4;
        screenOrSave = 'save';
        if printPlot
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
        else
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
        end
        choicePlotXMargin = .03;
        ssdMargin = 20;
        ylimArray = [];
        
    end
    
    
    
    
    nTrial = size(trialData.trialOutcome, 1);
    goColor = [0 0 0];
    stopColor = [1 0 0];
    
    
    
    % ***********************************************************************
    % Chronometric Function:     Proportion(Right) vs RT
    % ***********************************************************************
    
    if plotFlag
        ax(pRightvRT) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2) yAxesPosition(2, 2) axisWidth axisHeight]);
        cla
        hold(ax(pRightvRT), 'on')
         % Target RTs - Distractor RTs for go and stop trials
        %     ax(targDistRT) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 2) yAxesPosition(3, 2) axisWidth axisHeight]);
        %     hold(ax(targDistRT), 'on')
    end
    
    
    
    
    
    % COLLECT RT DATA PER CONDITION
    % ===========================================
    
    % Get correct go RTs
    
    goTargRT            = cell(1, length(pSignalArray));
    goDistRT            = cell(1, length(pSignalArray));
    go50TargRightRT     = cell(1);
    go50TargLeftRT      = cell(1);
    go50DistRightRT     = cell(1);
    go50DistLeftRT      = cell(1);
    nGoTrial            = nan(1, length(pSignalArray));
    nNonCancelTrial     = nan(1, length(pSignalArray));
    nCancelTrial        = nan(1, length(pSignalArray));
    
    if ~isempty(ssdArray) && options.doStops
        stopTargRT = cell(length(ssdArray), length(pSignalArray));
        stopDistRT = cell(length(ssdArray), length(pSignalArray));
        stop50TargRightRT = cell(length(ssdArray), 1);
        stop50TargLeftRT = cell(length(ssdArray), 1);
        stop50DistRightRT = cell(length(ssdArray), 1);
        stop50DistLeftRT = cell(length(ssdArray), 1);
    end
    % Get default trial selection options
    optSelect       = ccm_trial_selection;
    
    for iPropIndex = 1 : length(pSignalArray);
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
        
        
        
        
        % ****** GO TRIALS  ******
        optSelect.ssd       = 'none';
        
        
        
        if iPct ~= 50
            optSelect.targDir    	= kTargAngle;
            optSelect.outcome     = {'goCorrectTarget', 'targetHoldAbort'};
            iGoTargTrial          = ccm_trial_selection(trialData, optSelect);
            optSelect.outcome     = {'goCorrectDistractor', 'distractorHoldAbort'};
            iGoDistTrial = ccm_trial_selection(trialData, optSelect);
            %       iGoTargTrial(ismember(iGoTargTrial, outlierTrial)) = [];
            %       iGoDistTrial(ismember(iGoDistTrial, outlierTrial)) = [];
            
            goTargRT{iPropIndex}  = allRT(iGoTargTrial);
            goDistRT{iPropIndex}  = allRT(iGoDistTrial);
            nGoTrial(iPropIndex)  = length(iGoTargTrial) + length(iGoDistTrial);
        elseif iPct == 50
            optSelect.targDir   	= rightTargArray(kTarg);
            optSelect.outcome     = {'goCorrectTarget', 'targetHoldAbort'};
            go50TargRight         = ccm_trial_selection(trialData, optSelect);
            optSelect.targDir   	= leftTargArray(kTarg);
            go50TargLeft          = ccm_trial_selection(trialData, optSelect);
            optSelect.targDir    	= rightTargArray(kTarg);
            optSelect.outcome     = {'goCorrectDistractor', 'distractorHoldAbort'};
            go50DistRight         = ccm_trial_selection(trialData, optSelect);
            optSelect.targDir    	= leftTargArray(kTarg);
            go50DistLeft          = ccm_trial_selection(trialData, optSelect);
            
            %       go50TargRight(ismember(go50TargRight, outlierTrial)) = [];
            %       go50TargLeft(ismember(go50TargLeft, outlierTrial)) = [];
            %       go50DistRight(ismember(go50DistRight, outlierTrial)) = [];
            %       go50DistLeft(ismember(go50DistLeft, outlierTrial)) = [];
            
            go50TargRightRT{1}    = allRT(go50TargRight);
            go50TargLeftRT{1}     = allRT(go50TargLeft);
            go50DistRightRT{1}    = allRT(go50DistRight);
            go50DistLeftRT{1}     = allRT(go50DistLeft);
            
        end
        
        %     oddData = find(isnan(cell2mat(trialData.saccadeToTargetIndex)) & iGoTargIndices & ismember(cell2mat(trialData.targ1CheckerProp), pSignalArray));
        %     if oddData
        %         fprintf('%d trials to target %d are listed as %s but don''t have valid saccades to target:\n', length(oddData), iTarget, goOutcome)
        %         disp(oddData)
        %     end
        %     iGoTargIndices(oddData) = 0;
        
        
        %    end
        
        
        
        
        % ****** STOP TRIALS  ******
        
        % Stop Incorrect RTs
        if ~isempty(ssdArray) && options.doStops
            
            %       for iPropIndex = 1 : length(pSignalArray);
            %          iPct = pSignalArray(iPropIndex) * 100;
            %          optSelect.rightCheckerPct = iPct;
            %
            for jSSDIndex = 1 : length(ssdArray)
                jSSD = ssdArray(jSSDIndex);
                optSelect.ssd       = jSSD;
                if iPct ~= 50
                    optSelect.targDir       = kTargAngle;
                    optSelect.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
                    stopTargTrial           = ccm_trial_selection(trialData, optSelect);
                    optSelect.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
                    stopDistTrial           = ccm_trial_selection(trialData, optSelect);
                    %             stopTargTrial(ismember(stopTargTrial, outlierTrial)) = [];
                    %             stopDistTrial(ismember(stopDistTrial, outlierTrial)) = [];
                    iStopTargIndices        = zeros(nTrial, 1);
                    iStopDistIndices        = zeros(nTrial, 1);
                    iStopTargIndices(stopTargTrial) = 1;
                    iStopDistIndices(stopDistTrial) = 1;
                    
                    stopTargRT{jSSDIndex, iPropIndex} = allRT(stopTargTrial);
                    stopDistRT{jSSDIndex, iPropIndex} = allRT(stopDistTrial);
                    
                elseif iPct == 50
                    optSelect.targDir       = rightTargArray(kTarg);
                    optSelect.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
                    stop50TargRight         = ccm_trial_selection(trialData, optSelect);
                    optSelect.targDir       = leftTargArray(kTarg);
                    stop50TargLeft          = ccm_trial_selection(trialData, optSelect);
                    optSelect.targDir       = rightTargArray(kTarg);
                    optSelect.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
                    stop50DistRight         = ccm_trial_selection(trialData, optSelect);
                    optSelect.targDir       = leftTargArray(kTarg);
                    stop50DistLeft          = ccm_trial_selection(trialData, optSelect);
                    
                    %             stop50TargRight(ismember(stop50TargRight, outlierTrial)) = [];
                    %             stop50TargLeft(ismember(stop50TargLeft, outlierTrial)) = [];
                    %             stop50DistRight(ismember(stop50DistRight, outlierTrial)) = [];
                    %             stop50DistLeft(ismember(stop50DistLeft, outlierTrial)) = [];
                    
                    stop50TargRightRT{jSSDIndex} = allRT(stop50TargRight);
                    stop50TargLeftRT{jSSDIndex} = allRT(stop50TargLeft);
                    stop50DistRightRT{jSSDIndex} = allRT(stop50DistRight);
                    stop50DistLeftRT{jSSDIndex} = allRT(stop50DistLeft);
                end
            end
            
       
            %         oddData = find(isnan(cell2mat(trialData.saccadeToTargetIndex)) & iStopTargIndices & ismember(cell2mat(trialData.targ1CheckerProp), pSignalArray));
            %         if oddData
            %             fprintf('%d trials to target %d are listed as Incorrect Stops but don''t have valid saccades to target:\n', length(oddData), iTarget)
            %             disp(oddData)
            %         end
            %         iStopTargIndices(oddData) = 0;
            
            
        end
    end
    
    
    
    
%     goTargFastRT = cellfun(@(x) x(x < nanmedian(x)), goTargRT, 'uniformoutput', false);
    
    
    
    pSignalArrayLeft      = pSignalArray(pSignalArray < .5);
    pSignalArrayRight     = pSignalArray(pSignalArray > .5);
    fiftyPct        = pSignalArray(pSignalArray == .5);
    
    goLeftToTarg          = goTargRT(pSignalArray < .5);
    goRightToTarg         = goTargRT(pSignalArray > .5);
    goLeftToDist          = goDistRT(pSignalArray > .5);
    goRightToDist         = goDistRT(pSignalArray < .5);
    if ~isempty(fiftyPct)
        goLeft50TargPct  = go50TargLeftRT;
        goRight50TargPct  = go50TargRightRT;
        goLeft50DistPct  = go50DistLeftRT;
        goRight50DistPct  = go50DistRightRT;
    end
    
    goLeftToTargMean      = cellfun(@nanmean, goLeftToTarg);
    goRightToTargMean     = cellfun(@nanmean, goRightToTarg);
    goLeftToDistMean      = cellfun(@nanmean, goLeftToDist);
    goRightToDistMean     = cellfun(@nanmean, goRightToDist);
    if ~isempty(fiftyPct)
        goLeft50TargPctMean     = cellfun(@nanmean, goLeft50TargPct);
        goRight50TargPctMean = cellfun(@nanmean, goRight50TargPct);
        goLeft50DistPctMean = cellfun(@nanmean, goLeft50DistPct);
        goRight50DistPctMean = cellfun(@nanmean, goRight50DistPct);
    end
    goLeftToTargStd      = cellfun(@nanstd, goLeftToTarg);
    goRightToTargStd     = cellfun(@nanstd, goRightToTarg);
    goLeftToDistStd      = cellfun(@nanstd, goLeftToDist);
    goRightToDistStd     = cellfun(@nanstd, goRightToDist);
    if ~isempty(fiftyPct)
        goLeft50TargPctStd = cellfun(@nanstd, goLeft50TargPct);
        goRight50TargPctStd = cellfun(@nanstd, goRight50TargPct);
        goLeft50DistPctStd = cellfun(@nanstd, goLeft50DistPct);
        goRight50DistPctStd = cellfun(@nanstd, goRight50DistPct);
    end
    
    
    
    if ~isempty(ssdArray) && options.doStops
        stopLeftToTarg            = stopTargRT(:, pSignalArray < .5);
        stopRightToTarg           = stopTargRT(:, pSignalArray > .5);
        stopLeftToDist            = stopDistRT(:, pSignalArray > .5);
        stopRightToDist           = stopDistRT(:, pSignalArray < .5);
        if ~isempty(fiftyPct)
            stopLeft50TargPct  = stop50TargLeftRT;
            stopRight50TargPct  = stop50TargRightRT;
            stopLeft50DistPct  = stop50DistLeftRT;
            stopRight50DistPct  = stop50DistRightRT;
        end
        
        for iPropIndex = 1 : length(pSignalArrayLeft)
            stopLeftToTargMean(iPropIndex)        = nanmean(cell2mat(stopLeftToTarg(:, iPropIndex)));
            stopRightToDistMean(iPropIndex)       = nanmean(cell2mat(stopRightToDist(:, iPropIndex)));
            
            stopLeftToTargStd(iPropIndex)        = nanstd(cell2mat(stopLeftToTarg(:, iPropIndex)));
            stopRightToDistStd(iPropIndex)       = nanstd(cell2mat(stopRightToDist(:, iPropIndex)));
        end
        for iPropIndex = 1 : length(pSignalArrayRight)
            stopRightToTargMean(iPropIndex)       = nanmean(cell2mat(stopRightToTarg(:, iPropIndex)));
            stopLeftToDistMean(iPropIndex)        = nanmean(cell2mat(stopLeftToDist(:, iPropIndex)));
            
            stopRightToTargStd(iPropIndex)       = nanstd(cell2mat(stopRightToTarg(:, iPropIndex)));
            stopLeftToDistStd(iPropIndex)        = nanstd(cell2mat(stopLeftToDist(:, iPropIndex)));
        end
        
        
        if ~isempty(fiftyPct)
            stopLeft50TargPctMean = nanmean(cell2mat(stopLeft50TargPct));
            stopRight50TargPctMean = nanmean(cell2mat(stopRight50TargPct));
            stopLeft50DistPctMean = nanmean(cell2mat(stopLeft50DistPct));
            stopRight50DistPctMean = nanmean(cell2mat(stopRight50DistPct));
            
            stopLeft50TargPctStd = nanstd(cell2mat(stopLeft50TargPct));
            stopRight50TargPctStd = nanstd(cell2mat(stopRight50TargPct));
            stopLeft50DistPctStd = nanstd(cell2mat(stopLeft50DistPct));
            stopRight50DistPctStd = nanstd(cell2mat(stopRight50DistPct));
        end
    end
    
    if ~isempty(fiftyPct)
        pSignalArrayLeft = [pSignalArrayLeft; fiftyPct];
        pSignalArrayRight = [fiftyPct; pSignalArrayRight];
        
        goLeftToTarg = [goLeftToTarg, goLeft50TargPct];
        goRightToTarg = [goRight50TargPct, goRightToTarg];
        goLeftToTargMean = [goLeftToTargMean, goLeft50TargPctMean];
        goRightToTargMean = [goRight50TargPctMean, goRightToTargMean];
        goLeftToTargStd = [goLeftToTargStd, goLeft50TargPctStd];
        goRightToTargStd = [goRight50TargPctStd, goRightToTargStd];
        
        goLeftToDist = [goLeft50DistPct, goLeftToDist];
        goRightToDist = [goRightToDist, goRight50DistPct];
        goLeftToDistMean = [goLeft50DistPctMean, goLeftToDistMean];
        goRightToDistMean = [goRightToDistMean, goRight50DistPctMean];
        goLeftToDistStd = [goLeft50DistPctStd, goLeftToDistStd];
        goRightToDistStd = [goRightToDistStd, goRight50DistPctStd];
        
        if ~isempty(ssdArray) && options.doStops
            stopLeftToTarg = [stopLeftToTarg, stopLeft50TargPct];
            stopRightToTarg = [stopRight50TargPct, stopRightToTarg];
            stopLeftToTargMean = [stopLeftToTargMean, stopLeft50TargPctMean];
            stopRightToTargMean = [stopRight50TargPctMean, stopRightToTargMean];
            stopLeftToTargStd = [stopLeftToTargStd, stopLeft50TargPctStd];
            stopRightToTargStd = [stopRight50TargPctStd, stopRightToTargStd];
            
            stopLeftToDist = [stopLeft50DistPct, stopLeftToDist];
            stopRightToDist = [stopRightToDist, stopRight50DistPct];
            stopLeftToDistMean = [stopLeft50DistPctMean, stopLeftToDistMean];
            stopRightToDistMean = [stopRightToDistMean, stopRight50DistPctMean];
            stopLeftToDistStd = [stopLeft50DistPctStd, stopLeftToDistStd];
            stopRightToDistStd = [stopRightToDistStd, stopRight50DistPctStd];
        end
    end
    
    
    
    
    
    
    
    
    
    % ANOVA calculations
    % ==================================
    if ~isempty(ssdArray) && options.doStops
        anovaData = [];
        groupInh = {};
        groupSig = [];
        pSignalArray = [pSignalArrayLeft; pSignalArrayRight];
        goTarg = [goLeftToTarg, goRightToTarg];
        stopTarg = [stopLeftToTarg, stopRightToTarg];
        for i = 1 : size(goTarg, 2)
            anovaData = [anovaData; goTarg{i}];
            groupInh = [groupInh; repmat({'go'}, length(goTarg{i}), 1)];
            groupSig = [groupSig; repmat(i, length(goTarg{i}), 1)];
        end
        for i = 1 : size(stopTarg, 2)
            anovaData = [anovaData; cell2mat(stopTarg(:,i))];
            groupInh = [groupInh; repmat({'stop'}, length(cell2mat(stopTarg(:,i))), 1)];
            groupSig = [groupSig; repmat(i, length(cell2mat(stopTarg(:,i))), 1)];
        end
        [p,table,stats] = anovan(anovaData,{groupInh, groupSig}, 'display', 'off');
        fprintf('RT ANOVA:\nStop vs. Go: \t\tp = %d\nSignal Strength: \tp = %d\n', p(1), p(2))
        
    end
    
    
    Data(kTarg).ssdArray    = ssdArray;
    Data(kTarg).pSignalArrayLeft    = pSignalArrayLeft;
    Data(kTarg).pSignalArrayRight   = pSignalArrayRight;
    Data(kTarg).goLeftToTarg          = goLeftToTarg;
    Data(kTarg).goRightToTarg         = goRightToTarg;
    Data(kTarg).goLeftToDist          = goLeftToDist;
    Data(kTarg).goRightToDist         = goRightToDist;
    if options.doStops
        Data(kTarg).stopLeftToTarg        = stopLeftToTarg;
        Data(kTarg).stopRightToTarg       = stopRightToTarg;
        Data(kTarg).stopLeftToDist        = stopLeftToDist;
        Data(kTarg).stopRightToDist       = stopRightToDist;
    else
        Data(kTarg).stopLeftToTarg        = [];
        Data(kTarg).stopRightToTarg       = [];
        Data(kTarg).stopLeftToDist        = [];
        Data(kTarg).stopRightToDist       = [];
    end
    
    % PLOTTING
    % ===========================
    if plotFlag
        plot(ax(pRightvRT), pSignalArrayLeft, goLeftToTargMean, '-o', 'color', goColor, 'linewidth', 1, 'markerfacecolor', goColor, 'markeredgecolor', goColor)
        plot(ax(pRightvRT), pSignalArrayRight, goRightToTargMean, '-o', 'color', goColor, 'linewidth', 1, 'markerfacecolor', goColor, 'markeredgecolor', goColor)
        plot(ax(pRightvRT), pSignalArrayRight, goLeftToDistMean, 'o', 'color', goColor, 'markeredgecolor', goColor)
        plot(ax(pRightvRT), pSignalArrayLeft, goRightToDistMean, 'o', 'color', goColor, 'markeredgecolor', goColor)
        % errorbar(ax(pRightvRT), pSignalArrayLeft ,goLeftToTargMean, goLeftToTargStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
        % errorbar(ax(pRightvRT), pSignalArrayRight ,goRightToTargMean, goRightToTargStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
        % errorbar(ax(pRightvRT), pSignalArray(pSignalArray < .5), goLeftToDistMean, goLeftToDistStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
        % errorbar(ax(pRightvRT), pSignalArray(pSignalArray > .5), goRightToDistMean, goRightToDistStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
        %     plot(ax(pRightvRT), pSignalArray, goDistRT, 'o', 'color', goColor, 'markeredgecolor', goColor)
        %     plot(ax(targDistRT), pSignalArray, goTargRT - goDistRT, 'o', 'color', goColor, 'markerfacecolor', goColor, 'markeredgecolor', goColor)
        if ~isempty(ssdArray) && options.doStops
            %         plot(ax(pRightvRT), pSignalArrayLeft(1:length(stopLeftToTargMean)), stopLeftToTargMean, '-o', 'color', stopColor, 'linewidth', 1, 'markerfacecolor', stopColor, 'markeredgecolor', stopColor)
            %         plot(ax(pRightvRT), pSignalArrayRight(end+1-length(stopRightToTargMean) : end), stopRightToTargMean, '-o', 'color', stopColor, 'linewidth', 1, 'markerfacecolor', stopColor, 'markeredgecolor', stopColor)
            %         plot(ax(pRightvRT), pSignalArrayRight(end+1-length(stopLeftToDistMean) : end), stopLeftToDistMean, 'o', 'color', stopColor, 'markeredgecolor', stopColor)
            %         plot(ax(pRightvRT), pSignalArrayLeft(1:length(stopRightToDistMean)), stopRightToDistMean, 'o', 'color', stopColor, 'markeredgecolor', stopColor)
            plot(ax(pRightvRT), pSignalArrayLeft, stopLeftToTargMean, '-o', 'color', stopColor, 'linewidth', 1, 'markerfacecolor', stopColor, 'markeredgecolor', stopColor)
            plot(ax(pRightvRT), pSignalArrayRight, stopRightToTargMean, '-o', 'color', stopColor, 'linewidth', 1, 'markerfacecolor', stopColor, 'markeredgecolor', stopColor)
            plot(ax(pRightvRT), pSignalArrayRight, stopLeftToDistMean, 'o', 'color', stopColor, 'markeredgecolor', stopColor)
            plot(ax(pRightvRT), pSignalArrayLeft, stopRightToDistMean, 'o', 'color', stopColor, 'markeredgecolor', stopColor)
            % errorbar(ax(pRightvRT), pSignalArrayLeft ,stopLeftToTargMean, stopLeftToTargStd, '.' , 'linestyle' , 'none', 'color', stopColor, 'linewidth' , 1)
            % errorbar(ax(pRightvRT), pSignalArrayRight ,stopRightToTargMean, stopRightToTargStd, '.' , 'linestyle' , 'none', 'color', stopColor, 'linewidth' , 1)
            % errorbar(ax(pRightvRT), pSignalArray(pSignalArray < .5), stopLeftToDistMean, stopLeftToDistStd, '.' , 'linestyle' , 'none', 'color', stopColor, 'linewidth' , 1)
            % errorbar(ax(pRightvRT), pSignalArray(pSignalArray > .5), stopRightToDistMean, stopRightToDistStd, '.' , 'linestyle' , 'none', 'color', stopColor, 'linewidth' , 1)
            %         plot(ax(pRightvRT), pSignalArray, stopDistRT, 'o', 'color', stopColor, 'markeredgecolor', stopColor)
            %         plot(ax(targDistRT), pSignalArray, stopTargRT - stopDistRT, 'o', 'color', stopColor, 'markerfacecolor', stopColor, 'markeredgecolor', stopColor)
            % legend(ax(pRightvRT), {'Go Target', 'Stop Target', 'Go Distractor', 'Stop Distractor'}, 'location', 'southeast');
        end
        
        set(ax(pRightvRT), 'xtick', pSignalArray)
        set(ax(pRightvRT), 'xtickLabel', pSignalArray*100)
        set(get(ax(pRightvRT), 'ylabel'), 'String', 'RT')
        set(ax(pRightvRT),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
        set(ax(pRightvRT),'YLim',[min([goLeftToTargMean, goRightToTargMean, goLeftToDistMean, goRightToDistMean]) - 100, max([goLeftToTargMean, goRightToTargMean, goLeftToDistMean, goRightToDistMean]) + 25])
        %               set(ax(pRightvRT),'YLim',[200 600])
        plot(ax(pRightvRT), [.5 .5], ylim, '--k')
        
        %     set(ax(targDistRT), 'xtick', pSignalArray)
        %     set(ax(targDistRT), 'xtickLabel', pSignalArray*100)
        %     set(get(ax(targDistRT), 'ylabel'), 'String', 'dRT (Correct - Error)')
        %     set(ax(targDistRT),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
        %     plot(ax(targDistRT), [.5 .5], ylim, '--k')
        %     plot(ax(targDistRT), [pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin], [0 0], '--k')
        
        %         if printPlot
        %             localFigurePath = local_figure_path;
        %             print(figureHandle,[localFigurePath, sessionID, '_ccm_chronometric'],'-dpdf', '-r300')
        %         end
    end
    
end % kTarg

clear trialData


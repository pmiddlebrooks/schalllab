function [stopProbRespond, stopProbTarget, inhibitionFn, ssrt, goCorrectProbRight, stopIncorrectProbRight, goCorrectTargRT, stopIncorrectTargRT] = ccm_behavior(subjectID, sessionID)
%%
stopProbRespond = [];
stopProbTarget = [];
inhibitionFn = [];
ssrt = [];
stopIncorrectProbRight = [];
stopIncorrectTargRT = [];

% subjectID = 'bz';
% sessionID = '0905';
% Load the data
[dataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessionID);
% If the file hasn't already been copied to a local directory, do it now
if exist(localDataFile, 'file') ~= 2
    copyfile(dataFile, localDataPath)
end
load(localDataFile);

DO_STOPS = 1;

if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end


% Convert cells to doubles if necessary
trialData = cell_to_mat(trialData);



% sessionID = 'bp032n04';
plotFlag = 1;

% axes names
pRightvRT = 1;
pRightvProbRight = 2;
inhibition = 3;
SSDvPCorrect = 4;
ssrtPRight = 5;
targDistRT = 6;
pRightVsPCorrect = 7;




if plotFlag
    %     nPlotColumn = 2;
    %     nPlotRow = 3;
    %     [figureHandle, cmWidth, cmHeight] = standard_figure;
    nRow = 3;
    nColumn = 3;
    screenOrSave = 'screen';
    [figureHandle, axisWidth, axisHeight, xAxesPosition, yAxesPosition] = big_figure(nRow, nColumn, screenOrSave);
    
    choicePlotXMargin = .03;
    ssdMargin = 20;
    ylimArray = [];
    
    %     if errorFlag
    %         xError = fill([1:length(condx_avg) length(condx_avg):-1:1], [(condx_avg - condx_sem) condx_avg(end:-1:1) + condx_sem(end:-1:1)], [1 .4 .4]);
    %         set(xError, 'edgecolor', 'none');
    %         % alpha(xError, .3)
    %         yError = fill([1:length(condy_avg) length(condy_avg):-1:1], [(condy_avg - condy_sem) condy_avg(end:-1:1) + condy_sem(end:-1:1)], [.4 .4 1]);
    %         set(yError, 'edgecolor', 'none');
    %         % alpha(yError, .3)
    %     end
end




nTrial = size(trialData.trialOutcome, 1);
goCorrectColor = [0 0 0];
stopIncorrectColor = [1 0 0];
target1ProportionArray = unique(trialData.targ1CheckerProp);
target1ProportionArray(isnan(target1ProportionArray)) = [];

% Need to do a little SSD value adjusting, due to ms difference and 1-frame
% differences in SSD values
ssdArray = unique(trialData.stopSignalOn - trialData.responseCueOn);
ssdArray(isnan(ssdArray)) = [];
if ~isempty(ssdArray) && DO_STOPS
    a = diff(ssdArray);
    ssdArray(a == 1) = ssdArray(a == 1) + 1;
    ssdArray = unique(ssdArray);
    b = [ssdArray(1); diff(ssdArray)];
    ssdArray(b < 18) = []
end











% ***********************************************************************
% Inhibition Function:
%       &
% SSD vs. Proportion of Response trials
%       &
% SSD vs. Proportion(Correct Choice)
% ***********************************************************************

% If there weren't stop trials, skip all stop-related analyses
if ~isempty(ssdArray) && DO_STOPS
    if plotFlag
        % inhibition function
        ax(inhibition) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
        hold(ax(inhibition))
        % p(right) vs ssrt
        ax(ssrtPRight) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 1) yAxesPosition(2, 1) axisWidth axisHeight]);
        hold(ax(ssrtPRight))
        % SSD vs p(correct)
        ax(SSDvPCorrect) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 1) yAxesPosition(3, 1) axisWidth axisHeight]);
        hold(ax(SSDvPCorrect))
    end
    
    stopStopTotal    = cell(length(target1ProportionArray), length(ssdArray));
    stopIncorrectTotal  = cell(length(target1ProportionArray), length(ssdArray));
    goCorrectTargTotal  = cell(length(target1ProportionArray), 1);
    goCorrectTotalRT    = cell(length(target1ProportionArray), 1);
    stopProbRespond     = nan(length(target1ProportionArray), length(ssdArray));
    stopProbTarget      = nan(length(target1ProportionArray), length(ssdArray));
    inhibitionFn        = cell(length(target1ProportionArray), 1);
    goCorrectTargRT     = zeros(length(target1ProportionArray), 1);
    nStop               = nan(length(target1ProportionArray), length(ssdArray));
    ssrt                = nan(length(target1ProportionArray), 1);
    minColorGun = .25;
    maxColorGun = 1;
    for iPropIndex = 1 : length(target1ProportionArray);
        iPercent = target1ProportionArray(iPropIndex) * 100;
        
        % Determine color to use for plot based on which checkerboard color
        % proportion being used. Normalize the available color spectrum to do
        % it
        if iPercent == 50
            inhColor = [0 0 0];
        elseif iPercent < 50
            colorNorm = .5 - target1ProportionArray(1);
            colorProp = (.5 - target1ProportionArray(iPropIndex)) / colorNorm;
            colorGun = minColorGun + (maxColorGun - minColorGun) * colorProp;
            %         colorProp = (.7 - target1ProportionArray(iPropIndex)) / .7;
            inhColor = [0 colorGun colorGun];
        elseif iPercent > 50
            colorNorm = target1ProportionArray(end) - .5;
            colorProp = (target1ProportionArray(iPropIndex) - .5) / colorNorm;
            colorGun = minColorGun + (maxColorGun - minColorGun) * colorProp;
            %         colorProp = abs(.3 - target1ProportionArray(iPropIndex)) / .7;
            inhColor = [colorGun 0 colorGun];
        end
        
        for jSSDIndex = 1 : length(ssdArray)
            tic
            jSSD = ssdArray(jSSDIndex);
            
            % stop correct trials
            stopStopTrial = ccm_trial_selection(subjectID, sessionID,  {'stopCorrect'}, iPercent, jSSD, 'all');
            stopStopTotal{iPropIndex, jSSDIndex} = stopStopTrial;  % Keep track of totals for grand inhibition fnct
            
            % stop incorrect trials
            stopIncorrectTargetOutcome = {'targetHoldAbort', 'stopIncorrectTarget'};
            stopIncorrectDistractorOutcome = {'distractorHoldAbort', 'stopIncorrectDistractor'};
            
            stopIncorrectTargTrial = ccm_trial_selection(subjectID, sessionID,  stopIncorrectTargetOutcome, iPercent, jSSD, 'all');
            stopIncorrectDistTrial = ccm_trial_selection(subjectID, sessionID,  stopIncorrectDistractorOutcome, iPercent, jSSD, 'all');
            
            % stop incorrect trials for inhibition: do we want stop incorrect to the target or to
            % target/distractor?
            stopIncorrectTrial = union(stopIncorrectTargTrial, stopIncorrectDistTrial);
            stopIncorrectTotal{iPropIndex, jSSDIndex} = stopIncorrectTrial;  % Keep track of totals for grand inhibition fnct
            
            % Inhibition function data points:
            stopProbRespond(iPropIndex, jSSDIndex) = length(stopIncorrectTrial) / (length(stopStopTrial) + length(stopIncorrectTrial));
            nStop(iPropIndex, jSSDIndex) = length(stopStopTrial) + length(stopIncorrectTrial);
            
            % p(Correct choice) vs. SSD data points:
            stopProbTarget(iPropIndex, jSSDIndex) = length(stopIncorrectTargTrial) / (length(stopIncorrectTargTrial) + length(stopIncorrectDistTrial));
            toc
        end % iSSDIndex
        
        
        % Inhibition function data points
        iStopProbRespond = stopProbRespond(iPropIndex, :);
        keepSSD = ~isnan(iStopProbRespond);
        iStopProbRespond = iStopProbRespond(keepSSD);
        
        % p(Correct choice) vs. SSD data points:
        iStopProbTarget = stopProbTarget(iPropIndex, :);
        keepSSDChoice = ~isnan(iStopProbTarget);
        iStopProbTarget = iStopProbTarget(keepSSDChoice);
        iSSDArrayChoice = ssdArray(keepSSDChoice)';
        [p, s] = polyfit(iSSDArrayChoice, iStopProbTarget, 1);
        xVal = min(iSSDArrayChoice) : max(iSSDArrayChoice);
        
        
        % Inhibition function calculation
        iNStop = nStop(iPropIndex, :);
        iNStop = iNStop(keepSSD);
        iSSDArray = ssdArray(keepSSD);
        
        [fitParameters, lowestSSE] = Weibull(iSSDArray, iStopProbRespond, iNStop);
        %         [fitParameters, lowestSSE] = Weibull_fast(iSSDArray, iStopProbRespond, iNStop);
        timePoints = iSSDArray(1) : iSSDArray(end);
        inhibitionFn{iPropIndex} = weibull_curve(fitParameters, timePoints);
        
        
        % SSRT: get go RTs and number of stop trials (already have other
        % necessary variables)
        goCorrectOutcome = {'goCorrectTarget'};
        goCorrectTarget = ccm_trial_selection(subjectID, sessionID, goCorrectOutcome, iPercent, 'none', 'all');
        iGoCorrectTargIndices = zeros(nTrial, 1);
        iGoCorrectTargIndices(goCorrectTarget) = 1;
        
        oddData = find(isnan(trialData.saccToTargetIndex) & iGoCorrectTargIndices & ismember(trialData.targ1CheckerProp, target1ProportionArray));
        if oddData
            fprintf('%d trials to target %d are listed as %s but don''t have valid saccades to target:\n', length(oddData), iTarget, goCorrectOutcome)
            disp(oddData)
        end
        iGoCorrectTargIndices(oddData) = 0;
        
        if sum(iGoCorrectTargIndices)
            iGoCorrectTarg = find(iGoCorrectTargIndices);
            goCorrectTargTotal{iPropIndex} = iGoCorrectTarg;  % Keep track of totals for grand inhibition fnct
            
            responseOnset = trialData.responseOnset(iGoCorrectTarg);
            responseCueOn = trialData.responseCueOn(iGoCorrectTarg);
            iGoCorrectTargRT = responseOnset - responseCueOn;
            goCorrectTotalRT{iPropIndex} = iGoCorrectTargRT;
            goCorrectTargRT(iPropIndex) = mean(iGoCorrectTargRT);
        end
        iNStop = nStop(iPropIndex, :)';
        iNStop(iNStop == 0) = [];
        [SSRT both_SSRTs SSRT_c SSRT_r meanSSD] = get_SSRT(iSSDArray, iStopProbRespond, iNStop, iGoCorrectTargRT, fitParameters);
        %         ssrt(iPropIndex) = SSRT_c;
        ssrt(iPropIndex) = SSRT;
        
        
        if plotFlag
            plot(ax(inhibition), timePoints, inhibitionFn{iPropIndex}, 'color', inhColor, 'linewidth', 2)
            plot(ax(ssrtPRight), target1ProportionArray(iPropIndex), ssrt(iPropIndex), '.', 'markersize', 30, 'color', inhColor)
            plot(ax(SSDvPCorrect), xVal, p(1) * xVal + p(2), 'color', inhColor, 'linewidth', 2)
        end
        
    end % iPropIndex
    
    if plotFlag
        % inhibition function
        xlim(ax(inhibition), [ssdArray(1)-ssdMargin ssdArray(end)+ssdMargin])
        set(ax(inhibition), 'xtick', ssdArray)
        set(ax(inhibition), 'xtickLabel', ssdArray)
        ylim(ax(inhibition), [0 1]);
        set(get(ax(inhibition), 'ylabel'), 'String', 'p(Respond | stop)')
        % p(correct) vs ssrt
        xlim(ax(ssrtPRight), [target1ProportionArray(1) - choicePlotXMargin target1ProportionArray(end) + choicePlotXMargin])
        set(ax(ssrtPRight), 'xtick', target1ProportionArray)
        set(ax(ssrtPRight), 'xtickLabel', target1ProportionArray*100)
        set(get(ax(ssrtPRight), 'ylabel'), 'String', 'SSRT')
        plot(ax(ssrtPRight), [.5 .5], ylim, '--k')
        % ssd vs p(correct)
        xlim(ax(SSDvPCorrect), [ssdArray(1)-ssdMargin ssdArray(end)+ssdMargin])
        set(ax(SSDvPCorrect), 'xtick', ssdArray)
        set(ax(SSDvPCorrect), 'xtickLabel', ssdArray)
        ylim(ax(SSDvPCorrect), [-.05 1.05]);
        set(get(ax(SSDvPCorrect), 'ylabel'), 'String', 'p(Correct)')
        
    end % if plotflag
    
    
    
    
    
    % Also calculate a grand inhibition function and SSRT, across
    % discriminability levels
    
    % Grand Inhibition function
    nStopIncorrectGrand = zeros(length(ssdArray), 1);
    nStopCorrectGrand = zeros(length(ssdArray), 1);
    stopProbRespondGrand = nan(length(ssdArray), 1);
    nStopGrand = zeros(length(ssdArray), 1);
    for iSSDIndex = 1 : length(ssdArray)
        for jPropIndex = 1 : length(target1ProportionArray);
            nStopIncorrectGrand(iSSDIndex) = nStopIncorrectGrand(iSSDIndex) + length(stopIncorrectTotal{jPropIndex, iSSDIndex});
            nStopCorrectGrand(iSSDIndex) = nStopCorrectGrand(iSSDIndex) + length(stopStopTotal{jPropIndex, iSSDIndex});
        end
        nStopGrand(iSSDIndex) = nStopIncorrectGrand(iSSDIndex) + nStopCorrectGrand(iSSDIndex);
        % Inhibition function data points
        stopProbRespondGrand(iSSDIndex) = nStopIncorrectGrand(iSSDIndex) / nStopGrand(iSSDIndex);
    end
    keepSSD = ~isnan(stopProbRespondGrand);
    stopProbRespondGrand = stopProbRespondGrand(keepSSD);
    nStopGrand = nStopGrand(keepSSD);
    ssdArray = ssdArray(keepSSD);
    [fitParameters, lowestSSE] = Weibull(ssdArray, stopProbRespondGrand, nStopGrand);
    % [fitParameters, lowestSSE] = Weibull_fast(ssdArray, stopProbRespondGrand, nStopGrand);
    timePoints = ssdArray(1) : ssdArray(end);
    inhibitionFnGrand = weibull_curve(fitParameters, timePoints);
    
    if plotFlag
        plot(ax(inhibition), timePoints, inhibitionFnGrand, 'color', 'g', 'linewidth', 2)
        plot(ax(inhibition), ssdArray, stopProbRespondGrand, '.k', 'markersize', 25)
    end
    
    % Grand SSRT
    goCorrectGrandRT = [];
    for iPropIndex = 1 : length(target1ProportionArray)
        goCorrectGrandRT = [goCorrectGrandRT; goCorrectTotalRT{iPropIndex}];
    end
    [SSRT both_SSRTs SSRT_c SSRT_r meanSSD] = get_SSRT(ssdArray, stopProbRespondGrand, nStopGrand, goCorrectGrandRT, fitParameters);
    ssrtGrand = SSRT
    [ssdArray, nStopGrand, stopProbRespondGrand]
    %     SSRT
    %     SSRT_r
    %     meanSSD
    %     stopIncorrectTotal
    %     stopStopTotal
    %     goCorrectTargTotal
    
end % if ~isempty(ssdArray)






% ***********************************************************************
% Probability(Rightward response) vs Proportion(Correct Choice)
% ***********************************************************************
goCorrectTargTotal      = cell(length(target1ProportionArray), 1);
goCorrectDistTotal      = cell(length(target1ProportionArray), 1);
nGoCorrectTargTotal      = nan(length(target1ProportionArray), 1);
nGoCorrectDistTotal      = nan(length(target1ProportionArray), 1);
for iPropIndex = 1 : length(target1ProportionArray);
    iPercent = target1ProportionArray(iPropIndex) * 100;
    
    goCorrectDistOutcome =  {'goCorrectDistractor'};
    goCorrectDist = ccm_trial_selection(subjectID, sessionID, goCorrectDistOutcome, iPercent, 'none', 'all');
    iGoCorrectDistIndices = zeros(nTrial, 1);
    iGoCorrectDistIndices(goCorrectDist) = 1;
    
    oddData = find(isnan(trialData.saccToTargetIndex) & iGoCorrectDistIndices & ismember(trialData.targ1CheckerProp, target1ProportionArray));
    if oddData
        fprintf('%d trials are listed as %s but don''t have valid saccades to target:\n', length(oddData), goCorrectDistOutcome)
        disp(oddData)
    end
    iGoCorrectDistIndices(oddData) = 0;
    
    if sum(iGoCorrectDistIndices)
        iGoCorrectDist = find(iGoCorrectDistIndices);
        goCorrectDistTotal{iPropIndex} = iGoCorrectDist;  % Keep track of totals for grand inhibition fnct
        nGoCorrectDistTotal(iPropIndex) = length(iGoCorrectDist);  % Keep track of totals for grand inhibition fnct
    end
    
    
    goCorrectTargtOutcome =  {'goCorrectTarget'};
    goCorrectTarg = ccm_trial_selection(subjectID, sessionID, goCorrectTargtOutcome, iPercent, 'none', 'all');
    iGoCorrectTargIndices = zeros(nTrial, 1);
    iGoCorrectTargIndices(goCorrectTarg) = 1;
    
    oddData = find(isnan(trialData.saccToTargIndex) & iGoCorrectTargIndices & ismember(trialData.targ1CheckerProp, target1ProportionArray));
    if oddData
        fprintf('%d trials are listed as %s but don''t have valid saccades to target:\n', length(oddData), goCorrectTargtOutcome)
        disp(oddData)
    end
    iGoCorrectTargIndices(oddData) = 0;
    
    if sum(iGoCorrectTargIndices)
        iGoCorrectTarg = find(iGoCorrectTargIndices);
        goCorrectTargTotal{iPropIndex} = iGoCorrectTarg;  % Keep track of totals for grand inhibition fnct
        nGoCorrectTargTotal(iPropIndex) = length(iGoCorrectTarg);  % Keep track of totals for grand inhibition fnct
    end
end

pGoCorrectTarg = nGoCorrectTargTotal ./ (nGoCorrectTargTotal + nGoCorrectDistTotal);
pGoCorrectTarg(isnan(pGoCorrectTarg)) = 1;

if plotFlag
    % p(rightward response) vs p(correct choice)
    ax(pRightVsPCorrect) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 3) yAxesPosition(1, 3) axisWidth axisHeight]);
    hold(ax(pRightVsPCorrect))
    plot(ax(pRightVsPCorrect), target1ProportionArray, pGoCorrectTarg, 'color', 'k', 'linewidth', 2)
    
    set(ax(pRightVsPCorrect), 'xtick', target1ProportionArray)
    set(ax(pRightVsPCorrect), 'xtickLabel', target1ProportionArray*100)
    set(get(ax(pRightVsPCorrect), 'ylabel'), 'String', 'p(Correct)')
    set(ax(pRightVsPCorrect),'XLim',[target1ProportionArray(1) - choicePlotXMargin target1ProportionArray(end) + choicePlotXMargin])
    set(ax(pRightVsPCorrect),'YLim',[.45 1])
    plot(ax(pRightVsPCorrect), [.5 .5], ylim, '--k')
end

minColorGun = .25;
maxColorGun = 1;
if ~isempty(ssdArray) && DO_STOPS
    for jSSDIndex = 1 : length(ssdArray)
        % Determine color to use for plot based on which checkerboard color
        % proportion being used. Normalize the available color spectrum to do
        % it
        %             colorNorm = .5 - target1ProportionArray(1);
        %             colorProp = (.5 - target1ProportionArray(iPropIndex)) / colorNorm;
        ssdFrac = ssdArray(jSSDIndex) / ssdArray(end);
        colorGun = minColorGun + (maxColorGun - minColorGun) * ssdFrac;
        %         colorProp = (.7 - target1ProportionArray(iPropIndex)) / .7;
        ssdColor = [colorGun 0 0];
        
        plot(ax(pRightVsPCorrect), target1ProportionArray, stopProbTarget(:, jSSDIndex), 'color', ssdColor, 'linewidth', 2)
    end
end










% ***********************************************************************
% Psychometric Function: Proportion(Red Checker) vs Probability(go Right)
% ***********************************************************************

if plotFlag
    ax(pRightvProbRight) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
end


% Get correct go probabilities rightward
goCorrectProbRight = zeros(length(target1ProportionArray), 1);
for iPropIndex = 1 : length(target1ProportionArray);
    iPercent = target1ProportionArray(iPropIndex) * 100;
    
    % All go Correct trials
    goCorrectTarget = ccm_trial_selection(subjectID, sessionID, {'goCorrectTarget'}, iPercent, 'none', 'all');
    goCorrectDistractor = ccm_trial_selection(subjectID, sessionID, {'goCorrectDistractor'}, iPercent, 'none', 'all');
    goCorrect = union(goCorrectTarget, goCorrectDistractor);
    
    % Rightward go correct trials
    goCorrectRightTarget = ccm_trial_selection(subjectID, sessionID, {'goCorrectTarget'}, iPercent, 'none', 'right');
    goCorrectRightDistractor = ccm_trial_selection(subjectID, sessionID, {'goCorrectDistractor'}, iPercent, 'none', 'left');
    goCorrectRight = union(goCorrectRightTarget, goCorrectRightDistractor);
    
    goCorrectProbRight(iPropIndex) = length(goCorrectRight) / length(goCorrect) ;
end




% Get incorrect stop probabilities rightward
if ~isempty(ssdArray) && DO_STOPS
    stopIncorrectProbRight = zeros(length(target1ProportionArray), 1);
    for iPropIndex = 1 : length(target1ProportionArray);
        iPercent = target1ProportionArray(iPropIndex) * 100;
        
        % All stop incorrect trials
        stopStopOutcome = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectDistractor', 'distractorHoldAbort'};
        stopIncorrectTrial = ccm_trial_selection(subjectID, sessionID,  stopStopOutcome, iPercent, 'all', 'all');
        
        % All stop incorrect rightward trials
        stopTargetRight = ccm_trial_selection(subjectID, sessionID,  {'stopIncorrectTarget', 'targetHoldAbort'}, iPercent, 'all', 'right');
        stopDistractorRight = ccm_trial_selection(subjectID, sessionID,  {'stopIncorrectDistractor', 'distractorHoldAbort'}, iPercent, 'all', 'left');
        stopIncorrectTrialRight = union(stopTargetRight, stopDistractorRight);
        
        stopIncorrectProbRight(iPropIndex) = length(stopIncorrectTrialRight) / length(stopIncorrectTrial) ;
    end
end


if plotFlag
    plot(ax(pRightvProbRight), target1ProportionArray, goCorrectProbRight, '-o', 'color', goCorrectColor, 'linewidth', 2, 'markerfacecolor', goCorrectColor, 'markeredgecolor', goCorrectColor)
    hold on
    if ~isempty(ssdArray) && DO_STOPS
        plot(ax(pRightvProbRight), target1ProportionArray, stopIncorrectProbRight, '-o', 'color', stopIncorrectColor, 'linewidth', 2, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
    end
    set(ax(pRightvProbRight), 'xtick', target1ProportionArray)
    set(ax(pRightvProbRight), 'xtickLabel', target1ProportionArray*100)
    set(get(ax(pRightvProbRight), 'ylabel'), 'String', 'p(Right)')
    set(ax(pRightvProbRight),'XLim',[target1ProportionArray(1) - choicePlotXMargin target1ProportionArray(end) + choicePlotXMargin])
    plot(ax(pRightvProbRight), [.5 .5], ylim, '--k')
end








% ***********************************************************************
% Chronometric Function:     Proportion(Right) vs RT
% ***********************************************************************

if plotFlag
    ax(pRightvRT) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2) yAxesPosition(2, 2) axisWidth axisHeight]);
    hold(ax(pRightvRT))
    % Target RTs - Distractor RTs for go and stop trials
    ax(targDistRT) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 2) yAxesPosition(3, 2) axisWidth axisHeight]);
    hold(ax(targDistRT))
    
end


% Get correct go RTs
goCorrectOutcome =  {'goCorrectTarget'; 'targetHoldAbort'};
ssdRange = 'none';
targetHemifield = 'all';

goCorrectTargRT = cell(1, length(target1ProportionArray));
goCorrectDistRT = cell(1, length(target1ProportionArray));
for iPropIndex = 1 : length(target1ProportionArray);
    iPercent = target1ProportionArray(iPropIndex) * 100;
    
    if iPercent ~= 50
        goCorrectTarget = ccm_trial_selection(subjectID, sessionID, {'goCorrectTarget'; 'targetHoldAbort'}, iPercent, ssdRange, targetHemifield);
        goCorrectDistractor = ccm_trial_selection(subjectID, sessionID, {'goCorrectDistractor', 'distractorHoldAbort'}, iPercent, ssdRange, targetHemifield);
        iGoCorrectTargIndices = zeros(nTrial, 1);
        iGoCorrectDistIndices = zeros(nTrial, 1);
        iGoCorrectTargIndices(goCorrectTarget) = 1;
        iGoCorrectDistIndices(goCorrectDistractor) = 1;
        
        if sum(iGoCorrectTargIndices)
            iGoCorrectTarg = find(iGoCorrectTargIndices);
            responseOnset = trialData.responseOnset(iGoCorrectTarg);
            responseCueOn = trialData.responseCueOn(iGoCorrectTarg);
            iGoCorrectTargRT = responseOnset - responseCueOn;
            goCorrectTargRT{iPropIndex} = iGoCorrectTargRT;
        end
        if sum(iGoCorrectDistIndices)
            iGoCorrectDist = find(iGoCorrectDistIndices);
            responseOnset = trialData.responseOnset(iGoCorrectDist);
            responseCueOn = trialData.responseCueOn(iGoCorrectDist);
            iGoCorrectDistRT = responseOnset - responseCueOn;
            goCorrectDistRT{iPropIndex} = iGoCorrectDistRT;
        end
        
    elseif iPercent == 50
        goCorrect50Right = ccm_trial_selection(subjectID, sessionID, {'goCorrectTarget', 'targetHoldAbort', 'goCorrectDistractor', 'distractorHoldAbort'}, iPercent, ssdRange, 'right');
        goCorrect50Left = ccm_trial_selection(subjectID, sessionID, {'goCorrectTarget', 'targetHoldAbort', 'goCorrectDistractor', 'distractorHoldAbort'}, iPercent, ssdRange, 'left');
        iGoCorrect50RightIndices = zeros(nTrial, 1);
        iGoCorrect50LeftIndices = zeros(nTrial, 1);
        iGoCorrect50RightIndices(goCorrect50Right) = 1;
        iGoCorrect50LeftIndices(goCorrect50Left) = 1;
        
        if sum(iGoCorrect50RightIndices)
            iGoCorrect50Right = find(iGoCorrect50RightIndices);
            responseOnset = trialData.responseOnset(iGoCorrect50Right);
            responseCueOn = trialData.responseCueOn(iGoCorrect50Right);
            goCorrect50RightRT = responseOnset - responseCueOn;
        else
            goCorrect50RightRT = [];
        end
        if sum(iGoCorrect50LeftIndices)
            iGoCorrect50Left = find(iGoCorrect50LeftIndices);
            responseOnset = trialData.responseOnset(iGoCorrect50Left);
            responseCueOn = trialData.responseCueOn(iGoCorrect50Left);
            goCorrect50LeftRT = responseOnset - responseCueOn;
        else
            goCorrect50LeftRT = [];
        end
        
    end
    
    %     oddData = find(isnan(cell2mat(trialData.saccToTargIndex)) & iGoCorrectTargIndices & ismember(cell2mat(trialData.targ1CheckerProp), target1ProportionArray));
    %     if oddData
    %         fprintf('%d trials to target %d are listed as %s but don''t have valid saccades to target:\n', length(oddData), iTarget, goCorrectOutcome)
    %         disp(oddData)
    %     end
    %     iGoCorrectTargIndices(oddData) = 0;
    
    
end


% Get stop Incorrect RTs
ssdRange = 'all';
targetHemifield = 'all';
if ~isempty(ssdArray) && DO_STOPS
    stopIncorrectTargRT = cell(1, length(target1ProportionArray));
    stopIncorrectDistRT = cell(1, length(target1ProportionArray));
    
    for iPropIndex = 1 : length(target1ProportionArray);
        iPercent = target1ProportionArray(iPropIndex) * 100;
        
        if iPercent ~= 50
            stopIncorrectTargTrial = ccm_trial_selection(subjectID, sessionID,  {'stopIncorrectTarget', 'targetHoldAbort'}, iPercent, ssdRange, targetHemifield);
            stopIncorrectDistTrial = ccm_trial_selection(subjectID, sessionID,  {'stopIncorrectDistractor', 'distractorHoldAbort'}, iPercent, ssdRange, targetHemifield);
            iStopTargIndices = zeros(nTrial, 1);
            iStopDistIndices = zeros(nTrial, 1);
            iStopTargIndices(stopIncorrectTargTrial) = 1;
            iStopDistIndices(stopIncorrectDistTrial) = 1;
            
            if sum(iStopTargIndices)
                iStopIncorrectTarg  = find(iStopTargIndices);
                responseOnset = trialData.responseOnset(iStopIncorrectTarg);
                responseCueOn    = trialData.responseCueOn(iStopIncorrectTarg);
                iStopIncorrectTargRT = responseOnset - responseCueOn;
                stopIncorrectTargRT{iPropIndex} = iStopIncorrectTargRT;
            end
            if sum(iStopDistIndices)
                iStopIncorrectDist  = find(iStopDistIndices);
                responseOnset = trialData.responseOnset(iStopIncorrectDist);
                responseCueOn    = trialData.responseCueOn(iStopIncorrectDist);
                iStopIncorrectDistRT = responseOnset - responseCueOn;
                stopIncorrectDistRT{iPropIndex} = iStopIncorrectDistRT;
            end
        elseif iPercent == 50
            stopIncorrect50Right = ccm_trial_selection(subjectID, sessionID, {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectDistractor', 'distractorHoldAbort'}, iPercent, ssdRange, 'right');
            stopIncorrect50Left = ccm_trial_selection(subjectID, sessionID, {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectDistractor', 'distractorHoldAbort'}, iPercent, ssdRange, 'left');
            iStopIncorrect50RightIndices = zeros(nTrial, 1);
            iStopIncorrect50LeftIndices = zeros(nTrial, 1);
            iStopIncorrect50RightIndices(stopIncorrect50Right) = 1;
            iStopIncorrect50LeftIndices(stopIncorrect50Left) = 1;
            
            if sum(iStopIncorrect50RightIndices)
                iStopIncorrect50Right = find(iStopIncorrect50RightIndices);
                responseOnset = trialData.responseOnset(iStopIncorrect50Right);
                responseCueOn = trialData.responseCueOn(iStopIncorrect50Right);
                stopIncorrect50RightRT = responseOnset - responseCueOn;
            else
                stopIncorrect50RightRT = [];
            end
            if sum(iStopIncorrect50LeftIndices)
                iStopIncorrect50Left = find(iStopIncorrect50LeftIndices);
                responseOnset = trialData.responseOnset(iStopIncorrect50Left);
                responseCueOn = trialData.responseCueOn(iStopIncorrect50Left);
                stopIncorrect50LeftRT = responseOnset - responseCueOn;
            else
                stopIncorrect50LeftRT = [];
            end
            
        end
        
        
        
        %         oddData = find(isnan(cell2mat(trialData.saccToTargIndex)) & iStopTargIndices & ismember(cell2mat(trialData.targ1CheckerProp), target1ProportionArray));
        %         if oddData
        %             fprintf('%d trials to target %d are listed as Incorrect Stops but don''t have valid saccades to target:\n', length(oddData), iTarget)
        %             disp(oddData)
        %         end
        %         iStopTargIndices(oddData) = 0;
        
        
    end
end


targetLeft     = target1ProportionArray(target1ProportionArray < .5);
targetRight    = target1ProportionArray(target1ProportionArray > .5);
fiftyPercent    = target1ProportionArray(target1ProportionArray == .5);

goLeftToTarg          = goCorrectTargRT(target1ProportionArray < .5);
goRightToTarg         = goCorrectTargRT(target1ProportionArray > .5);
goLeftToDist          = goCorrectDistRT(target1ProportionArray < .5);
goRightToDist         = goCorrectDistRT(target1ProportionArray > .5);
goLeft50Percent  = goCorrect50LeftRT;
goRight50Percent  = goCorrect50RightRT;

goLeftToTargMean      = cellfun(@mean, goLeftToTarg);
goRightToTargMean     = cellfun(@mean, goRightToTarg);
goLeftToDistMean      = cellfun(@mean, goLeftToDist);
goRightToDistMean     = cellfun(@mean, goRightToDist);
goLeft50PercentMean = mean(goLeft50Percent);
goRight50PercentMean = mean(goRight50Percent);

if ~isempty(ssdArray) && DO_STOPS
    stopLeftToTarg            = stopIncorrectTargRT(:, target1ProportionArray < .5);
    stopRightToTarg           = stopIncorrectTargRT(:, target1ProportionArray > .5);
    stopLeftToDist            = stopIncorrectDistRT(:, target1ProportionArray < .5);
    stopRightToDist           = stopIncorrectDistRT(:, target1ProportionArray > .5);
    stopLeft50Percent  = stopIncorrect50LeftRT;
    stopRight50Percent  = stopIncorrect50RightRT;
    
    stopLeftToTargMean        = cellfun(@mean, stopLeftToTarg);
    stopRightToTargMean       = cellfun(@mean, stopRightToTarg);
    stopLeftToDistMean        = cellfun(@mean, stopLeftToDist);
    stopRightToDistMean       = cellfun(@mean, stopRightToDist);
    stopLeft50PercentMean = mean(stopLeft50Percent);
    stopRight50PercentMean = mean(stopRight50Percent);
end

if ~isempty(fiftyPercent)
    targetLeft = [targetLeft; fiftyPercent];
    targetRight = [fiftyPercent; targetRight];
    goLeftToTargMean = [goLeftToTargMean, goLeft50PercentMean];
    goRightToTargMean = [goRight50PercentMean, goRightToTargMean];
    if ~isempty(ssdArray) && DO_STOPS
        stopLeftToTargMean = [stopLeftToTargMean, stopLeft50PercentMean];
        stopRightToTargMean = [stopRight50PercentMean, stopRightToTargMean];
    end
end



if plotFlag
    plot(ax(pRightvRT), targetLeft, goLeftToTargMean, '-o', 'color', goCorrectColor, 'linewidth', 2, 'markerfacecolor', goCorrectColor, 'markeredgecolor', goCorrectColor)
    plot(ax(pRightvRT), targetRight, goRightToTargMean, '-o', 'color', goCorrectColor, 'linewidth', 2, 'markerfacecolor', goCorrectColor, 'markeredgecolor', goCorrectColor)
    plot(ax(pRightvRT), target1ProportionArray(target1ProportionArray < .5), goLeftToDistMean, 'o', 'color', goCorrectColor, 'markeredgecolor', goCorrectColor)
    plot(ax(pRightvRT), target1ProportionArray(target1ProportionArray > .5), goRightToDistMean, 'o', 'color', goCorrectColor, 'markeredgecolor', goCorrectColor)
    %     plot(ax(pRightvRT), target1ProportionArray, goCorrectDistRT, 'o', 'color', goCorrectColor, 'markeredgecolor', goCorrectColor)
    %     plot(ax(targDistRT), target1ProportionArray, goCorrectTargRT - goCorrectDistRT, 'o', 'color', goCorrectColor, 'markerfacecolor', goCorrectColor, 'markeredgecolor', goCorrectColor)
    if ~isempty(ssdArray) && DO_STOPS
        plot(ax(pRightvRT), targetLeft, stopLeftToTargMean, '-o', 'color', stopIncorrectColor, 'linewidth', 2, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
        plot(ax(pRightvRT), targetRight, stopRightToTargMean, '-o', 'color', stopIncorrectColor, 'linewidth', 2, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
        plot(ax(pRightvRT), target1ProportionArray(target1ProportionArray < .5), stopLeftToDistMean, 'o', 'color', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
        plot(ax(pRightvRT), target1ProportionArray(target1ProportionArray > .5), stopRightToDistMean, 'o', 'color', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
        %         plot(ax(pRightvRT), target1ProportionArray, stopIncorrectDistRT, 'o', 'color', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
        %         plot(ax(targDistRT), target1ProportionArray, stopIncorrectTargRT - stopIncorrectDistRT, 'o', 'color', stopIncorrectColor, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
    end
    
    set(ax(pRightvRT), 'xtick', target1ProportionArray)
    set(ax(pRightvRT), 'xtickLabel', target1ProportionArray*100)
    set(get(ax(pRightvRT), 'ylabel'), 'String', 'RT')
    set(ax(pRightvRT),'XLim',[target1ProportionArray(1) - choicePlotXMargin target1ProportionArray(end) + choicePlotXMargin])
    set(ax(pRightvRT),'YLim',[min([goLeftToTargMean, goRightToTargMean, goLeftToDistMean, goRightToDistMean]) - 50, max([goLeftToTargMean, goRightToTargMean, goLeftToDistMean, goRightToDistMean]) + 50])
    plot(ax(pRightvRT), [.5 .5], ylim, '--k')
    
    set(ax(targDistRT), 'xtick', target1ProportionArray)
    set(ax(targDistRT), 'xtickLabel', target1ProportionArray*100)
    set(get(ax(targDistRT), 'ylabel'), 'String', 'dRT (Correct - Error)')
    set(ax(targDistRT),'XLim',[target1ProportionArray(1) - choicePlotXMargin target1ProportionArray(end) + choicePlotXMargin])
    plot(ax(targDistRT), [.5 .5], ylim, '--k')
    plot(ax(targDistRT), [target1ProportionArray(1) - choicePlotXMargin target1ProportionArray(end) + choicePlotXMargin], [0 0], '--k')
end


print(gcf, ['~/schalllab/tempfigures/',sessionID, '_Behavior'], '-dpdf')
delete(localDataFile);




function post_session_psychometric(trialData)


plotFlag = 1;
nTrial = size(trialData.trialOutcome, 1);
targ1PropArray = unique(trialData.targ1CheckerProp);
 targ1PropArray(isnan(targ1PropArray)) = [];

% ***********************************************************************
% Psychometric Function: Proportion(Red Checker) vs Probability(go Right)
% ***********************************************************************

% if plotFlag
%     ax(pRightvProbRight) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
% end

% Get correct go probabilities rightward
goCorrectProbRight = zeros(length(targ1PropArray), 1);
for iPropIndex = 1 : length(targ1PropArray);   
    
    % All go Correct trials
    trialLogical = ones(nTrial, 1);
    outcomeList = zeros(nTrial, 1);
    outcomeArray = {'goCorrectTarget', 'goCorrectDistractor'};
    for iOutcomeIndex = 1 : length(outcomeArray)
        iOutcome = outcomeArray{iOutcomeIndex};        
        outcomeList = outcomeList + strcmp(trialData.trialOutcome, iOutcome);
    end
    
    trialLogical = trialLogical & outcomeList;
    trialLogical = trialLogical & ismember(trialData.targ1CheckerProp, targ1PropArray(iPropIndex));
    trialLogical = trialLogical & isnan(trialData.ssd);
    goCorrect = find(trialLogical);
    
    
    targAngle = trialData.targAngle;
    % Correct to the right
    trialLogical = ones(nTrial, 1);
    outcomeList = zeros(nTrial, 1);
    outcomeArray = {'goCorrectTarget'};
    for iOutcomeIndex = 1 : length(outcomeArray)
        iOutcome = outcomeArray{iOutcomeIndex};        
        outcomeList = outcomeList + strcmp(trialData.trialOutcome, iOutcome);
    end
    
    trialLogical = trialLogical & outcomeList;
    trialLogical = trialLogical & ismember(trialData.targ1CheckerProp, targ1PropArray(iPropIndex));
    trialLogical = trialLogical & isnan(trialData.ssd);

    targetTrial = ((targAngle > 270) & (targAngle <= 360)) | ...
        ((targAngle >= 0) & (targAngle < 90)) | ...
        ((targAngle > -90) & (targAngle < 0)) | ...
        ((targAngle >= -360) & (targAngle < -270));
    goRightTarget = trialLogical & targetTrial;
    
    % Error to the right
    trialLogical = ones(nTrial, 1);
    outcomeList = zeros(nTrial, 1);
    outcomeArray = {'goCorrectDistractor'};
    for iOutcomeIndex = 1 : length(outcomeArray)
        iOutcome = outcomeArray{iOutcomeIndex};        
        outcomeList = outcomeList + strcmp(trialData.trialOutcome, iOutcome);
    end
    
    trialLogical = trialLogical & outcomeList;
    trialLogical = trialLogical & ismember(trialData.targ1CheckerProp, targ1PropArray(iPropIndex));
    trialLogical = trialLogical & isnan(trialData.ssd);

    distractorTrial = ((targAngle > 90) & (targAngle <= 270)) | ...
        ((targAngle < -90) & (targAngle > -270));
    goRightDistractor = trialLogical & distractorTrial;

    
    goCorrectRight = find(goRightTarget | goRightDistractor);
    
    goCorrectProbRight(iPropIndex) = length(goCorrectRight) / length(goCorrect) ;
end




% % Get incorrect stop probabilities rightward
% if ~isempty(ssdArray)
%     stopIncorrectProbRight = zeros(length(targ1PropArray), 1);
%     for iPropIndex = 1 : length(targ1PropArray);
%         iPercent = targ1PropArray(iPropIndex) * 100;
%
%         % All stop incorrect trials
%         stopStopOutcome = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectDistractor', 'distractorHoldAbort'};
%         stopIncorrectTrial = ccm_trial_selection(monkey, sessionID,  stopStopOutcome, iPercent, 'all', 'all');
%
%         % All stop incorrect rightward trials
%         stopTargetRight = ccm_trial_selection(monkey, sessionID,  {'stopIncorrectTarget', 'targetHoldAbort'}, iPercent, 'all', 'right');
%         stopDistractorRight = ccm_trial_selection(monkey, sessionID,  {'stopIncorrectDistractor', 'distractorHoldAbort'}, iPercent, 'all', 'left');
%         stopIncorrectTrialRight = union(stopTargetRight, stopDistractorRight);
%
%         stopIncorrectProbRight(iPropIndex) = length(stopIncorrectTrialRight) / length(stopIncorrectTrial) ;
%     end
% end


if plotFlag
    choicePlotXMargin = .05;
    goCorrectColor = 'k';
    plot(targ1PropArray, goCorrectProbRight, '-o', 'color', goCorrectColor, 'linewidth', 2, 'markerfacecolor', goCorrectColor, 'markeredgecolor', goCorrectColor)
    hold on
    %     plot(ax(pRightvProbRight), targ1PropArray, goCorrectProbRight, '-o', 'color', goCorrectColor, 'linewidth', 2, 'markerfacecolor', goCorrectColor, 'markeredgecolor', goCorrectColor)
    %     if ~isempty(ssdArray)
    %         plot(ax(pRightvProbRight), targ1PropArray, stopIncorrectProbRight, '-o', 'color', stopIncorrectColor, 'linewidth', 2, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
    %     end
    %     set(ax(pRightvProbRight), 'xtick', targ1PropArray)
    %     set(ax(pRightvProbRight), 'xtickLabel', targ1PropArray*100)
    %     set(get(ax(pRightvProbRight), 'ylabel'), 'String', 'p(Right)')
    %     set(ax(pRightvProbRight),'XLim',[targ1PropArray(1) - choicePlotXMargin targ1PropArray(end) + choicePlotXMargin])
    %     plot(ax(pRightvProbRight), [.5 .5], ylim, '--k')
    set(gca, 'xtick', targ1PropArray)
    set(gca, 'xtickLabel', targ1PropArray*100)
    set(get(gca, 'ylabel'), 'String', 'p(Right)')
    set(gca, 'XLim',[targ1PropArray(1) - choicePlotXMargin targ1PropArray(end) + choicePlotXMargin])
    plot([.5 .5], ylim, '--k')
end

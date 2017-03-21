function [Data] = ccm_session_behavior(subjectID, sessionID, options)
%%

% Set default options or return a default options structure
if nargin < 3
    options.collapseSignal   	= false;
    options.collapseTarg        = false;
    options.include50           = false;
    options.doStops             = true;
    options.USE_TWO_COLORS      = false;
    
    options.plotFlag            = true;
    options.printPlot           = true;
    options.figureHandle      	= 400;
    
    % Return just the default options struct if no input
    if nargin == 0
        Data           = options;
        return
    end
end
plotFlag = options.plotFlag;
printPlot = options.printPlot;
collapseTarg = options.collapseTarg;

% specifiy minimal variables to load:
variables = [ccm_min_vars, 'trialOnset', 'trialDuration','rewardDuration'];


% Load the data
[trialData, SD, ExtraVar] = load_data(subjectID, sessionID, variables);
ssdArray = ExtraVar.ssdArray;
pSignalArray = unique(trialData.targ1CheckerProp);
targAngleArray = unique(trialData.targAngle);
nTrial = size(trialData.trialOutcome, 1);

if ~strcmp(SD.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end


figureHandle = 9239;
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
        
    end


% axes names
pRightVsPCorrect = 7;
rewardRateMA = 7;









if sum(~isnan(trialData.ssd)) == 0
    options.doStops = false;
end


if options.doStops
    
    % ***********************************************************************
    % Inhibition Function:
    %       &
    % SSD vs. Proportion of Response trials
    %       &
    % SSD vs. Proportion(Correct Choice)
    % ***********************************************************************
    
    
    disp('Stopping data')
    optInh                  = ccm_options;
    optInh.collapseTarg     = collapseTarg;
    optInh.printPlot        = printPlot;
    optInh.plotFlag         = plotFlag;
    optInh.figureHandle     = figureHandle;
    optInh.USE_TWO_COLORS      = options.USE_TWO_COLORS;

    dataInh                 = ccm_inhibition(subjectID, sessionID, optInh);
    
    stopRespondProb     = dataInh.stopRespondProb;
    nStop               = dataInh.nStop;
    nStopStop           = dataInh.nStopStop;
    nStopTarg           = dataInh.nStopTarg;
    nStopDist           = dataInh.nStopDist;
    stopTargetProb      = dataInh.stopTargetProb;
    inhibitionFn        = dataInh.inhibitionFn;
    ssrtGrand          	= dataInh.ssrtGrand;
    ssrtMean          	= dataInh.ssrtMean;
    ssrtIntegration   	= dataInh.ssrtIntegration;
    ssrtIntegrationWeighted 	= dataInh.ssrtIntegrationWeighted;
    ssrtIntegrationSimple   	= dataInh.ssrtIntegrationSimple;
    stopRespondProbGrand        = dataInh.stopRespondProbGrand;
    inhibitionFnGrand           = dataInh.inhibitionFnGrand;
    ssrtCollapseGrand           = dataInh.ssrtCollapseGrand;
    ssrtCollapseIntegrationWeighted 	= dataInh.ssrtCollapseIntegrationWeighted;
    ssrtCollapseIntegration             = dataInh.ssrtCollapseIntegration;
    ssrtCollapseMean                    = dataInh.ssrtCollapseMean;
    
end










% % ***********************************************************************
% Reward Rate during the session
% % ***********************************************************************
sessionDuration = floor((trialData.trialOnset(end) + trialData.trialDuration(end)) / 1000/60);
rewardRate = nan(sessionDuration, 1);

rewardDuration = cellfun(@sum, trialData.rewardDuration);

maMinute = 5;

for i = maMinute : sessionDuration
    iMsBegin = (i - maMinute) * 60000;
earliestTrial = find(trialData.trialOnset >= iMsBegin, 1);
    iMsEnd = i * 60000;
latestTrial = find(trialData.trialOnset + trialData.trialDuration <= iMsEnd, 1, 'last');

rewardRate(i) = sum(rewardDuration(earliestTrial : latestTrial)) / maMinute;

end

if plotFlag
    % p(rightward response) vs p(correct choice)
    ax(rewardRateMA) = axes('units', 'centimeters', 'position', [xAxesPosition(3,1) yAxesPosition(3,1) axisWidth axisHeight]);
    hold(ax(rewardRateMA), 'on')
    plot(ax(rewardRateMA), rewardRate, 'color', 'k', 'linewidth', 2)

%     set(ax(rewardRateMA), 'xtick', pSignalArray)
%     set(ax(rewardRateMA), 'xtickLabel', pSignalArray*100)
    set(get(ax(rewardRateMA), 'ylabel'), 'String', 'Reward (ms/min')
    xlim(ax(rewardRateMA), [0 length(rewardRate)])
    set(get(ax(rewardRateMA), 'xlabel'), 'String', 'Minutes')

end % plotFlag










%
%
% % ***********************************************************************
% % Probability(Rightward response) vs Proportion(Correct Choice)
% % ***********************************************************************
goTargTotal      = cell(length(pSignalArray), 1);
goDistTotal      = cell(length(pSignalArray), 1);
nGoTargTotal      = nan(length(pSignalArray), 1);
nGoDistTotal      = nan(length(pSignalArray), 1);
% Get default trial selection options
selectOpt = ccm_trial_selection;

for iPropIndex = 1 : length(pSignalArray);
    iPct = pSignalArray(iPropIndex) * 100;
    selectOpt.rightCheckerPct = iPct;
    selectOpt.plotFlag    = plotFlag;
    selectOpt.figureHandle    = figureHandle;
    
    selectOpt.ssd = 'none';
    selectOpt.outcome     = {'goCorrectDistractor'};
    goDist = ccm_trial_selection(trialData, selectOpt);
    iGoDistIndices = zeros(nTrial, 1);
    iGoDistIndices(goDist) = 1;
    
    oddData = find(isnan(trialData.saccToTargIndex) & iGoDistIndices & ismember(trialData.targ1CheckerProp, pSignalArray));
    if oddData
        fprintf('%d trials are listed as %s but don''t have valid saccades to target:\n', length(oddData), selectOpt.outcome{:})
        %         disp(oddData)
        %         trialData(goDist,:)
    end
    iGoDistIndices(oddData) = 0;
    
    if sum(iGoDistIndices)
        iGoDist = find(iGoDistIndices);
        goDistTotal{iPropIndex} = iGoDist;  % Keep track of totals for grand inhibition fnct
        nGoDistTotal(iPropIndex) = length(iGoDist);  % Keep track of totals for grand inhibition fnct
    end
    
    
    selectOpt.outcome     = {'goCorrectTarget'};
    goTarg = ccm_trial_selection(trialData, selectOpt);
    iGoTargIndices = zeros(nTrial, 1);
    iGoTargIndices(goTarg) = 1;
    
    oddData = find(isnan(trialData.saccToTargIndex) & iGoTargIndices & ismember(trialData.targ1CheckerProp, pSignalArray));
    if oddData
        fprintf('%d trials are listed as %s but don''t have valid saccades to target:\n', length(oddData), selectOpt.outcome{:})
        %         disp(oddData)
        %         trialData(goTarg,:)
    end
    iGoTargIndices(oddData) = 0;
    
    if sum(iGoTargIndices)
        iGoTarg = find(iGoTargIndices);
        goTargTotal{iPropIndex} = iGoTarg;  % Keep track of totals for grand inhibition fnct
        nGoTargTotal(iPropIndex) = length(iGoTarg);  % Keep track of totals for grand inhibition fnct
    end
end

goTargetProb = nGoTargTotal ./ (nGoTargTotal + nGoDistTotal);
goTargetProb(isnan(goTargetProb)) = 1;



% if plotFlag
%     % p(rightward response) vs p(correct choice)
%     ax(pRightVsPCorrect) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 3) yAxesPosition(1, 3) axisWidth axisHeight]);
%     hold(ax(pRightVsPCorrect), 'on')
%     plot(ax(pRightVsPCorrect), pSignalArray, goTargetProb, 'color', 'k', 'linewidth', 2)
%
%     set(ax(pRightVsPCorrect), 'xtick', pSignalArray)
%     set(ax(pRightVsPCorrect), 'xtickLabel', pSignalArray*100)
%     set(get(ax(pRightVsPCorrect), 'ylabel'), 'String', 'p(Correct)')
%     set(ax(pRightVsPCorrect),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
%     set(ax(pRightVsPCorrect),'YLim',[.45 1])
%     plot(ax(pRightVsPCorrect), [.5 .5], ylim, '--k')
%
%     minColorGun = .25;
%     maxColorGun = 1;
%     if ~isempty(ssdArray) && DO_STOPS
%         for jSSDIndex = 1 : length(ssdArray)
%             % Determine color to use for plot based on which checkerboard color
%             % proportion being used. Normalize the available color spectrum to do
%             % it
%             ssdFrac = ssdArray(jSSDIndex) / ssdArray(end);
%             colorGun = minColorGun + (maxColorGun - minColorGun) * ssdFrac;
%             ssdColor = [colorGun 0 0];
%
%             plot(ax(pRightVsPCorrect), pSignalArray, stopTargetProb(:, jSSDIndex), 'color', ssdColor, 'linewidth', 2)
%         end
%     end
%
% end % plotFlag
%









% ***********************************************************************
% Psychometric Function: Proportion(Red Checker) vs Probability(go Right)
% ***********************************************************************
disp('Psychometric data')
optPsy = ccm_psychometric;
optPsy.collapseTarg = collapseTarg;
optPsy.printPlot = printPlot;
optPsy.doStops = options.doStops;
optPsy.USE_TWO_COLORS = options.USE_TWO_COLORS;

optPsy.plotFlag    = plotFlag;
optPsy.figureHandle    = figureHandle;
dataPsych = ccm_psychometric(subjectID, sessionID, optPsy);
% [nGo, nGoRight, nStopIncorrect, nStopIncorrectRight, goRightLogical, goRightSignalStrength, stopRightLogical, stopRightSignalStrength] = ccm_psychometric(subjectID, sessionID, plotFlag);
nGo                 = dataPsych.nGo;
nGoRight            = dataPsych.nGoRight;
nStopIncorrect      = dataPsych.nStopIncorrect;
nStopIncorrectRight = dataPsych.nStopIncorrectRight;
goRightLogical      = dataPsych.goRightLogical;
goRightSignalStrength = dataPsych.goRightSignalStrength;
stopRightLogical    = dataPsych.stopRightLogical;
stopRightSignalStrength = dataPsych.stopRightSignalStrength;
goPsychFn           =   dataPsych.goPsychFn;
%















% ***********************************************************************
% Chronometric Function:     Proportion(Right) vs RT
% ***********************************************************************
% [signalStrengthLeft, signalStrengthRight ...
%     goLeftToTarg, goRightToTarg, goLeftToDist, goRightToDist ...
%     stopLeftToTarg, stopRightToTarg, stopLeftToDist, stopRightToDist] ...
disp('Chronometric data')
optChr = ccm_chronometric;
optChr.collapseTarg = collapseTarg;
optChr.printPlot = printPlot;
optChr.plotFlag    = plotFlag;
optChr.figureHandle    = figureHandle;
optChr.doStops = options.doStops;
optChr.USE_TWO_COLORS = options.USE_TWO_COLORS;

dataChron = ccm_chronometric(subjectID, sessionID, optChr);

pSignalArrayLeft 	= dataChron.pSignalArrayLeft;
pSignalArrayRight	= dataChron.pSignalArrayRight;
goLeftToTarg      	= dataChron.goLeftToTarg;
goRightToTarg     	= dataChron.goRightToTarg;
goLeftToDist       	= dataChron.goLeftToDist;
goRightToDist      	= dataChron.goRightToDist;
stopLeftToTarg    	= dataChron.stopLeftToTarg;
stopRightToTarg    	= dataChron.stopRightToTarg;
stopLeftToDist     	= dataChron.stopLeftToDist;
stopRightToDist    	= dataChron.stopRightToDist;








for i = 1 : size(dataChron, 1)
    
    Data(i).pSignalArray    = pSignalArray;
    Data(i).ssdArray        = ssdArray;
    
    if options.doStops
    Data(i).stopRespondProb      = dataInh(i).stopRespondProb;
    Data(i).nStop                = dataInh(i).nStop;
    Data(i).nStopStop            = dataInh(i).nStopStop;
    Data(i).nStopTarg            = dataInh(i).nStopTarg;
    Data(i).nStopDist            = dataInh(i).nStopDist;
    % Data(i).goTargetProb         = dataInh(i).goTargetProb;
    Data(i).stopTargetProb       = dataInh(i).stopTargetProb;
    Data(i).inhibitionFn         = dataInh(i).inhibitionFn;
    Data(i).ssrtGrand            = dataInh(i).ssrtGrand;
    Data(i).ssrtMean             = dataInh(i).ssrtMean;
    Data(i).ssrtIntegration      = dataInh(i).ssrtIntegration;
    Data(i).ssrtIntegrationWeighted = dataInh(i).ssrtIntegrationWeighted;
    Data(i).ssrtIntegrationSimple = dataInh(i).ssrtIntegrationSimple;
    Data(i).ssrtCollapseGrand   	= dataInh(i).ssrtCollapseGrand;
    Data(i).ssrtCollapseIntegrationWeighted     	= dataInh(i).ssrtCollapseIntegrationWeighted;
    Data(i).ssrtCollapseIntegration  	= dataInh(i).ssrtCollapseIntegration;
    Data(i).ssrtCollapseMean    	= dataInh(i).ssrtCollapseMean;
    Data(i).stopRespondProbGrand = dataInh(i).stopRespondProbGrand;
    Data(i).inhibitionFnGrand    = dataInh(i).inhibitionFnGrand;
    end
    
    Data(i).goPsychFn         	= dataPsych(i).goPsychFn;
    Data(i).nGo                  = dataPsych(i).nGo;
    Data(i).nGoRight             = dataPsych(i).nGoRight;
    Data(i).nStopIncorrect       = dataPsych(i).nStopIncorrect;
    Data(i).nStopIncorrectRight  = dataPsych(i).nStopIncorrectRight;
    Data(i).goRightLogical       = dataPsych(i).goRightLogical;
    Data(i).goRightSignalStrength = dataPsych(i).goRightSignalStrength;
    Data(i).stopRightLogical     = dataPsych(i).stopRightLogical;
    Data(i).stopRightSignalStrength = dataPsych(i).stopRightSignalStrength;
    
    
    Data(i).pSignalArrayLeft = dataChron(i).pSignalArrayLeft;
    Data(i).pSignalArrayRight = dataChron(i).pSignalArrayRight;
    Data(i).goLeftToTarg    = dataChron(i).goLeftToTarg;
    Data(i).goRightToTarg   = dataChron(i).goRightToTarg;
    Data(i).goLeftToDist    = dataChron(i).goLeftToDist;
    Data(i).goRightToDist   = dataChron(i).goRightToDist;
    Data(i).stopLeftToTarg  = dataChron(i).stopLeftToTarg;
    Data(i).stopRightToTarg = dataChron(i).stopRightToTarg;
    Data(i).stopLeftToDist  = dataChron(i).stopLeftToDist;
    Data(i).stopRightToDist = dataChron(i).stopRightToDist;
end


      h=axes('Position', [0 0 1 1], 'Visible', 'Off');
      titleString = sprintf('\n\n%s \t %s', subjectID,sessionID);
      text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')

if printPlot && plotFlag   
    print(figureHandle,fullfile(local_figure_path, subjectID, [sessionID,'_ccm_behavior.pdf']),'-dpdf', '-r300')
end
% delete(localDataFile);



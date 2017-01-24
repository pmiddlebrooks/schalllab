function data = ccm_latency_probability_function(subjectID, sessionID, plotFlag)



DO_STOPS = 1;


% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
pSignalArray = ExtraVar.pSignalArray;
ssdArray = ExtraVar.ssdArray;



if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end


% Set up figure
if plotFlag
    nPlotColumn = 2;
    % Plot as many rows as there are discriminatory levels (use the maximum
    % from the individual targets), plus an extra row for collapsed RTs
    nPlotRow = 3;
    figureHandleLP = 1003;
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition]       = standard_figure(nPlotRow, nPlotColumn, 'portrait', figureHandleLP);
end




nTrial = size(trialData.trialOutcome, 1);
goColor = [0 0 0];
stopColor = [1 0 0];




plotTitle = {'Target Right', 'Target Left', 'All Trials', 'Left/Right Collapsed'};
pGoTarg        = nan(4, length(pSignalArray));
pGoDist    = nan(4, length(pSignalArray));
pStopTarg      = nan(4, length(pSignalArray));
pStopDist  = nan(4, length(pSignalArray));
rtGoTarg          = cell(4, length(pSignalArray));
rtGoDist      = cell(4, length(pSignalArray));
rtStopTarg        = cell(4, length(pSignalArray));
rtStopDist    = cell(4, length(pSignalArray));


% Make LP functions for each target separately
for iTarg = 1 : 4
    switch iTarg
        case 1
            iPlotRow = 1;
            iPropRange = pSignalArray(pSignalArray >= .5);
            iPropRange = flipud(iPropRange); % Rearrage order from easiest to hardest
            iPlotColumn = 2;
        case 2
            iPlotRow = 1;
            iPropRange = pSignalArray(pSignalArray <= .5);
            iPlotColumn = 1;
        case 3  % all of the signal strengths together
            iPlotRow = 3;
            iPlotColumn = 2;
            iPropRight = pSignalArray(pSignalArray >= .51);
            iPropRight = flipud(iPropRight); % Rearrage order from easiest to hardest
            iPropLeft = pSignalArray(pSignalArray <= .49);
            if ~ismember(.5, pSignalArray)
                iPropRange = [iPropLeft; iPropRight];
            elseif ismember(.5, pSignalArray)
                iPropRange = [iPropLeft; [.5; .5]; iPropRight];
            end
        case 4  % binned difficulty levels
            iPlotRow = 2;
            iPlotColumn = 1;
            iPropRight = pSignalArray(pSignalArray >= .51);
            iPropRight = flipud(iPropRight); % Rearrage order from easiest to hardest
            iPropLeft = pSignalArray(pSignalArray <= .49);
            iPropRange = [iPropLeft, iPropRight];
            if ismember(.5, pSignalArray)
                iPropRange = [iPropRange; [.5 .5]];
            end
    end
    
    
        % Get default trial selection options
        selectOpt       = ccm_trial_selection;

        for jPropIndex = 1 : size(iPropRange, 1)
        jPct = iPropRange(jPropIndex,:) * 100;
        selectOpt.rightCheckerPct = jPct;
        
        
        % Get correct go RT distribution
        
        selectOpt.ssd       = 'none';
        
        
        selectOpt.outcome     = {'goCorrectTarget', 'targetHoldAbort'};
        if iTarg ~= 3 || (iTarg == 3 && jPct ~= 50)
            selectOpt.targDir    	= 'collapse';
            ijGoTargetTrial = ccm_trial_selection(trialData, selectOpt);
            % If it'ss the first .5 signal strength, get only the left target trials
        elseif iTarg == 3 && jPct  == 50 && iPropRange(jPropIndex+1,:) * 100 == 50
            selectOpt.targDir    	= 'left';
            ijGoTargetTrial = ccm_trial_selection(trialData, selectOpt);
            % If it'ss the first .5 signal strength, get only the right target trials
        elseif iTarg == 3 && jPct  == 50 && iPropRange(jPropIndex+1,:) * 100 ~= 50
            selectOpt.targDir    	= 'right';
            ijGoTargetTrial = ccm_trial_selection(trialData, selectOpt);
        end
        ijGoTargetTrialIndices = zeros(nTrial, 1);
        ijGoTargetTrialIndices(ijGoTargetTrial) = 1;
        % Get rid of any trials in which the response to target index was
        % greater than the available reponses from a given trial
        mismatch = (trialData.saccToTargIndex - cellfun(@length, trialData.saccBegin)) > 0;
        ijGoTargetTrialIndices(mismatch) = 0;
        
        
        selectOpt.outcome     = {'goCorrectDistractor', 'distractorHoldAbort'};
        if iTarg ~= 3 || (iTarg == 3 && jPct ~= 50)
            selectOpt.targDir    	= 'collapse';
            ijGoDistractorTrial = ccm_trial_selection(trialData, selectOpt);
            % If it'ss the first .5 signal strength, get only the left target trials
        elseif iTarg == 3 && jPct  == 50 && iPropRange(jPropIndex+1,:) * 100 == 50
            selectOpt.targDir    	= 'left';
            ijGoDistractorTrial = ccm_trial_selection(trialData, selectOpt);
            % If it'ss the first .5 signal strength, get only the right target trials
        elseif iTarg == 3 && jPct  == 50 && iPropRange(jPropIndex+1,:) * 100 ~= 50
            selectOpt.targDir    	= 'right';
            ijGoDistractorTrial = ccm_trial_selection(trialData, selectOpt);
        end
        ijGoDistractorTrialIndices = zeros(nTrial, 1);
        ijGoDistractorTrialIndices(ijGoDistractorTrial) = 1;
        % Get rid of any trials in which the response to target index was
        % greater than the available reponses from a given trial
        ijGoDistractorTrialIndices(mismatch) = 0;
        
        
        pGoTarg(iTarg, jPropIndex) = sum(ijGoTargetTrialIndices) / (sum(ijGoTargetTrialIndices) + sum(ijGoDistractorTrialIndices));
        pGoDist(iTarg, jPropIndex) = 1 - pGoTarg(iTarg, jPropIndex);
        
        
        if sum(ijGoTargetTrialIndices)
%             rtGoTarg{iTarg, jPropIndex} = trialData.responseOnset(find(ijGoTargetTrialIndices)) - trialData.responseOnset(find(ijGoTargetTrialIndices));
%             rtGoDist{iTarg, jPropIndex} = trialData.responseOnset(find(ijGoDistractorTrialIndices)) - trialData.responseOnset(find(ijGoDistractorTrialIndices));
            rtGoTarg{iTarg, jPropIndex} = trialData.rt(find(ijGoTargetTrialIndices));
            rtGoDist{iTarg, jPropIndex} = trialData.rt(find(ijGoDistractorTrialIndices));
        end
        
        
        
        
        
        
        
        
        if DO_STOPS
            % Get incorrect stop RT distribution
            
            % Get default trial selection options
            selectOpt       = ccm_trial_selection;
            selectOpt.ssd       = 'collapse';
            
            selectOpt.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
            if iTarg ~= 3 || (iTarg == 3 && jPct ~= 50)
                selectOpt.targDir       = 'collapse';
                ijStopTargetTrial = ccm_trial_selection(trialData, selectOpt);
                % If it'ss the first .5 signal strength, get only the left target trials
            elseif iTarg == 3 && jPct  == 50 && iPropRange(jPropIndex+1,:) * 100 == 50
                selectOpt.targDir       = 'left';
                ijStopTargetTrial = ccm_trial_selection(trialData, selectOpt);
                % If it'ss the first .5 signal strength, get only the right target trials
            elseif iTarg == 3 && jPct  == 50 && iPropRange(jPropIndex+1,:) * 100 ~= 50
                selectOpt.targDir       = 'right';
                ijStopTargetTrial = ccm_trial_selection(trialData, selectOpt);
            end
            ijStopTargetTrialIndices = zeros(nTrial, 1);
            ijStopTargetTrialIndices(ijStopTargetTrial) = 1;
            % Get rid of any trials in which the response to target index was
            % greater than the available reponses from a given trial
            mismatch = (trialData.saccToTargIndex - cellfun(@length, trialData.saccBegin)) > 0;
            ijStopTargetTrialIndices(mismatch) = 0;
            
            
            selectOpt.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
            if iTarg ~= 3 || (iTarg == 3 && jPct ~= 50)
                selectOpt.targDir       = 'collapse';
                ijStopDistractorTrial = ccm_trial_selection(trialData, selectOpt);
                % If it'ss the first .5 signal strength, get only the left target trials
            elseif iTarg == 3 && jPct  == 50 && iPropRange(jPropIndex+1,:) * 100 == 50
                selectOpt.targDir       = 'left';
                ijStopDistractorTrial = ccm_trial_selection(trialData, selectOpt);
                % If it'ss the first .5 signal strength, get only the right target trials
            elseif iTarg == 3 && jPct  == 50 && iPropRange(jPropIndex+1,:) * 100 ~= 50
                selectOpt.targDir       = 'right';
                ijStopDistractorTrial = ccm_trial_selection(trialData, selectOpt);
            end
            ijStopDistractorTrialIndices = zeros(nTrial, 1);
            ijStopDistractorTrialIndices(ijStopDistractorTrial) = 1;
            % Get rid of any trials in which the response to target index was
            % greater than the available reponses from a given trial
            ijStopDistractorTrialIndices(mismatch) = 0;
            
            
            pStopTarg(iTarg, jPropIndex) = sum(ijStopTargetTrialIndices) / (sum(ijStopTargetTrialIndices) + sum(ijStopDistractorTrialIndices));
            pStopDist(iTarg, jPropIndex) = 1 - pStopTarg(iTarg, jPropIndex);
            
            if sum(ijStopTargetTrialIndices)
%                 rtStopTarg{iTarg, jPropIndex} = trialData.responseOnset(find(ijStopTargetTrialIndices)) - trialData.responseOnset(find(ijStopTargetTrialIndices));
%                 rtStopDist{iTarg, jPropIndex} = trialData.responseOnset(find(ijStopDistractorTrialIndices)) - trialData.responseOnset(find(ijStopDistractorTrialIndices));
                rtStopTarg{iTarg, jPropIndex} = trialData.rt(find(ijStopTargetTrialIndices));
                rtStopDist{iTarg, jPropIndex} = trialData.rt(find(ijStopDistractorTrialIndices));
            end
        end
        
        
        
        
        
    end % jProp
    
    
    rtGoTargMean{iTarg}         = cellfun(@mean, rtGoTarg(iTarg, :));
    rtGoDistMean{iTarg}     = cellfun(@mean, rtGoDist(iTarg, :));
    rtStopTargMean{iTarg}       = cellfun(@mean, rtStopTarg(iTarg, :));
    rtStopDistMean{iTarg}   = cellfun(@mean, rtStopDist(iTarg, :));
    
    rtGoTargStd{iTarg}         = cellfun(@std, rtGoTarg(iTarg, :));
    rtGoDistStd{iTarg}     = cellfun(@std, rtGoDist(iTarg, :));
    rtStopTargStd{iTarg}       = cellfun(@std, rtStopTarg(iTarg, :));
    rtStopDistStd{iTarg}   = cellfun(@std, rtStopDist(iTarg, :));
    
    
    
    if plotFlag
        % RT Distribution
        figure(figureHandleLP)
        if iTarg < 3
            ax(iPlotRow, iPlotColumn) = axes('units', 'centimeters', 'position', [xAxesPosition(iPlotRow, iPlotColumn) yAxesPosition(iPlotRow, iPlotColumn) axisWidth axisHeight]);
        elseif iTarg > 2
            ax(iPlotRow, iPlotColumn) = axes('units', 'centimeters', 'position', [mean(xAxesPosition(iPlotRow, :)) yAxesPosition(iPlotRow, iPlotColumn) axisWidth axisHeight]);
        end
        hold(ax(iPlotRow, iPlotColumn), 'on')
        xlim(ax(iPlotRow, iPlotColumn), [0 1])
        title(ax(iPlotRow, iPlotColumn), plotTitle{iTarg});
        %             set(ax(iPlotRow, iPlotColumn), 'xlim', [0 1])
        %             set(ax(iPlotRow, iPlotColumn), 'ylim', [0 1])
        %             if sum(ijGoTargetTrialIndices) > 1
        markerSize = 10;
        plot(ax(iPlotRow, iPlotColumn), pGoTarg(iTarg, :), rtGoTargMean{iTarg}, 'o', 'markeredgecolor', goColor, 'markerfacecolor', goColor, 'markersize', markerSize)
        plot(ax(iPlotRow, iPlotColumn), pGoDist(iTarg, :), rtGoDistMean{iTarg}, 'o', 'markeredgecolor', goColor, 'markersize', markerSize)
        plot(ax(iPlotRow, iPlotColumn), pStopTarg(iTarg, :), rtStopTargMean{iTarg}, 'o', 'markeredgecolor', stopColor, 'markerfacecolor', stopColor, 'markersize', markerSize)
        plot(ax(iPlotRow, iPlotColumn), pStopDist(iTarg, :), rtStopDistMean{iTarg}, 'o', 'markeredgecolor', stopColor, 'markersize', markerSize)
        
        errorbar(ax(iPlotRow, iPlotColumn), pGoTarg(iTarg, :), rtGoTargMean{iTarg}, rtGoTargStd{iTarg}, 'linestyle' , 'none', 'color', 'k')
        errorbar(ax(iPlotRow, iPlotColumn), pGoDist(iTarg, :), rtGoDistMean{iTarg}, rtGoDistStd{iTarg}, 'linestyle' , 'none', 'color', 'k')
        errorbar(ax(iPlotRow, iPlotColumn), pStopTarg(iTarg, :), rtStopTargMean{iTarg}, rtStopTargStd{iTarg}, 'linestyle' , 'none', 'color', 'k')
        errorbar(ax(iPlotRow, iPlotColumn), pStopDist(iTarg, :), rtStopDistMean{iTarg}, rtStopDistStd{iTarg}, 'linestyle' , 'none', 'color', 'k')
    end
    
    
end % iTarg

data.pGoTarg = pGoTarg;
data.pGoDist = pGoDist;
data.rtGoTarg = rtGoTarg;
data.rtGoDist = rtGoDist;
data.pStopTarg = pStopTarg;
data.pStopDist = pStopDist;
data.rtStopTarg = rtStopTarg;
data.rtStopDist = rtStopDist;

end

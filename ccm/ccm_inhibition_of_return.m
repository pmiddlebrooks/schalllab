function ccm_inhibition_of_return(subjectArray, sessionArray, options)
%
% function ccm_inhibition_of_return(subjectArray, sessionArray, options)
%
% input:
%   subjectArray:
%             a single subject, e.g. {'Broca'}
%             a list of subjects, e.g. {broca', xena'}, etc
%   sessionArray: a cell input
%             a single session (if subjectArray has single subject)e.g. {'bp111n01'}
%             a list of single sessions (if subjectArray has single subject)e.g. {'bp111n01', 'bp112n02'}
%             a single batch session term that refers to a collection of sesssions data file (e.g. {'behavior1'}
%
%   options: A structure with various ways to select/organize data: If
%   ccm_inhibition_of_return.m is called without input arguments, the default
%   options structure is returned. options has the following fields with
%   possible values (default listed first):
%
%    options.dataType = 'neuron', 'lfp', 'erp';
%
%    options.figureHandle   = 1000;
%    options.printPlot      = false, true;
%    options.plotFlag       = true, false;
%    options.collapseSignal = false, true;
%     options.collapseTarg         = false, true;
%    options.doStops        = true, false;
%    options.filterData 	= false, true;
%    options.stopHz         = 50, <any number, above which signal is filtered;
%    options.normalize      = false, true;
%    options.unitArray      = {'spikeUnit17a'},'each', units want to analyze
%
%
%   Include/exclude aborted trials: This will affect the number of paired
%   and triplet trials that make it into analyses, since aborts between
%   trials may or may not count as successive trials. (this seems not to
%   matter though: deleteAborts = true vs. false
%
%   Analyze data across sessions, taking the mean across sessions, or
%   analyze with all data collapsed (as if one big session). Also doesn't
%   alter the results much. acrossSession = true vs. false


if nargin < 3
    %    options.subjectArray     = {'broca','xena'};
    %    options.sessionArray     = {'behavior'};
    options.deleteAborts     = false;
    options.printPlot        = true;
    options.plotFlag         = true;
    options.plotData            = 'both';
    if nargin == 0
        data = options;
        return
    end
end

% subjectArray    = options.subjectArray;
% sessionArray    = options.sessionArray;
plotFlag        = options.plotFlag;
printPlot       = options.printPlot;
plotData        = options.plotData;
deleteAborts    = options.deleteAborts;

figureHandle = 48;

rrColor = 'k';

acrossSession = true;
adjustForMean = false;

firstOutcomeArray = {'goCorrectTarget', 'goCorrectDistractor', 'stopIncorrectTarget'};
% firstOutcomeArray = {'goCorrectTarget'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     LOOP THROUGH SUBJECTS AND PLOT PAIRED IOR AND RTS ON ONE GRAPH EACH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(subjectArray)
    iSubject = subjectArray{i};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    LOAD DATA AND SET VARIABLES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iSubject = subjectArray{i};
    
    switch sessionArray
        case 'concat'
            iFile = fullfile('local_data',iSubject,strcat(iSubject,'RT.mat'));
            load(iFile)  % Loads trialData into workspace (and SessionData)
        case {'behavior1', 'behavior2', 'neural1', 'neural2'}
            iFile = fullfile('local_data',iSubject,strcat(iSubject,'_',sessionArray,'.mat'));
            load(iFile)  % Loads trialData into workspace (and SessionData)
        otherwise
            % Load the data
            [trialData, SessionData, ExtraVar] = load_data(iSubject, sessionArray);
            pSignalArray    = ExtraVar.pSignalArray;
            targAngleArray	= ExtraVar.targAngleArray;
            ssdArray        = ExtraVar.ssdArray;
            nSSD            = length(ssdArray);
            
            if ~strcmp(SessionData.taskID, 'ccm')
                fprintf('Not a chioce countermanding saccade session, try again\n')
                return
            end
    end
    
    % Get rid of extraneous RTs
    [rt, outlierTrial] = truncate_rt(trialData.rt, 100, 1200);
    trialData.rt(outlierTrial) = nan;
    
    
    if deleteAborts
        selectOpt = ccm_options;
        selectOpt.outcome = {...
            'goCorrectTarget', 'goCorrectDistractor', ...
            'stopCorrect', ...
            'targetHoldAbort', 'distractorHoldAbort', ...
            'stopIncorrectTarget', 'stopIncorrectDistractor'};
        validTrial = ccm_trial_selection(trialData, selectOpt);
        trialData = trialData(validTrial,:);
    end
    
    
    % Treat data differently if analyzing across sessions vs collapsed
    % sessions:
    % excludeTrialTriplet: Find session switch trials so we don't process them as if a new
    % session was a continuation of one big session
    if acrossSession
        if ismember(sessionArray, {'concat','behavior1', 'behavior2', 'neural1', 'neural2'})
            nSession = max(trialData.sessionTag);
        else
            nSession = 1;
        end
        excludeTrialTriplet = [];
    else
        nSession = 1;
        excludeTrialPair = find(diff(trialData.sessionTag) < 0);
        excludeTrialTriplet = [excludeTrialPair; excludeTrialPair-1]; % Exclude last 2 trials of a session as possible beginning trials in triplets
    end
    
    
    % Initialize vectors for RTs
    rtNs          = nan(nSession, 1);
    
    % SomeOutcome -> Correct Pairs
    saccRR1     = nan(nSession, 1);
    saccRR2     = nan(nSession, 1);
    saccRRIor 	= nan(nSession, 1);
    saccRL1     = nan(nSession, 1);
    saccRL2     = nan(nSession, 1);
    saccRLIor 	= nan(nSession, 1);
    saccLL1     = nan(nSession, 1);
    saccLL2     = nan(nSession, 1);
    saccLLIor 	= nan(nSession, 1);
    saccLR1     = nan(nSession, 1);
    saccLR2     = nan(nSession, 1);
    saccLRIor 	= nan(nSession, 1);
    
    nStimRR     = nan(nSession, 1);
    nStimRL     = nan(nSession, 1);
    nStimLL     = nan(nSession, 1);
    nStimLR     = nan(nSession, 1);
    
    
    
    for k = 1 : length(firstOutcomeArray)
        
        for j = 1 : nSession
            
            if acrossSession && ismember(sessionArray, {'concat','behavior1', 'behavior2', 'neural1', 'neural2'})
                jTD = trialData(trialData.sessionTag == j, :);
            else
                jTD = trialData;
            end
            
            % Get overall mean rt for correct choices (on no-stop trials)
            opt = ccm_options;
            opt.outcome = {'goCorrectTarget'};%, 'goCorrectDistractor'};
            opt.ssd = 'none';
            nsTrial = ccm_trial_selection(jTD, opt);
            rtNs(j) = nanmean(jTD.rt(nsTrial));
            
            
            
            
            % Initialize structure with 2 levels
            Opt2(1) = opt;
            Opt2(2) = opt;
            
            
            
            % Set the first trial outcome for the pairs:
            Opt2(1).outcome = {firstOutcomeArray{k}};
            if strcmp(firstOutcomeArray{k}, 'stopIncorrectTarget')
                Opt2(1).ssd = 'collapse';
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %      Correct/Error Choices -> Correct Choices
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %--------------------------------------------------------------------------
            % saccR -> saccR
            
            Opt2(1).responseDir = {'right'};
            
            Opt2(2).responseDir = {'right'};
            
            saccRRTrial = ccm_trial_sequence(jTD, Opt2);
            saccRRTrial = setxor(saccRRTrial, excludeTrialTriplet);
            saccRR1(j) = nanmean(jTD.rt(saccRRTrial));
            saccRR2(j) = nanmean(jTD.rt(saccRRTrial + 1));
            saccRRIor(j) = nanmean(saccRR2(j) - saccRR1(j));
            %         [h,p,ci,stats] = ttest2(jTD.rt(saccRRTrial), jTD.rt(saccRRTrial+1));
            
            
            %--------------------------------------------------------------------------
            % saccR -> saccL
            
            Opt2(1).responseDir = {'right'};
            
            Opt2(2).responseDir = {'left'};
            
            saccRLTrial = ccm_trial_sequence(jTD, Opt2);
            saccRLTrial = setxor(saccRLTrial, excludeTrialTriplet);
            saccRL1(j) = nanmean(jTD.rt(saccRLTrial));
            saccRL2(j) = nanmean(jTD.rt(saccRLTrial + 1));
            saccRLIor(j) = nanmean(saccRL2(j) - saccRL1(j));
            %         [h,p,ci,stats] = ttest2(jTD.rt(saccRLTrial), jTD.rt(saccRLTrial+1));
            
            %--------------------------------------------------------------------------
            % saccL -> saccL
            
            Opt2(1).responseDir = {'left'};
            
            Opt2(2).responseDir = {'left'};
            
            saccLLTrial = ccm_trial_sequence(jTD, Opt2);
            saccLLTrial = setxor(saccLLTrial, excludeTrialTriplet);
            saccLL1(j) = nanmean(jTD.rt(saccLLTrial));
            saccLL2(j) = nanmean(jTD.rt(saccLLTrial + 1));
            saccLLIor(j) = nanmean(saccLL2(j) - saccLL1(j));
            %         [h,p,ci,stats] = ttest2(jTD.rt(saccLLTrial), jTD.rt(saccLLTrial+1));
            
            
            %--------------------------------------------------------------------------
            % saccL -> saccR
            
            Opt2(1).responseDir = {'left'};
            
            Opt2(2).responseDir = {'right'};
            
            saccLRTrial = ccm_trial_sequence(jTD, Opt2);
            saccLRTrial = setxor(saccLRTrial, excludeTrialTriplet);
            saccLR1(j) = nanmean(jTD.rt(saccLRTrial));
            saccLR2(j) = nanmean(jTD.rt(saccLRTrial + 1));
            saccLRIor(j) = nanmean(saccLR2(j) - saccLR1(j));
            %         [h,p,ci,stats] = ttest2(jTD.rt(saccLRTrial), jTD.rt(saccLRTrial+1));
            
            
            
            
            
            
            
        end % for j = 1 : nSession
        
        
        
        fprintf('%s paired t-tests\n', firstOutcomeArray{k})
        [h,p,ci,stats] = ttest2(saccRR1, saccRR2);
        fprintf('R->R:\tp = %0.3f\n', p)
        [h,p,ci,stats] = ttest2(saccRL1, saccRL2);
        fprintf('R->L:\tp = %0.3f\n', p)
        [h,p,ci,stats] = ttest2(saccLL1, saccLL2);
        fprintf('L->L:\tp = %0.3f\n', p)
        [h,p,ci,stats] = ttest2(saccLR1, saccLR2);
        fprintf('L->R:\tp = %0.3f\n', p)
        
        
        if plotFlag
            figureHandle = figureHandle + 1;
            % ********************************************************
            %       Correct Choice -> Correct Choice
            % ********************************************************
            
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(2, 1, figureHandle);
            clf
            axRT = axes('units', 'centimeters', 'position', [xAxesPosition(1) yAxesPosition(1) axisWidth axisHeight]);
            hold(axRT, 'all')
            title(sprintf('%s', firstOutcomeArray{k}))
            axDiff = axes('units', 'centimeters', 'position', [xAxesPosition(1) yAxesPosition(2) axisWidth axisHeight/2]);
            hold(axDiff, 'all')
           
            
            % Plot RTs from all the sessions
            % ********************************************************
            
            minRT = min([nanmean(saccRR1), nanmean(saccRR2), nanmean(saccRL1), nanmean(saccRL2), nanmean(saccLL1), nanmean(saccLL2), nanmean(saccLR1), nanmean(saccLR2)]);
            maxRT = max([nanmean(saccRR1), nanmean(saccRR2), nanmean(saccRL1), nanmean(saccRL2), nanmean(saccLL1), nanmean(saccLL2), nanmean(saccLR1), nanmean(saccLR2)]);
            axes(axRT)
            %             set(ax1, 'ylim', [.9 * minRT, 1.1 * maxRT])
            set(axRT, 'ylim', [150 320])
            plot(axRT, [0 6], [nanmean(rtNs), nanmean(rtNs)], '--', 'color', 'k')
            xData = [.85 1.15 1.85 2.15 3.85 4.15 4.85 5.15];
            yData = [nanmean(saccLL1), nanmean(saccLL2), nanmean(saccLR1), nanmean(saccLR2), nanmean(saccRR1), nanmean(saccRR2), nanmean(saccRL1), nanmean(saccRL2)];
            bar(xData, yData, 1)
            
            % Error bars
            errorData = [std(saccLL1), std(saccLL2), std(saccLR1), std(saccLR2), std(saccRR1), std(saccRR2), std(saccRL1), std(saccRL2)];
            errorbar(axRT, xData, yData,  errorData, 'linestyle' , 'none', 'color', 'k', 'linewidth' , 2)
            
            % Plot IORs from all the sessions
            % ********************************************************
            axes(axDiff)
            iMarkerSize = 15;
            set(axDiff, 'ylim', [-120 120])
            plot(axDiff, [0 6], [0 0], '--', 'color', 'k')
            
            xData = [1 2 4 5];
            yData = [nanmean(saccLLIor), nanmean(saccLRIor), nanmean(saccRRIor), nanmean(saccRLIor)];
            plot(axDiff, xData, yData, 'o', 'color', rrColor, 'markerSize', iMarkerSize, 'markerFaceColor', rrColor)
            
            % Error bars
            errorData = [std(saccLLIor), std(saccLRIor), std(saccRRIor), std(saccRLIor)];
            errorbar(axDiff, xData, yData,  errorData, 'linestyle' , 'none', 'color', 'k', 'linewidth' , 2)
            
                print(figureHandle,[local_figure_path, sprintf('%s_IOR_%s_%s', iSubject, firstOutcomeArray{k}, sessionArray)],'-dpdf', '-r300')
           
            
        end % if plotFlag
    end % for k = 1 : length(firstOutcomeArray}
end % for i = 1 : nSubject



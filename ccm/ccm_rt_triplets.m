function data = ccm_rt_triplets(subjectArray, sessionArray, Opt)
%
% function data = ccm_rt_triplets(subjectArray, sessionArray, opt)
%
% Replicating Nelson et al 2010 and Emeric et al. 2007 for Choice
% countermanding data.
% These analyses can be done a few ways:
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
%   Opt: A structure with various ways to select/organize data: If
%   ccm_rt_triplets.m is called without input arguments, the default
%   Opt structure is returned. Opt has the following fields with
%   possible values (default listed first):
%
%
%    Opt.figureHandle   = 1000;
%    Opt.Opt.printPlot      = false, true;
%    Opt.plotFlag       = true, false;
%    Opt.collapseSignal = false, true;
%     Opt.collapseTarg         = false, true;
%    Opt.doStops        = true, false;
%    Opt.filterData 	= false, true;
%    Opt.stopHz         = 50, <any number, above which signal is filtered;
%    Opt.normalize      = false, true;
%    Opt.unitArray      = {'spikeUnit17a'},'each', units want to analyze


%   Include/exclude aborted trials: This will affect the number of paired
%   and triplet trials that make it into analyses, since aborts between
%   trials may or may not count as successive trials. (this seems not to
%   matter though: Opt.deleteAborts = true vs. false
%
%   Analyze data across sessions, taking the mean across sessions, or
%   analyze with all data collapsed (as if one big session). Also doesn't
%   alter the results much. acrossSession = true vs. false



% %% Analyze RTs surrounding fixation aborts
%
% load('local_data/broca/brocaRT.mat')
%
% % Fixation aborts only
% fa = strcmp(trialData.trialOutcome, 'fixationAbort');
% prefa = [fa(2:end); false];
% postfa = [false; fa(1:end-1)];
% nanmean(trialData.rt(prefa))
% nanmean(trialData.rt(postfa))

% %% Analyze RTs surrounding any/all aborts
%
% ab = strcmp(trialData.trialOutcome, ('fixationAbort')) | ...
%     strcmp(trialData.trialOutcome, ('choiceStimulusAbort')) | ...
%     strcmp(trialData.trialOutcome, ('noFixation')) | ...
%     strcmp(trialData.trialOutcome, ('saccadeAbort')) | ...
%     strcmp(trialData.trialOutcome, ('targetHoldAbort')) | ...
%     strcmp(trialData.trialOutcome, ('distractorHoldAbort'));
%
% prefa = [ab(2:end); false];
% postfa = [false; ab(1:end-1)];
% nanmean(trialData.rt(prefa))
% nanmean(trialData.rt(postfa))

% Replicating Nelson et al 2010 and Emeric et al. 2007 for Choice
% countermanding data.
% These analyses can be done a few ways:
%
%   Include/exclude aborted trials: This will affect the number of paired
%   and triplet trials that make it into analyses, since aborts between
%   trials may or may not count as successive trials. (this seems not to
%   matter though: Opt.deleteAborts = true vs. false
%
%   Analyze data across sessions, taking the mean across sessions, or
%   analyze with all data collapsed (as if one big session). Also doesn't
%   alter the results much. acrossSession = true vs. false

%% triplet analysis:

% Keeping aborted trials in
% without respect to choice difficulty. As a first, dont' remove any aborted
% trials. This will greatly reduce the data, but is a more valid test

% Get default Options
if nargin < 3
    Opt = ccm_options;
    Opt.deleteAborts     = true;
end

acrossSession = true;
adjustForMean = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Opt.plotFlag
nsColor     = [0 .8 0];
cColor      = [1 0 .5];
ncColor     = [150 50 0] ./ 255;
eColor      = 'b';

figureHandle = 466;
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(1, 2, figureHandle);
clf
axPair = axes('units', 'centimeters', 'position', [xAxesPosition(1,1) yAxesPosition(1,1) axisWidth*16/20 axisHeight * .4]);
hold(axPair, 'all')
axTrip = axes('units', 'centimeters', 'position', [xAxesPosition(1,2)*.85 yAxesPosition(1,2) axisWidth*24/20 axisHeight * .4]);
hold(axTrip, 'all')
set(axTrip, 'yTickLabel', [])
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     LOOP THROUGH SUBJECTS AND PLOT TRIPLET RTS ON ONE GRAPH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(subjectArray)
    
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
    
    
    
    if Opt.deleteAborts
        selectOpt = ccm_trial_selection;
        selectOpt.outcome = {...
            'goCorrectTarget', 'goCorrectDistractor', ...
            'stopCorrect', ...
            'targetHoldAbort', 'distractorHoldAbort', ...
            'stopIncorrectTarget', 'stopIncorrectDistractor'};
        validTrial = ccm_trial_selection(trialData, selectOpt);
        trialData = trialData(validTrial,:);
    end
    
    % Find session switch trials so we don't process them as if a new
    % session was a continuation of one big session
    %     if strcmp(iSubject, 'human')
    %         excludeTrialPair = find(diff(trialData.sessionTag) < 0);
    %     else
    %         excludeTrialPair = find(diff(trialData.trial) < 0);
    %     end
    %     excludeTrialTriplet = [excludeTrialPair; excludeTrialPair-1]; % Exclude last 2 trials of a session as possible beginning trials in triplets
    
    
    
    
    % Treat data differently if analyzing across sessions vs collapsed
    % sessions:
    % excludeTrialTriplet: Find session switch trials so we don't process them as if a new
    % session was a continuation of one big session
    if acrossSession
        if ismember(sessionArray, {'concat','behavior1', 'behavior2', 'neural1', 'neural2'})
            nSession = max(trialData.sessionTag);
        else
            nSession = length(sessionArrya);
        end
        excludeTrialTriplet = [];
    else
        nSession = 1;
        excludeTrialPair = find(diff(trialData.sessionTag) < 0);
        excludeTrialTriplet = [excludeTrialPair; excludeTrialPair-1]; % Exclude last 2 trials of a session as possible beginning trials in triplets
    end
    
    
    
    
    % Initialize vectors for per-session RT means
    % Overall session RTs
    rtNs          = nan(nSession, 1);
    
    % Pairs
    rtNsNs1     = nan(nSession, 1);
    rtNsNs2     = nan(nSession, 1);
    rtCNs1      = nan(nSession, 1);
    rtCNs2      = nan(nSession, 1);
    rtNcNs1     = nan(nSession, 1);
    rtNcNs2     = nan(nSession, 1);
    rtENs1      = nan(nSession, 1);
    rtENs2      = nan(nSession, 1);
    
    % Triplets: no-stop correct choices -> various outcomes -> no-stop correct choices
    nNsCNs      = nan(nSession, 1);
    rtNsCNs1    = nan(nSession, 1);
    rtNsCNs3    = nan(nSession, 1);
    nNsNcNs     = nan(nSession, 1);
    rtNsNcNs1   = nan(nSession, 1);
    rtNsNcNs2   = nan(nSession, 1);
    rtNsNcNs3   = nan(nSession, 1);
    nNsENs      = nan(nSession, 1);
    rtNsENs1    = nan(nSession, 1);
    rtNsENs2    = nan(nSession, 1);
    rtNsENs3    = nan(nSession, 1);
    nNsNsNs     = nan(nSession, 1);
    rtNsNsNs1   = nan(nSession, 1);
    rtNsNsNs2   = nan(nSession, 1);
    rtNsNsNs3   = nan(nSession, 1);
    
    % Triplets that end in no-stop error: no-stop correct choices -> various outcomes -> no-stop error choices
    nNsCNse      = nan(nSession, 1);
    rtNsCNse1    = nan(nSession, 1);
    rtNsCNse3    = nan(nSession, 1);
    nNsNcNse     = nan(nSession, 1);
    rtNsNcNse1   = nan(nSession, 1);
    rtNsNcNse2   = nan(nSession, 1);
    rtNsNcNse3   = nan(nSession, 1);
    nNsENse      = nan(nSession, 1);
    rtNsENse1    = nan(nSession, 1);
    rtNsENse2    = nan(nSession, 1);
    rtNsENse3    = nan(nSession, 1);
    nNsNsNse     = nan(nSession, 1);
    rtNsNsNse1   = nan(nSession, 1);
    rtNsNsNse2   = nan(nSession, 1);
    rtNsNsNse3   = nan(nSession, 1);
    
    for j = 1 : nSession
        
        if acrossSession && ismember(sessionArray, {'concat','behavior1', 'behavior2', 'neural1', 'neural2'})
            jTD = trialData(trialData.sessionTag == j, :);
        else
            jTD = trialData;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %      PARSE DATA AND CALCULATE METRICS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iSubject = subjectArray{i};
        % Total mean No-stop RT
        opt = ccm_options;
        opt.outcome = {'goCorrectTarget', 'goCorrectDistractor'};
        opt.ssd = 'none';
        nsTrial = ccm_trial_selection(jTD, opt);
        rtNs(j) = nanmean(jTD.rt(nsTrial));
        
        
        %--------------------------------------------------------------------------
        %       PAIRS
        %--------------------------------------------------------------------------
        Opt2(1) = opt; % Initialize structure with 2 levels
        Opt2(2) = opt;
        
        % 2nd trial is always goCorrectTarget
        Opt2(2).outcome = {'goCorrectTarget'};
        Opt2(2).ssd = 'none';
        
        %--------------------------------------------------------------------------
        % NS -> NS
        %         disp('NoStop - NoStop')
        
        Opt2(1).outcome = {'goCorrectTarget'};
        Opt2(1).ssd = 'none';
        
        rtNsNsTrial = ccm_trial_sequence(jTD, Opt2);
        rtNsNsTrial = setxor(rtNsNsTrial, excludeTrialTriplet);
        rtNsNs1(j) = nanmean(jTD.rt(rtNsNsTrial));
        rtNsNs2(j) = nanmean(jTD.rt(rtNsNsTrial + 1));
        [h,p,ci,stats] = ttest2(jTD.rt(rtNsNsTrial), jTD.rt(rtNsNsTrial+1));
        
        %--------------------------------------------------------------------------
        % C -> NS
        %         disp('Canceled - NoStop')
        
        Opt2(1).outcome = {'stopCorrect'};
        Opt2(1).ssd = 'any';
              
        CNsTrial = ccm_trial_sequence(jTD, Opt2);
        CNsTrial = setxor(CNsTrial, excludeTrialTriplet);
        rtCNs1(j) = nanmean(jTD.rt(CNsTrial));
        rtCNs2(j) = nanmean(jTD.rt(CNsTrial + 1));
        [h,p,ci,stats] = ttest2(jTD.rt(CNsTrial), jTD.rt(CNsTrial+1));
        
        %--------------------------------------------------------------------------
        % NC -> NS
        %         disp('Noncanceled - NoStop')
        
        Opt2(1).outcome = {'stopIncorrectTarget','stopIncorrectDistractor','targetHoldAbort','distractorHoldAbort'};
        Opt2(1).ssd = 'any';
        
        rtNcNsTrial = ccm_trial_sequence(jTD, Opt2);
        rtNcNsTrial = setxor(rtNcNsTrial, excludeTrialTriplet);
        rtNcNs1(j) = nanmean(jTD.rt(rtNcNsTrial));
        rtNcNs2(j) = nanmean(jTD.rt(rtNcNsTrial + 1));
        [h,p,ci,stats] = ttest2(jTD.rt(rtNcNsTrial), jTD.rt(rtNcNsTrial+1));
        
        %--------------------------------------------------------------------------
        % E -> NS    No-stop Error Choice -> No-stop Correct Choice
        %         disp('Error - NoStop')
        
        Opt2(1).outcome = {'goCorrectDistractor'};
        Opt2(1).ssd = 'none';
        
        rtENsTrial = ccm_trial_sequence(jTD, Opt2);
        rtENsTrial = setxor(rtENsTrial, excludeTrialTriplet);
        rtENs1(j) = nanmean(jTD.rt(rtENsTrial));
        rtENs2(j) = nanmean(jTD.rt(rtENsTrial + 1));
        [h,p,ci,stats] = ttest2(jTD.rt(rtENsTrial), jTD.rt(rtENsTrial+1));
        
        
        
        
        
        
        
        
        
        
        
        
        %--------------------------------------------------------------------------
        %       TRIPLETS
        %--------------------------------------------------------------------------
        Opt3(1) = opt; % Initialize structure with 3 levels
        Opt3(2) = opt;
        Opt3(3) = opt;
        
        
        %--------------------------------------------------------------------------
        % NS -> C -> NS
        %         disp('NoStop - Canceled - NoStop')
        
        Opt3(1).outcome     = {'goCorrectTarget'};
        Opt3(1).ssd         = 'none';
        
        Opt3(2).outcome     = {'stopCorrect'};
        Opt3(2).ssd         = 'any';
        
        Opt3(3).outcome     = {'goCorrectTarget'};
        Opt3(3).ssd         = 'none';
        
        rtNsCNsTrial        = ccm_trial_sequence(jTD, Opt3);
        rtNsCNsTrial        = setxor(rtNsCNsTrial, excludeTrialTriplet);
        nNsCNs(j)           = length(rtNsCNsTrial);
        rtNsCNs1(j)         = nanmean(jTD.rt(rtNsCNsTrial));
        rtNsCNs3(j)         = nanmean(jTD.rt(rtNsCNsTrial + 2));
        [h,p,ci,stats]      = ttest2(jTD.rt(rtNsCNsTrial), jTD.rt(rtNsCNsTrial+2));
        
        
        
        %--------------------------------------------------------------------------
        % NS -> NC -> NS
        %          disp('NoStop - NonCanceled - NoStop')
        
        Opt3(1).outcome     = {'goCorrectTarget'};
        Opt3(1).ssd         = 'none';
        
        %         Opt3(2).outcome = {'stopIncorrectTarget', 'stopIncorrectDistractor', 'targetHoldAbort', 'distractorHoldAbort'};
        Opt3(2).outcome     = {'stopIncorrectTarget', 'targetHoldAbort'};
        Opt3(2).ssd         = 'any';
        
        % Opt3(3).outcome = {'goCorrectTarget', 'goCorrectDistractor'};
        Opt3(3).outcome     = {'goCorrectTarget'};
        Opt3(3).ssd         = 'none';
        
        rtNsNcNsTrial       = ccm_trial_sequence(jTD, Opt3);
        rtNsNcNsTrial       = setxor(rtNsNcNsTrial, excludeTrialTriplet);
        nNsNcNs(j)        	= length(rtNsNcNsTrial);
        rtNsNcNs1(j)        = nanmean(jTD.rt(rtNsNcNsTrial));
        rtNsNcNs2(j)        = nanmean(jTD.rt(rtNsNcNsTrial + 1));
        rtNsNcNs3(j)        = nanmean(jTD.rt(rtNsNcNsTrial + 2));
        [h,p,ci,stats]      = ttest2(jTD.rt(rtNsNcNsTrial), jTD.rt(rtNsNcNsTrial+2));
        
        
        
        %--------------------------------------------------------------------------
        % NS -> Error -> NS
        %         disp('NoStop - Choice Error - NoStop')
        
        Opt3(1).outcome     = {'goCorrectTarget'};
        Opt3(1).ssd         = 'none';
        
        Opt3(2).outcome  	= {'goCorrectDistractor'};
        Opt3(2).ssd         = 'none';
        
        Opt3(3).outcome     = {'goCorrectTarget'};
        Opt3(3).ssd         = 'none';
        
        rtNsENsTrial        = ccm_trial_sequence(jTD, Opt3);
        rtNsENsTrial        = setxor(rtNsENsTrial, excludeTrialTriplet);
        nNsENs(j)        	= length(rtNsENsTrial);
        rtNsENs1(j)         = nanmean(jTD.rt(rtNsENsTrial));
        rtNsENs2(j)         = nanmean(jTD.rt(rtNsENsTrial + 1));
        rtNsENs3(j)         = nanmean(jTD.rt(rtNsENsTrial + 2));
        [h,p,ci,stats]      = ttest2(jTD.rt(rtNsENsTrial), jTD.rt(rtNsENsTrial+2));
        
        
        
        %--------------------------------------------------------------------------
        % NS -> NS -> NS
        %          disp('NoStop - Choice Error - NoStop')
        
        Opt3(1).outcome     = {'goCorrectTarget'};
        Opt3(1).ssd         = 'none';
        
        Opt3(2).outcome     = {'goCorrectTarget'};
        Opt3(2).ssd         = 'none';
        
        Opt3(3).outcome     = {'goCorrectTarget'};
        Opt3(3).ssd         = 'none';
        
        rtNsNsNsTrial       = ccm_trial_sequence(jTD, Opt3);
        rtNsNsNsTrial       = setxor(rtNsNsNsTrial, excludeTrialTriplet);
        nNsNsNs(j)        	= length(rtNsNsNsTrial);
        rtNsNsNs1(j)        = nanmean(jTD.rt(rtNsNsNsTrial));
        rtNsNsNs2(j)        = nanmean(jTD.rt(rtNsNsNsTrial + 1));
        rtNsNsNs3(j)        = nanmean(jTD.rt(rtNsNsNsTrial + 2));
        [h,p,ci,stats]      = ttest2(jTD.rt(rtNsNsNsTrial), jTD.rt(rtNsNsNsTrial+2));
        
        
        
        
        
        
        % ERROR ANALYSIS
        
        
        
        %--------------------------------------------------------------------------
        % NS -> C -> NSe (no-stop choice errors
        %         disp('NoStop - Canceled - NoStopError')
        
        Opt3(1).outcome     = {'goCorrectTarget'};
        Opt3(1).ssd         = 'none';
        
        Opt3(2).outcome     = {'stopCorrect'};
        Opt3(2).ssd         = 'any';
        
        Opt3(3).outcome     = {'goCorrectDistractor'};
        Opt3(3).ssd         = 'none';
        
        rtNsCNseTrial       = ccm_trial_sequence(jTD, Opt3);
        rtNsCNseTrial       = setxor(rtNsCNseTrial, excludeTrialTriplet);
        nNsCNse(j)          = length(rtNsCNseTrial);
        rtNsCNse1(j)        = nanmean(jTD.rt(rtNsCNseTrial));
        rtNsCNse3(j)        = nanmean(jTD.rt(rtNsCNseTrial + 2));
        [h,p,ci,stats]      = ttest2(jTD.rt(rtNsCNseTrial), jTD.rt(rtNsCNseTrial+2));
        
        
        %--------------------------------------------------------------------------
        % NS -> NC -> NSe
        %          disp('NoStop - NonCanceled - NoStopError')
        
        Opt3(1).outcome     = {'goCorrectTarget'};
        Opt3(1).ssd         = 'none';
        
        %         Opt3(2).outcome = {'stopIncorrectTarget', 'stopIncorrectDistractor', 'targetHoldAbort', 'distractorHoldAbort'};
        Opt3(2).outcome     = {'stopIncorrectTarget', 'targetHoldAbort'};
        Opt3(2).ssd         = 'any';
        
        % Opt3(3).outcome = {'goCorrectTarget', 'goCorrectDistractor'};
        Opt3(3).outcome     = {'goCorrectDistractor'};
        Opt3(3).ssd         = 'none';
        
        rtNsNcNseTrial       = ccm_trial_sequence(jTD, Opt3);
        rtNsNcNseTrial       = setxor(rtNsNcNseTrial, excludeTrialTriplet);
        nNsNcNse(j)        	= length(rtNsNcNseTrial);
        rtNsNcNse1(j)        = nanmean(jTD.rt(rtNsNcNseTrial));
        rtNsNcNse2(j)        = nanmean(jTD.rt(rtNsNcNseTrial + 1));
        rtNsNcNse3(j)        = nanmean(jTD.rt(rtNsNcNseTrial + 2));
        [h,p,ci,stats]      = ttest2(jTD.rt(rtNsNcNseTrial), jTD.rt(rtNsNcNseTrial+2));
        
        
        
        %--------------------------------------------------------------------------
        % NS -> Error -> NSe
        %         disp('NoStop - Choice Error - NoStopError')
        
        Opt3(1).outcome     = {'goCorrectTarget'};
        Opt3(1).ssd         = 'none';
        
        Opt3(2).outcome  	= {'goCorrectDistractor'};
        Opt3(2).ssd         = 'none';
        
        Opt3(3).outcome     = {'goCorrectDistractor'};
        Opt3(3).ssd         = 'none';
        
        rtNsENseTrial        = ccm_trial_sequence(jTD, Opt3);
        rtNsENseTrial        = setxor(rtNsENseTrial, excludeTrialTriplet);
        nNsENse(j)        	= length(rtNsENseTrial);
        rtNsENse1(j)         = nanmean(jTD.rt(rtNsENseTrial));
        rtNsENse2(j)         = nanmean(jTD.rt(rtNsENseTrial + 1));
        rtNsENse3(j)         = nanmean(jTD.rt(rtNsENseTrial + 2));
        [h,p,ci,stats]      = ttest2(jTD.rt(rtNsENseTrial), jTD.rt(rtNsENseTrial+2));
        
        
        
        %--------------------------------------------------------------------------
        % NS -> NS -> NSe
        %          disp('NoStop - Choice Error - NoStopError')
        
        Opt3(1).outcome     = {'goCorrectTarget'};
        Opt3(1).ssd         = 'none';
        
        Opt3(2).outcome     = {'goCorrectTarget'};
        Opt3(2).ssd         = 'none';
        
        Opt3(3).outcome     = {'goCorrectDistractor'};
        Opt3(3).ssd         = 'none';
        
        rtNsNsNseTrial       = ccm_trial_sequence(jTD, Opt3);
        rtNsNsNseTrial       = setxor(rtNsNsNseTrial, excludeTrialTriplet);
        nNsNsNse(j)        	= length(rtNsNsNseTrial);
        rtNsNsNs1(j)        = nanmean(jTD.rt(rtNsNsNseTrial));
        rtNsNsNs2(j)        = nanmean(jTD.rt(rtNsNsNseTrial + 1));
        rtNsNsNs3(j)        = nanmean(jTD.rt(rtNsNsNseTrial + 2));
        [h,p,ci,stats]      = ttest2(jTD.rt(rtNsNsNseTrial), jTD.rt(rtNsNsNseTrial+2));
        
        
        
        
        %   	plot([1 21], [rtNs(j) rtNs(j)], '--', 'color', colorArray{i})
        %     plot([1:2], [rtNsNs1(j) rtNsNs2(j)], '--o', 'color', colorArray{i})
        %     plot([4], [rtCNs2(j)], '--o', 'color', colorArray{i})
        %     plot([5:6], [rtNcNs1(j) rtNcNs2(j)], '--o', 'color', colorArray{i})
        %     plot([7:8], [rtENs1(j) rtENs2(j)], '--o', 'color', colorArray{i})
        %     plot([10 11 12], [rtNsNsNs1(j) rtNsNsNs2(j) rtNsNsNs3(j)], '-o', 'color', colorArray{i})
        %     plot([13 15], [rtNsCNs1(j) rtNsCNs3(j)], '-o', 'color', colorArray{i})
        %     plot([16 17 18], [rtNsNcNs1(j) rtNsNcNs2(j) rtNsNcNs3(j)], '-o', 'color', colorArray{i})
        %     plot([19 20 21], [rtNsENs1(j) rtNsENs2(j) rtNsENs3(j)], '-o', 'color', colorArray{i})
        %
        
        
    end % for j = 1 : nSession
    

    
    
    if Opt.plotFlag
   
    pNsCNse  	= nNsCNs ./ (nNsCNse + nNsCNs);
    pNsNcNse    = nNsNcNs ./ (nNsNcNse + nNsNcNs);
    pNsENse     = nNsENs ./ (nNsENse + nNsENs);
    pNsNsNse    = nNsNsNs ./ (nNsNsNse + nNsNsNs);
    %     pNsCNse  	= sum(nNsCNs) / sum([nNsCNse;nNsCNs]);
    %     pNsNcNse    = sum(nNsNcNs) / sum([nNsNcNse;nNsNcNs]);
    %     pNsENse     = sum(nNsENs) / sum([nNsENse;nNsENs]);
    %     pNsNsNse    = sum(nNsNsNs) / sum([nNsNsNse;nNsNsNs]);
    fprintf('%s accuracy after trial type:\n', subjectArray{i})
    fprintf('No-stop Correct:\t%0.3f\n', pNsNsNse)
    fprintf('Canceled:\t\t%0.3f\n', pNsCNse)
    fprintf('Nonanceled:\t\t%0.3f\n', pNsNcNse)
    fprintf('No-stop Error:\t\t%0.3f\n', pNsENse)
    probFig = figure(784);
    clf
    hold on;
    yData = [mean(pNsNsNse), mean(pNsCNse), mean(pNsNcNse), mean(pNsENse)];
    bar(yData)
    errorbar(1:4, yData, [sem(pNsNsNse) sem(pNsCNse) sem(pNsNcNse) sem(pNsENse)], 'linestyle' , 'none', 'color', 'k', 'linewidth' , 2)
    ylim([.7 .9])
    print(probFig, [local_figure_path, sprintf('%s_%s_prob_correct_after_trial', subjectArray{1}, sessionArray)],'-dpdf', '-r300')
    
    
    
    
    
    
    % ********************************************************
    % Adjust each session's RTs to the mean overall no-stop Rt and plot
    % ********************************************************
    switch iSubject
        case 'human'
            set(axPair, 'yLim',[550 750])
            set(axTrip, 'yLim',[550 750])
        case 'broca'
            yLimit = [150 400];
            yLimit = [200 300];
            set(axPair, 'yLim',yLimit)
            set(axTrip, 'yLim',yLimit)
        case 'xena'
            set(axPair, 'yLim',[200 450])
            set(axTrip, 'yLim',[200 450])
    end
    set(axPair, 'xLim', [.5 8.5])
    set(axPair, 'xtick', [1.5 3.5 5.5 7.5])
    set(axPair, 'xticklabel', {'NS-NS','C-NS','NC-NS','E-NS'})
    
    set(axTrip, 'xLim', [.5 12.5])
    set(axTrip, 'xtick', [2 5 8 11])
    set(axTrip, 'xticklabel', {'NS-NS-NS','NS-C-NS','NS-NC-NS','NS-E-NS'})
    
    
    
    jMarkererSize = 6;
    jColor = [.5 .5 .5];
    %     jColor = 'k';
    for j = 1 : nSession
        
        if adjustForMean
            meanAdjust = rtNs(j) - nanmean(rtNs);
        else
            meanAdjust = 0;
        end
        
        %   	plot([1 21], [rtNs(j)-meanAdjust rtNs(j)-meanAdjust], '--', 'color', jColor)
        plot(axPair, [1:2], [rtNsNs1(j)-meanAdjust rtNsNs2(j)-meanAdjust], '-o', 'color', jColor, 'markerSize', jMarkererSize, 'markerFaceColor', jColor, 'markerEdgeColor', jColor)
        plot(axPair, [4], [rtCNs2(j)-meanAdjust], '-o', 'color', jColor, 'markerSize', jMarkererSize, 'markerFaceColor', jColor, 'markerEdgeColor', jColor)
        plot(axPair, [5:6], [rtNcNs1(j)-meanAdjust rtNcNs2(j)-meanAdjust], '-o', 'color', jColor, 'markerSize', jMarkererSize, 'markerFaceColor', jColor, 'markerEdgeColor', jColor)
        plot(axPair, [7:8], [rtENs1(j)-meanAdjust rtENs2(j)-meanAdjust], '-o', 'color', jColor, 'markerSize', jMarkererSize, 'markerFaceColor', jColor, 'markerEdgeColor', jColor)
        
        
        plot(axTrip, [1 2 3], [rtNsNsNs1(j)-meanAdjust rtNsNsNs2(j)-meanAdjust rtNsNsNs3(j)-meanAdjust], '-o', 'color', jColor, 'markerSize', jMarkererSize, 'markerFaceColor', jColor, 'markerEdgeColor', jColor)
        plot(axTrip, [4 6], [rtNsCNs1(j)-meanAdjust rtNsCNs3(j)-meanAdjust], '-o', 'color', jColor, 'markerSize', jMarkererSize, 'markerFaceColor', jColor, 'markerEdgeColor', jColor)
        plot(axTrip, [7 8 9], [rtNsNcNs1(j)-meanAdjust rtNsNcNs2(j)-meanAdjust rtNsNcNs3(j)-meanAdjust], '-o', 'color', jColor, 'markerSize', jMarkererSize, 'markerFaceColor', jColor, 'markerEdgeColor', jColor)
        plot(axTrip, [10 11 12], [rtNsENs1(j)-meanAdjust rtNsENs2(j)-meanAdjust rtNsENs3(j)-meanAdjust], '-o', 'color', jColor, 'markerSize', jMarkererSize, 'markerFaceColor', jColor, 'markerEdgeColor', jColor)
    end % for j = 1 : nSession
    
    
    
    % ********************************************************
    % Plot mean rts from all the sessions
    % ********************************************************
    iMarkerSize = 15;
    grayFill = [.9 .9 .9];
    
    % Pairs
    axes(axPair)
    switch sessionArray
        case {'concat', 'behavior1', 'behavior2', 'neural1', 'neural2'}
            h = fill([0 9 9 0], [nanmean(rtNs)-sem(rtNs) nanmean(rtNs)-sem(rtNs) nanmean(rtNs)+sem(rtNs) nanmean(rtNs)+sem(rtNs)], grayFill);
        otherwise
            h = fill([0 9 9 0], [nanmean(rtNs)-sem(jTD.rt(nsTrial)) nanmean(rtNs)-sem(jTD.rt(nsTrial)) nanmean(rtNs)+sem(jTD.rt(nsTrial)) nanmean(rtNs)+sem(jTD.rt(nsTrial))], grayFill);
    end
    set(h, 'edgecolor', 'none');
    
    plot([0 9], [nanmean(rtNs) nanmean(rtNs)], '--', 'color', 'k')
    
    plot(axPair, [1:2], [nanmean(rtNsNs1) nanmean(rtNsNs2)], '-o', 'color', nsColor, 'markerSize', iMarkerSize, 'markerFaceColor', nsColor)
    plot(axPair, [4], [nanmean(rtCNs2)], '-o', 'color', cColor, 'markerSize', iMarkerSize, 'markerFaceColor', cColor)
    plot(axPair, [5:6], [nanmean(rtNcNs1) nanmean(rtNcNs2)], '-o', 'color', ncColor, 'markerSize', iMarkerSize, 'markerFaceColor', ncColor)
    plot(axPair, [7:8], [nanmean(rtENs1) nanmean(rtENs2)], '-o', 'color', eColor, 'markerSize', iMarkerSize, 'markerFaceColor', eColor)
    
    % Triplets
    axes(axTrip)
    switch sessionArray
        case {'concat', 'behavior1', 'behavior2', 'neural1', 'neural2'}
            h = fill([0 13 13 0], [nanmean(rtNs)-sem(rtNs) nanmean(rtNs)-sem(rtNs) nanmean(rtNs)+sem(rtNs) nanmean(rtNs)+sem(rtNs)], grayFill);
        otherwise
            h = fill([0 13 13 0], [nanmean(rtNs)-sem(jTD.rt(nsTrial)) nanmean(rtNs)-sem(jTD.rt(nsTrial)) nanmean(rtNs)+sem(jTD.rt(nsTrial)) nanmean(rtNs)+sem(jTD.rt(nsTrial))], grayFill);
    end
    set(h, 'edgecolor', 'none');
    
    plot(axTrip, [0 13], [nanmean(rtNs) nanmean(rtNs)], '--', 'color', 'k')
    plot(axTrip, [1 2 3], [nanmean(rtNsNsNs1) nanmean(rtNsNsNs2) nanmean(rtNsNsNs3)], '-o', 'color', nsColor, 'markerSize', iMarkerSize, 'markerFaceColor', nsColor)
    plot(axTrip, [4 6], [nanmean(rtNsCNs1) nanmean(rtNsCNs3)], '-o', 'color', cColor, 'markerSize', iMarkerSize, 'markerFaceColor', cColor)
    plot(axTrip, [7 8 9], [nanmean(rtNsNcNs1) nanmean(rtNsNcNs2) nanmean(rtNsNcNs3)], '-o', 'color', ncColor, 'markerSize', iMarkerSize, 'markerFaceColor', ncColor)
    plot(axTrip, [10 11 12], [nanmean(rtNsENs1) nanmean(rtNsENs2) nanmean(rtNsENs3)], '-o', 'color', eColor, 'markerSize', iMarkerSize, 'markerFaceColor', eColor)
    
    fprintf('RT paired and triplet t-tests:\n')
    [h,p,ci,stats]      = ttest2(rtNsNs1, rtNsNs2);
    fprintf('NsNs:\t%0.3f\n', p)
    [h,p,ci,stats]      = ttest2(rtNcNs1, rtNcNs2);
    fprintf('NcNs:\t%0.3f\n', p)
    [h,p,ci,stats]      = ttest2(rtENs1, rtENs2);
    fprintf('ENs:\t%0.3f\n', p)
    [h,p,ci,stats]      = ttest2(rtNsNsNs1, rtNsNsNs3);
    fprintf('NsNsNs:\t%0.3f\n', p)
    [h,p,ci,stats]      = ttest2(rtNsCNs1, rtNsCNs3);
    fprintf('NsCNs:\t%0.3f\n', p)
    [h,p,ci,stats]      = ttest2(rtNsNcNs1, rtNsNcNs3);
    fprintf('NsNcNs:\t%0.3f\n', p)
    [h,p,ci,stats]      = ttest2(rtNsENs1, rtNsENs3);
    fprintf('NsENs:\t%0.3f\n', p)
    
    end
    
if Opt.printPlot
    print(figureHandle,[local_figure_path, subjectArray{1}, '_rt_triplets_',sessionArray, '.pdf'],'-dpdf', '-r300')
end
    
    
end





if Opt.plotFlag
switch sessionArray
    case {'concat','behavior1', 'behavior2', 'neural1', 'neural2'}
    otherwise
        [rt, outlierTrial] = truncate_rt(trialData.rt, 100, 1200);
        trialData.rt(outlierTrial) = nan;
        iMarkerSize = 10;
        figureHandle = 26;
        [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(3, 1, figureHandle);
        clf
        ax2 = axes('units', 'centimeters', 'position', [xAxesPosition(1) yAxesPosition(1) axisWidth axisHeight]);
        hold(ax2, 'all')
        ax3 = axes('units', 'centimeters', 'position', [xAxesPosition(1) yAxesPosition(2) axisWidth axisHeight]);
        hold(ax3, 'all')
        plot(ax3, trialData.rt, '-', 'color', [0 0 0])
        %         plot(ax2, trialData.rt, '-', 'color', [.3 .3 .3])
        %         plot(ax2, smooth(trialData.rt), '-', 'color', [1 0 0])
        plot(ax2, rtNsNsNsTrial+2, trialData.rt(rtNsNsNsTrial+2), 'o', 'markerSize', iMarkerSize, 'markerFaceColor', nsColor, 'markerEdgeColor', nsColor)
        plot(ax2, rtNsCNsTrial+2, trialData.rt(rtNsCNsTrial+2), 'o', 'markerSize', iMarkerSize, 'markerFaceColor', cColor, 'markerEdgeColor', cColor)
        plot(ax2, rtNsNcNsTrial+2, trialData.rt(rtNsNcNsTrial+2), 'o', 'markerSize', iMarkerSize, 'markerFaceColor', ncColor, 'markerEdgeColor', ncColor)
        plot(ax2, rtNsENsTrial+2, trialData.rt(rtNsENsTrial+2), 'o', 'markerSize', iMarkerSize, 'markerFaceColor', eColor, 'markerEdgeColor', eColor)
        if Opt.printPlot
            print(figureHandle,[local_figure_path, subjectArray{1}, '_rt_triplets_',sessionArray, '.pdf'],'-dpdf', '-r300')
        end
end
end


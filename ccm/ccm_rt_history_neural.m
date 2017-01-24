function Data = ccm_rt_history_neural(subjectID, sessionID, Opt)

% triplet analysis:
% Replicating Pouget et al
% These analyses can be done a few ways:
%
%   Include/exclude aborted trials: This will affect the number of paired
%   and triplet trials that make it into analyses, since aborts between
%   trials may or may not count as successive trials. (this seems not to
%   matter though: deleteAborts = true vs. false
%
%   Analyze data across sessions, taking the mean across sessions, or
%   analyze with all data collapsed (as if one big session). Also doesn't
%   alter the results much. acrossSession = true vs. false
%
% Keeping aborted trials in
% without respect to choice difficulty. As a first, dont' remove any aborted
% trials. This will greatly reduce the data, but is a more valid test
%
%    Opt.dataType = 'neuron' (default) 'lfp', 'erp';
%
%    Opt.figureHandle   = 1000;
%    Opt.printPlot      = false, true;
%    Opt.plotFlag       = true, false;
%    Opt.howProcess      = how to step through the list of units
%                                 'each' to plot all,
%                                 'step' (default): step through to see one
%                                 plot at a time, pausing between
%                                 'print' (default): step through each plot
%                                 individually, printing each to file
%    Opt.unitArray      = units you want to analyze (default gets filled
%                                   with all available).
%                                 {'spikeUnit17a'}, input a cell of units for a list of individaul units
%     Opt.responseDir  = the angle of target to which a response was made
%                       (on the last trial of the sequence)
%                       {'collapse','left','right','up','down','upLeft','upRight',downLeft','downRight'<list of possible angels, e.g. [45 -45]
%     Opt.pairTriplet  = perform the analysis on RT pairs (default) or triplets
%     Opt.epochName  = perform the analysis on which epoch?
%     Opt.epochWindow  = [-300:300] (default). Perform the analysis on which epoch window relative to epoch alignment?

if nargin < 3;
    Opt = ccm_options;
    Opt.responseDir = {'left', 'right'};
    Opt.pairTriplet = 'pair';
end
Opt.figureHandle = 66;

Data = struct;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA AND SET VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            switch Opt.dataType
                case 'neuron'
                    addVar = 'spikeData';
                case 'lfp'
                    addVar = 'lfpData';
                case 'eeg'
                    addVar = 'eegData';
            end
            variables = [ccm_min_vars, addVar];
            
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID, variables);
if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a chioce countermanding saccade session, try again\n')
    return
end
pSignalArray    = ExtraVar.pSignalArray;
targAngleArray	= ExtraVar.targAngleArray;
ssdArray        = ExtraVar.ssdArray;
nSSD            = length(ssdArray);




% Set defaults
switch Opt.dataType
    case 'neuron'
        dataArray     = SessionData.spikeUnitArray;
    case 'lfp'
        chNum = SessionData.lfpChannel;
        dataArray 	= num2cell(SessionData.lfpChannel);
        dataArray   = cellfun(@(x) sprintf('lfp_%s', num2str(x, '%02d')), dataArray, 'uniformoutput', false);
    case 'erp'
        dataArray     = eeg_electrode_map(subjectID);
end
if isempty(Opt.unitArray)
    Opt.unitArray     = dataArray;
end

% Make sure user input a Opt.dataType that was recorded during the session
dataTypePossible = {'neuron', 'lfp', 'erp'};
if ~sum(strcmp(Opt.dataType, dataTypePossible))
    fprintf('%s Is not a valid data type \n', Opt.dataType)
    return
end
if isempty(Opt.unitArray)
    fprintf('Session %s apparently does not contain %s data \n', sessionID, Opt.dataType)
    return
end





% CONSTANTS
MIN_RT          = 120;
MAX_RT          = 1200;
STD_MULTIPLE    = 3;
DELETE_ABORTS   = true;
N_UNIT          = length(Opt.unitArray);

cropWindow  = -1000 : 1500;  % used to extract a semi-small portion of signal for each epoch/alignemnt
baseWindow 	= -149 : 0;   % To baseline-shift the eeg signals, relative to event alignment index;

Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;

% Kernel.method = 'gaussian';
% Kernel.sigma = 15;

% Get rid of trials with outlying RTs
[allRT, rtOutlierTrial] = truncate_rt(trialData.rt, MIN_RT, MAX_RT, STD_MULTIPLE);
trialData.rt(rtOutlierTrial) = nan;




if DELETE_ABORTS
    %     selectOpt = ccm_trial_selection;
    %     selectOpt.outcome = {...
    %         'goCorrectTarget', 'goCorrectDistractor', ...
    %         'stopCorrect', ...
    %         'targetHoldAbort', 'distractorHoldAbort', ...
    %         'stopIncorrectTarget', 'stopIncorrectDistractor'};
    %     validTrial = ccm_trial_selection(trialData, selectOpt);
    %     trialData = trialData(validTrial,:);
    selectOpt = ccm_trial_selection;
    selectOpt.outcome = {...
        'noFixation', 'fixationAbort'};
    invalidTrial = ccm_trial_selection(trialData, selectOpt);
    trialData(invalidTrial,:) = [];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsColor     = [0 .7 0];
cColor      = 'r';
ncColor     = 'k';
eColor      = 'b';

nRow        = 1;
nCol        = length(Opt.responseDir);
lineWidth   = 3;

for kDataInd = 1 : N_UNIT
    
    switch Opt.dataType
        case 'neuron'
            [a, kUnit] = ismember(Opt.unitArray{kDataInd}, SessionData.spikeUnitArray);
            yLimit      = [0 80];
        case 'lfp'
            [a, kUnit] = ismember(chNum(kDataInd), SessionData.lfpChannel);
            yLimit      = [-.1 .1];
        case 'erp'
            [a, kUnit] = ismember(Opt.unitArray{kDataInd}, eeg_electrode_map(subjectID));
            yLimit      = [-.1 .1];
    end
    
    
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nCol, Opt.figureHandle);
     axisWidth = axisWidth * .95;
    clf
%             title(sprintf('%s_%s',sessionID, Opt.unitArray{kDataInd}))
    
    
    
    
    
    
    
    
    
    for dirInd = 1 : length(Opt.responseDir)
        iResponseDir = Opt.responseDir{dirInd};
        
        ax(1,dirInd) = axes('units', 'centimeters', 'position', [xAxesPosition(1,dirInd) yAxesPosition(1,dirInd) axisWidth axisHeight]);
        hold(ax(1,dirInd), 'on')
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %      PARSE DATA AND CALCULATE METRICS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Total mean No-stop RT
        opt = ccm_trial_selection;
        opt.outcome = {'goCorrectTarget', 'goCorrectDistractor'};
        opt.ssd = 'none';
        nsTrial = ccm_trial_selection(trialData, opt);
        rtNs = nanmean(trialData.rt(nsTrial));
        
        
        
        
        switch Opt.pairTriplet
            case 'pair'
                addTrial = 1;
                %--------------------------------------------------------------------------
                %       PAIRS
                %--------------------------------------------------------------------------
                Opt2(1) = opt; % Initialize structure with 2 levels
                Opt2(2) = opt;
                
                % 2nd trial is always goCorrectTarget
                Opt2(2).outcome = {'goCorrectTarget'};
                Opt2(2).ssd = 'none';
                Opt2(2).responseDir = iResponseDir;
                
                
                %--------------------------------------------------------------------------
                % NS -> NS
                %         disp('NoStop - NoStop')
                
                Opt2(1).outcome     = {'goCorrectTarget'};
                Opt2(1).ssd         = 'none';
                
                rtNsTrial = ccm_trial_sequence(trialData, Opt2);
                
                %--------------------------------------------------------------------------
                % C -> NS
                %         disp('Canceled - NoStop')
                
                Opt2(1).outcome = {'stopCorrect'};
                Opt2(1).ssd = 'any';
                
                rtCTrial = ccm_trial_sequence(trialData, Opt2);
                  
                
                %--------------------------------------------------------------------------
                % NC -> NS
                %         disp('Noncanceled - NoStop')
                
                Opt2(1).outcome = {'stopIncorrectTarget','stopIncorrectDistractor','targetHoldAbort','distractorHoldAbort'};
                Opt2(1).ssd = 'any';
                
                rtNcTrial = ccm_trial_sequence(trialData, Opt2);
                
                
                %--------------------------------------------------------------------------
                % E -> NS    No-stop Error Choice -> No-stop Correct Choice
                %         disp('Error - NoStop')
                
                Opt2(1).outcome = {'goCorrectDistractor'};
                Opt2(1).ssd = 'none';
                
                rtETrial = ccm_trial_sequence(trialData, Opt2);
                  
                
            case 'triplet'
                addTrial = 2;
                %--------------------------------------------------------------------------
                %       TRIPLETS
                %--------------------------------------------------------------------------
                Opt3(1) = opt; % Initialize structure with 3 levels
                Opt3(2) = opt;
                Opt3(3) = opt;
                
                
                %--------------------------------------------------------------------------
                % NS -> NS -> NS
                %          disp('NoStop - Choice Error - NoStop')
                
                Opt3(1).outcome     = {'goCorrectTarget'};
                Opt3(1).ssd         = 'none';
                
                Opt3(2).outcome     = {'goCorrectTarget'};
                Opt3(2).ssd         = 'none';
                
                Opt3(3).outcome     = {'goCorrectTarget'};
                Opt3(3).ssd         = 'none';
                Opt3(3).responseDir = iResponseDir;
                
                rtNsTrial       = ccm_trial_sequence(trialData, Opt3);
                
                
                %--------------------------------------------------------------------------
                % NS -> C -> NS
                %         disp('NoStop - Canceled - NoStop')
                
                Opt3(1).outcome     = {'goCorrectTarget'};
                Opt3(1).ssd         = 'none';
                
                Opt3(2).outcome     = {'stopCorrect'};
                Opt3(2).ssd         = 'any';
                
                Opt3(3).outcome     = {'goCorrectTarget'};
                Opt3(3).ssd         = 'none';
                Opt3(3).responseDir = iResponseDir;
                
                rtCTrial        = ccm_trial_sequence(trialData, Opt3);
                
                
                
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
                Opt3(3).responseDir = iResponseDir;
                
                rtNcTrial       = ccm_trial_sequence(trialData, Opt3);
                
                
                
                %--------------------------------------------------------------------------
                % NS -> Error -> NS
                %         disp('NoStop - Choice Error - NoStop')
                
                Opt3(1).outcome     = {'goCorrectTarget'};
                Opt3(1).ssd         = 'none';
                
                Opt3(2).outcome  	= {'goCorrectDistractor'};
                Opt3(2).ssd         = 'none';
                
                Opt3(3).outcome     = {'goCorrectTarget'};
                Opt3(3).ssd         = 'none';
                Opt3(3).responseDir = iResponseDir;
                
                rtETrial        = ccm_trial_sequence(trialData, Opt3);
                
                
        end
        
        
        
        
        
        %--------------------------------------------------------------------------
        % NS -> (NS) -> NS
        
        alignList = trialData.(Opt.epochName)(rtNsTrial+addTrial);
        
        switch Opt.dataType
            case 'neuron'
                [rasNs, alignNs] = spike_to_raster(trialData.spikeData(rtNsTrial+1, kUnit), alignList);
                signalNs = spike_density_function(rasNs, Kernel);
            case {'lfp', 'erp'}
                [signalNs, alignNs] 	= align_signals(trialData.lfpData(rtNsTrial+1, kUnit), alignList, cropWindow);
                satTrial                = signal_reject_saturate(signalNs, 'alignIndex', alignNs);
                signalNs(satTrial,:)     = [];
        end
        
        
        
        
        %--------------------------------------------------------------------------
        % NS -> (C) -> NS
        
        
        
        
        alignList = trialData.(Opt.epochName)(rtCTrial+addTrial);
        
        switch Opt.dataType
            case 'neuron'
                [rasC, alignC] = spike_to_raster(trialData.spikeData(rtCTrial+1, kUnit), alignList);
                signalC = spike_density_function(rasC, Kernel);
            case {'lfp', 'erp'}
                [signalC, alignC] 	= align_signals(trialData.lfpData(rtCTrial+1, kUnit), alignList, cropWindow);
                satTrial                = signal_reject_saturate(signalC, 'alignIndex', alignC);
                signalC(satTrial,:)     = [];
        end
        
        
        
        %--------------------------------------------------------------------------
        % NS -> (NC) -> NS
        alignList = trialData.(Opt.epochName)(rtNcTrial+addTrial);
        
        switch Opt.dataType
            case 'neuron'
                [rasNc, alignNc] = spike_to_raster(trialData.spikeData(rtNcTrial+1, kUnit), alignList);
                signalNc = spike_density_function(rasNc, Kernel);
            case {'lfp', 'erp'}
                [signalNc, alignNc] 	= align_signals(trialData.lfpData(rtNcTrial+1, kUnit), alignList, cropWindow);
                satTrial                = signal_reject_saturate(signalNc, 'alignIndex', alignNc);
                signalNc(satTrial,:)     = [];
        end
        
        
        %--------------------------------------------------------------------------
        % NS -> (Error) -> NS
        alignList = trialData.(Opt.epochName)(rtETrial+addTrial);
        
        switch Opt.dataType
            case 'neuron'
                [rasE, alignE] = spike_to_raster(trialData.spikeData(rtETrial+1, kUnit), alignList);
                signalE = spike_density_function(rasE, Kernel);
            case {'lfp', 'erp'}
                [signalE, alignE] 	= align_signals(trialData.lfpData(rtETrial+1, kUnit), alignList, cropWindow);
                satTrial                = signal_reject_saturate(signalE, 'alignIndex', alignE);
                signalE(satTrial,:)     = [];
        end
        
        
        
        
        
        
        % print the figure if we're stepping through
        if Opt.plotFlag
            switch Opt.dataType
                case 'neuron'
                    sigMax = max([nanmean(signalNs(:, Opt.epochWindow + alignNs)), ...
                        nanmean(signalC(:, Opt.epochWindow + alignNs)), ...
                        nanmean(signalE(:, Opt.epochWindow + alignNs)), ...
                        nanmean(signalNc(:, Opt.epochWindow + alignNs))]);
                    yLimit = [0 max(40, ceil(sigMax * 1.2) + 1)];
                case {'erp','lfp'}
                    yLimit      = [-.1 .1];
            end
            set(ax(1,dirInd), 'ylim', yLimit)
            plot(ax(dirInd), nanmean(signalNc(:, Opt.epochWindow + alignNc)), 'color', ncColor, 'linewidth', 1)
            plot(ax(dirInd), nanmean(signalE(:, Opt.epochWindow + alignE)), 'color', eColor, 'linewidth', 1)
            plot(ax(dirInd), nanmean(signalC(:, Opt.epochWindow + alignC)), 'color', cColor, 'linewidth', lineWidth)
            plot(ax(dirInd), nanmean(signalNs(:, Opt.epochWindow + alignNs)), 'color', nsColor, 'linewidth', lineWidth)
            legend('NsCNs','NsNcNs','NsENs','NsNs')
            plot(ax(dirInd), [-Opt.epochWindow(1)+1 -Opt.epochWindow(1)+1], [0 yLimit(2)*.8], 'k')
            title(sprintf('%s',Opt.responseDir{dirInd}))
        end
        
    end
    
    
    
    if Opt.printPlot
        localFigurePath = local_figure_path;
        print(Opt.figureHandle,[localFigurePath, sprintf('%s_rt_history_neural_%s_%s_%s_%s', sessionID, Opt.unitArray{kDataInd}, Opt.pairTriplet, Opt.dataType, Opt.epochName)],'-dpdf', '-r300')
    end
    
    
    
    howProcess = Opt.howProcess;
    switch howProcess
        case {'step','print'}
            % end
            
            if strcmp(howProcess, 'step')
                pause
            end
            clear Data
        otherwise
            Opt.figureHandle = Opt.figureHandle + 1;
            
    end
end








%
%
%
%
%
%
%
%
%
%
%
% % ERROR ANALYSIS
%
%
%
% %--------------------------------------------------------------------------
% % NS -> C -> NSe (no-stop choice errors
% %         disp('NoStop - Canceled - NoStopError')
%
% Opt3(1).outcome     = {'goCorrectTarget'};
% Opt3(1).ssd         = 'none';
%
% Opt3(2).outcome     = {'stopCorrect'};
% Opt3(2).ssd         = 'any';
%
% Opt3(3).outcome     = {'goCorrectDistractor'};
% Opt3(3).ssd         = 'none';
%
% rtNsCNseTrial       = ccm_trial_sequence(trialData, Opt3);
% rtNsCNseTrial       = setxor(rtNsCNseTrial, excludeTrialTriplet);
% nNsCNse          = length(rtNsCNseTrial);
% rtNsCNse1        = nanmean(trialData.rt(rtNsCNseTrial));
% rtNsCNse3        = nanmean(trialData.rt(rtNsCNseTrial + 2));
% [h,p,ci,stats]      = ttest2(trialData.rt(rtNsCNseTrial), trialData.rt(rtNsCNseTrial+2));
%
%
% %--------------------------------------------------------------------------
% % NS -> NC -> NSe
% %          disp('NoStop - NonCanceled - NoStopError')
%
% Opt3(1).outcome     = {'goCorrectTarget'};
% Opt3(1).ssd         = 'none';
%
% %         Opt3(2).outcome = {'stopIncorrectTarget', 'stopIncorrectDistractor', 'targetHoldAbort', 'distractorHoldAbort'};
% Opt3(2).outcome     = {'stopIncorrectTarget', 'targetHoldAbort'};
% Opt3(2).ssd         = 'any';
%
% % Opt3(3).outcome = {'goCorrectTarget', 'goCorrectDistractor'};
% Opt3(3).outcome     = {'goCorrectDistractor'};
% Opt3(3).ssd         = 'none';
%
% rtNsNcNseTrial       = ccm_trial_sequence(trialData, Opt3);
% rtNsNcNseTrial       = setxor(rtNsNcNseTrial, excludeTrialTriplet);
% nNsNcNse        	= length(rtNsNcNseTrial);
% rtNsNcNse1        = nanmean(trialData.rt(rtNsNcNseTrial));
% rtNsNcNse2        = nanmean(trialData.rt(rtNsNcNseTrial + 1));
% rtNsNcNse3        = nanmean(trialData.rt(rtNsNcNseTrial + 2));
% [h,p,ci,stats]      = ttest2(trialData.rt(rtNsNcNseTrial), trialData.rt(rtNsNcNseTrial+2));
%
%
%
% %--------------------------------------------------------------------------
% % NS -> Error -> NSe
% %         disp('NoStop - Choice Error - NoStopError')
%
% Opt3(1).outcome     = {'goCorrectTarget'};
% Opt3(1).ssd         = 'none';
%
% Opt3(2).outcome  	= {'goCorrectDistractor'};
% Opt3(2).ssd         = 'none';
%
% Opt3(3).outcome     = {'goCorrectDistractor'};
% Opt3(3).ssd         = 'none';
%
% rtNsENseTrial        = ccm_trial_sequence(trialData, Opt3);
% rtNsENseTrial        = setxor(rtNsENseTrial, excludeTrialTriplet);
% nNsENse        	= length(rtNsENseTrial);
% rtNsENse1         = nanmean(trialData.rt(rtNsENseTrial));
% rtNsENse2         = nanmean(trialData.rt(rtNsENseTrial + 1));
% rtNsENse3         = nanmean(trialData.rt(rtNsENseTrial + 2));
% [h,p,ci,stats]      = ttest2(trialData.rt(rtNsENseTrial), trialData.rt(rtNsENseTrial+2));
%
%
%
% %--------------------------------------------------------------------------
% % NS -> NS -> NSe
% %          disp('NoStop - Choice Error - NoStopError')
%
% Opt3(1).outcome     = {'goCorrectTarget'};
% Opt3(1).ssd         = 'none';
%
% Opt3(2).outcome     = {'goCorrectTarget'};
% Opt3(2).ssd         = 'none';
%
% Opt3(3).outcome     = {'goCorrectDistractor'};
% Opt3(3).ssd         = 'none';
%
% rtNsNseTrial       = ccm_trial_sequence(trialData, Opt3);
% rtNsNseTrial       = setxor(rtNsNseTrial, excludeTrialTriplet);
% nNsNsNse        	= length(rtNsNseTrial);
% rtNsNs1        = nanmean(trialData.rt(rtNsNseTrial));
% rtNsNs2        = nanmean(trialData.rt(rtNsNseTrial + 1));
% rtNsNs3        = nanmean(trialData.rt(rtNsNseTrial + 2));
% [h,p,ci,stats]      = ttest2(trialData.rt(rtNsNseTrial), trialData.rt(rtNsNseTrial+2));
%
%
%
%
%
% pNsCNse  	= sum(nNsCNse) / sum([nNsCNse;nNsCNsL]);
% pNsNcNse    = sum(nNsNcNse) / sum([nNsNcNse;nNsNcNs]);
% pNsENse     = sum(nNsENse) / sum([nNsENse;nNsENs]);
% pNsNsNse    = sum(nNsNsNse) / sum([nNsNsNse;nNsNsNs]);
% fprintf('%s error probability after trial type:\n', subjectArray{i})
% fprintf('Canceled:\t\t%0.3f\n', pNsCNse)
% fprintf('Nonanceled:\t\t%0.3f\n', pNsNcNse)
% fprintf('No-stop Error:\t\t%0.3f\n', pNsENse)
% fprintf('No-stop Correct:\t%0.3f\n', pNsNsNse)
%
%
%
%
%
% % Plot
% % ylim([250 350])
% plot([1:2], [nanmean(rtNs1) nanmean(rtNs2)], '--o', 'color', colorArray{i})
% plot([4], [nanmean(rtC2)], '--o', 'color', colorArray{i})
% plot([5:6], [nanmean(rtNc1) nanmean(rtNc2)], '--o', 'color', colorArray{i})
% plot([7:8], [nanmean(rtE1) nanmean(rtE2)], '--o', 'color', colorArray{i})


% colorArray = {'k'};
% i = 1;
%
% [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(1, 1, Opt.figureHandle + 100);
% clf
% axes('units', 'centimeters', 'position', [xAxesPosition(1,1) yAxesPosition(1,1) axisWidth axisHeight]);
%          hold on
%
% plot([1 22], [nanmean(rtNs) nanmean(rtNs)], '--', 'color', 'k')
% plot([10 11 12], [nanmean(rtNsNs1) nanmean(rtNsNs2) nanmean(rtNsNs3)], '-o', 'color', nsColor)
% plot([13 15], [nanmean(rtNsCNs1) nanmean(rtNsCNs3)], '-o', 'color', cColor)
% plot([16 17 18], [nanmean(rtNsNcNs1) nanmean(rtNsNcNs2) nanmean(rtNsNcNs3)], '-o', 'color', ncColor)
% plot([19 20 21], [nanmean(rtNsENs1) nanmean(rtNsENs2) nanmean(rtNsENs3)], '-o', 'color', eColor)
%
% %
% %
% %
% %
% %
% %
% if strcmp(iSubject, 'human')
%     ylim([550 750])
% else
%     ylim([0 80])
% end
% xlim([0 22])
% set(gca, 'xtick', [1.5 3.5 5.5 7.5 11 14 17 20])
% set(gca, 'xticklabel', {'NS-NS','C-NS','NC-NS','E-NS','NS-NS-NS','NS-C-NS','NS-NC-NS','NS-E-NS'})
% % % legend({'Broca','Xena'})
% %
% %
% savePlot = 'y';
% %     input('save?', 's');
% if strcmp(savePlot, 'y')

end


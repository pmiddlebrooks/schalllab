%% triplet analysis:
% Keeping aborted trials in
% without respect to choice difficulty. As a first, dont' remove any aborted
% trials. This will greatly reduce the data, but is a more valid test

colorArray = {'b','r'};
figure(1);
clf
hold all;

% monkeyArray = {'broca','xena'};
monkeyArray = {'broca'};
deleteAborts = false;

jTrode = 3;
alignEpoch = 'targOn';
window = -99:250;
for i = 1 : length(monkeyArray)
    iMonkey = monkeyArray{i};
    
    iFile = fullfile('local_data', iMonkey, strcat(iMonkey,'EEG.mat'));
    
    % load('local_data/broca/brocaRT.mat')
    % load('local_data/xena/xenaRT.mat')
    % monkeyB = monkeyX;
    
    load(iFile)
    
    if deleteAborts
        selectOpt = ccm_trial_selection;
        selectOpt.outcome = {...
            'goCorrectTarget', 'goCorrectDistractor', ...
            'stopCorrect', ...
            'targetHoldAbort', 'distractorHoldAbort', ...
            'stopIncorrectTarget', 'stopIncorrectDistractor'};
        validTrial = ccm_trial_selection(trialData, selectOpt);
        trialData = trialData(validTrial,:);
    end
    
    
    
    
    % Total mean No-stop RT
    sOpt = ccm_trial_selection;
    sOpt.outcome = {'goCorrectTarget', 'goCorrectDistractor'};
    sOpt.ssd = 'none';
    nsTrial = ccm_trial_selection(trialData, sOpt);
    ns = nanmean(trialData.rt(nsTrial));
    signalMatrix = trialData.eegData(nsTrial, jTrode);
    alignmentTimeList = trialData.(alignEpoch)(nsTrial);
    [alignNs, alignmentIndex] = align_signals(signalMatrix, alignmentTimeList, window);
    alignNs = signal_baseline_correct(alignNs, 1:-window(1), 1);
    
    
    % NS -> NS
    sOpt(1) = ccm_trial_selection;
    % sOpt(1).outcome = {'goCorrectTarget', 'goCorrectDistractor'};
    sOpt(1).outcome = {'goCorrectTarget'};
    sOpt(1).ssd = 'none';
    
    sOpt(2) = ccm_trial_selection;
    % sOpt(2).outcome = {'goCorrectTarget', 'goCorrectDistractor'};
    sOpt(2).outcome = {'goCorrectTarget'};
    sOpt(2).ssd = 'none';
    
    nsNsTrial = ccm_trial_sequence(trialData, sOpt);
    disp('NoStop - NoStop')
    nsNs1 = nanmean(trialData.rt(nsNsTrial));
    nsNs2 = nanmean(trialData.rt(nsNsTrial + 1));
    [h,p,ci,stats] = ttest2(trialData.rt(nsNsTrial), trialData.rt(nsNsTrial+1));
    signalMatrix = trialData.eegData(nsNsTrial+1, jTrode);
    alignmentTimeList = trialData.(alignEpoch)(nsNsTrial+1);
    [alignNsNs, alignmentIndex] = align_signals(signalMatrix, alignmentTimeList, window);
    alignNsNs = signal_baseline_correct(alignNsNs, 1:-window(1), 1);
    
    
    
    
    % NS -> C -> NS
    sOpt(1) = ccm_trial_selection;
    % sOpt(2).outcome = {'goCorrectTarget', 'goCorrectDistractor'};
    sOpt(1).outcome = {'goCorrectTarget'};
    sOpt(1).ssd = 'none';
    
    sOpt(2) = ccm_trial_selection;
    sOpt(2).outcome = {'stopCorrect'};
    sOpt(2).ssd = 'any';
    
    sOpt(3) = ccm_trial_selection;
    % sOpt(3).outcome = {'goCorrectTarget', 'goCorrectDistractor'};
    sOpt(3).outcome = {'goCorrectTarget'};
    sOpt(3).ssd = 'none';
    
    disp('NoStop - Canceled - NoStop')
    nsCNsTrial = ccm_trial_sequence(trialData, sOpt);
    nsCNs1 = nanmean(trialData.rt(nsCNsTrial));
    nsCNs2 = nanmean(trialData.rt(nsCNsTrial + 2));
    [h,p,ci,stats] = ttest2(trialData.rt(nsCNsTrial), trialData.rt(nsCNsTrial+2));
    signalMatrix = trialData.eegData(nsCNsTrial+2, jTrode);
    alignmentTimeList = trialData.(alignEpoch)(nsCNsTrial+2);
    [alignNsCNs, alignmentIndex] = align_signals(signalMatrix, alignmentTimeList, window);
    alignNsCNs = signal_baseline_correct(alignNsCNs, 1:-window(1), 1);
    
    
    
    % NS -> NC -> NS
    sOpt(1) = ccm_trial_selection;
    % sOpt(1).outcome = {'goCorrectTarget', 'goCorrectDistractor'};
    sOpt(1).outcome = {'goCorrectTarget'};
    sOpt(1).ssd = 'none';
    
    sOpt(2) = ccm_trial_selection;
    sOpt(2).outcome = {'stopIncorrectTarget', 'stopIncorrectDistractor', 'targetHoldAbort', 'distractorHoldAbort'};
    sOpt(2).outcome = {'stopIncorrectTarget', 'targetHoldAbort'};
    sOpt(2).ssd = 'any';
    
    sOpt(3) = ccm_trial_selection;
    % sOpt(3).outcome = {'goCorrectTarget', 'goCorrectDistractor'};
    sOpt(3).outcome = {'goCorrectTarget'};
    sOpt(3).ssd = 'none';
    
    disp('NoStop - NonCanceled - NoStop')
    nsNcNsTrial = ccm_trial_sequence(trialData, sOpt);
    nsNcNs1 = nanmean(trialData.rt(nsNcNsTrial));
    nsNcNs3 = nanmean(trialData.rt(nsNcNsTrial + 1));
    nsNcNs2 = nanmean(trialData.rt(nsNcNsTrial + 2));
    [h,p,ci,stats] = ttest2(trialData.rt(nsNcNsTrial), trialData.rt(nsNcNsTrial+2));
    signalMatrix = trialData.eegData(nsNcNsTrial+2, jTrode);
    alignmentTimeList = trialData.(alignEpoch)(nsNcNsTrial+2);
    [alignNsNcNs, alignmentIndex] = align_signals(signalMatrix, alignmentTimeList, window);
    alignNsNcNs = signal_baseline_correct(alignNsNcNs, 1:-window(1), 1);
    
    
    
    
    % NS -> Error -> NS
    sOpt(1) = ccm_trial_selection;
    sOpt(1).outcome = {'goCorrectTarget'};
    sOpt(1).ssd = 'none';
    
    sOpt(2) = ccm_trial_selection;
    sOpt(2).outcome = {'goCorrectDistractor'};
    sOpt(2).ssd = 'none';
    
    sOpt(3) = ccm_trial_selection;
    sOpt(3).outcome = {'goCorrectTarget'};
    sOpt(3).ssd = 'none';
    
    disp('NoStop - Choice Error - NoStop')
    nsENsTrial = ccm_trial_sequence(trialData, sOpt);
    nsENs1 = nanmean(trialData.rt(nsENsTrial));
    nsENs3 = nanmean(trialData.rt(nsENsTrial + 1));
    nsENs2 = nanmean(trialData.rt(nsENsTrial + 2));
    [h,p,ci,stats] = ttest2(trialData.rt(nsENsTrial), trialData.rt(nsENsTrial+2));
    signalMatrix = trialData.eegData(nsENsTrial+2, jTrode);
    alignmentTimeList = trialData.(alignEpoch)(nsENsTrial+2);
    [alignNsENs, alignmentIndex] = align_signals(signalMatrix, alignmentTimeList, window);
    alignNsENs = signal_baseline_correct(alignNsENs, 1:-window(1), 1);
    
    
    
    % NS -> NS -> NS
    sOpt(1) = ccm_trial_selection;
    sOpt(1).outcome = {'goCorrectTarget'};
    sOpt(1).ssd = 'none';
    
    sOpt(2) = ccm_trial_selection;
    sOpt(2).outcome = {'goCorrectTarget'};
    sOpt(2).ssd = 'none';
    
    sOpt(3) = ccm_trial_selection;
    sOpt(3).outcome = {'goCorrectTarget'};
    sOpt(3).ssd = 'none';
    
    disp('NoStop - Choice Error - NoStop')
    nsNsNsTrial = ccm_trial_sequence(trialData, sOpt);
    nsNsNs1 = nanmean(trialData.rt(nsNsNsTrial));
    nsNsNs3 = nanmean(trialData.rt(nsNsNsTrial + 1));
    nsNsNs2 = nanmean(trialData.rt(nsNsNsTrial + 2));
    [h,p,ci,stats] = ttest2(trialData.rt(nsNsNsTrial), trialData.rt(nsNsNsTrial+2));
    signalMatrix = trialData.eegData(nsNsNsTrial+2, jTrode);
    alignmentTimeList = trialData.(alignEpoch)(nsNsNsTrial+2);
    [alignNsNsNs, alignmentIndex] = align_signals(signalMatrix, alignmentTimeList, window);
    alignNsNsNs = signal_baseline_correct(alignNsNsNs, 1:-window(1), 1);
    
    
    
    
    % Plot
    % ylim([250 350])
    plot([1 10], [ns ns], '--', 'color', colorArray{i})
    plot([1 2], [nsNs1 nsNs2], '-o', 'color', colorArray{i})
    plot([3 4], [nsCNs1 nsCNs2], '-o', 'color', colorArray{i})
    plot([5 5.5 6], [nsNcNs1 nsNcNs3 nsNcNs2], '-o', 'color', colorArray{i})
    plot([7 7.5 8], [nsENs1 nsENs3 nsENs2], '-o', 'color', colorArray{i})
    plot([9 9.5 10], [nsNsNs1 nsNsNs3 nsNsNs2], '-o', 'color', colorArray{i})
end
ylim([200 350])
set(gca, 'xtick', [1 2 3 3.5 4 5 5.5 6 7 7.5 8 9 9.5 10])
set(gca, 'xticklabel', {'NS','NS','NS','C','NS','NS','NC','NS','NS','E','NS','NS','NS','NS'})
% legend({'Broca','Xena'})
%%
figure(9)
clf
hold on;
plot(nanmean(alignNs, 1), 'k--')
plot(nanmean(alignNsNs, 1), 'k')
plot(nanmean(alignNsCNs, 1), 'r')
plot(nanmean(alignNsNcNs, 1), 'b')
plot(nanmean(alignNsENs, 1), 'g')
plot(nanmean(alignNsNsNs, 1), 'k-.')


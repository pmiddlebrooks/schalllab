function data = ccm_n2pc(subjectID, sessionID, varargin)

%%

% subjectID = 'Broca';
% sessionID = 'bp063n01';
% sessionID = 'bp081n01';
% subjectID = 'Xena';
% sessionID = 'xp056n01';
% Set defaults
alignEvent = 'checkerOn';
figureHandle = 654;
plotFlag = 1;
printPlot = 0;
filterData = false;
for i = 1 : 2 : length(varargin)
    switch varargin{i}
        case 'plotFlag'
            plotFlag = varargin{i+1};
        case 'printPlot'
            printPlot = varargin{i+1};
        case 'figureHandle'
            figureHandle = varargin{i+1};
        case 'alignEvent'
            alignEvent = varargin{i+1};
        case 'filterData'
            filterData = varargin{i+1};
        otherwise
    end
end


% Constants
preSaccRejMS    = 20; % Used to reject EEG data 20 ms before a saccade begins
n2pcEpoch       = [125 250];
baseWindow      = -150 : -99;   % relative to start of subject fixating
rejectWindow    = baseWindow(2) : n2pcEpoch(2);
satThreshold    = 50; % How many ms will we allow a signal to be saturated?
displayRange    = [-100, 500];  % relative to
displayRange    = displayRange - 200;
stopHz          = 50;  % What is the cutoff frequency for filterins the signal?

electrodeArray = eeg_electrode_map(subjectID);
nElectrode = length(electrodeArray);

switch lower(subjectID)
    case 'broca'
        lElectrodeName = 'o1';
        rElectrodeName = 'o2';
    case 'xena'
        lElectrodeName = 'f1';
        rElectrodeName = 'f2';
end
lElectrodeChannel = find(strcmp(lElectrodeName, electrodeArray));
rElectrodeChannel = find(strcmp(rElectrodeName, electrodeArray));





% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID, [ccm_min_vars, 'eegData']);
pSignalArray = ExtraVar.pSignalArray;
ssdArray = ExtraVar.ssdArray;

% make sure there are as many electrodes as we expect- otherwise, alert
% user and return with empty data
if size(trialData.trialOutcome.eegData, 2) < nElectrode - 1
    fprintf('%s %s should have %d electrodes but seems to have %d for the session', subjectID, sessionID, nElectrode, size(trialData.trialOutcome.eegData, 1));
    data = [];
    return
end


alignTime = trialData.(alignEvent);












% Set up plots
if plotFlag
    cMap = ccm_colormap(pSignalArray);
    leftColor = cMap(1,:) .* .8;
    rightColor = cMap(end,:) .* .8;
    diffColor = [.4 .4 .4];
    
    nPlotColumn = 2;
    nPlotRow = 2;
    yAxLim = [-.03 .03];
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition]   = standard_figure(nPlotRow, nPlotColumn, 'portrait', figureHandle);
    
    % axes
    lE = 1;
    rE = 2;
    ax(lE) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
    hold(ax(lE), 'on')
    ax(rE) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
    hold(ax(rE), 'on')
end




% *********    RIGHT NOW ONLY DOING GO TARGET TRIALS-- NEED TO EXPAND TO DISTRACTOR TRIALS AND STOP TRIALS


% Perform the anlalyses on go trials to the (correct) target
% Get default trial selection options
selectOpt       = ccm_trial_selection;
selectOpt.ssd       = 'none';
selectOpt.rightCheckerPct = 'collapse';
selectOpt.outcome     = {'goCorrectTarget'};



% Trials to analyze
selectOpt.targDir   	= 'right';
rGoTargTrial = ccm_trial_selection(trialData, selectOpt);
selectOpt.targDir   	= 'left';
lGoTargTrial = ccm_trial_selection(trialData, selectOpt);

rGoTargTrial(isnan(alignTime(rGoTargTrial)) | isnan(trialData.rt(rGoTargTrial))) = [];
lGoTargTrial(isnan(alignTime(lGoTargTrial)) | isnan(trialData.rt(lGoTargTrial))) = [];

% Vectors of alignment times for the EEG data
rAlignTime = alignTime(rGoTargTrial);
lAlignTime = alignTime(lGoTargTrial);

% lSaccTime = trialData.responseOnset(lGoTargTrial);
% rSaccTime = trialData.responseOnset(rGoTargTrial);

nLeftTrial = length(lGoTargTrial);
nRightTrial = length(rGoTargTrial);
%%



% EEG data (aligned on checker onset) from the relevant trials and the correct channel
[o1RightTargEEG, rightTargAlign]    = align_signals(trialData.eegData(rGoTargTrial, lElectrodeChannel), rAlignTime);
[o2RightTargEEG, ~]                 = align_signals(trialData.eegData(rGoTargTrial, rElectrodeChannel), rAlignTime);
[o2LeftTargEEG, leftTargAlign]      = align_signals(trialData.eegData(lGoTargTrial, rElectrodeChannel), lAlignTime);
[o1LeftTargEEG, ~]                  = align_signals(trialData.eegData(lGoTargTrial, lElectrodeChannel), lAlignTime);

%
%%
rtR     = trialData.rt(rGoTargTrial);
rtL     = trialData.rt(lGoTargTrial);
rtMeanR = nanmean(rtR);
rtMeanL = nanmean(rtL);


% Loop through right-target trials
% ----------------------------------------
rSatTrial = zeros(nRightTrial, 1);
for iTrial = 1 : nRightTrial
    
    % Truncate signal just before saccade
    if strcmp(alignEvent, 'checkerOn')
        o1RightTargEEG(iTrial, rightTargAlign + trialData.rt(rGoTargTrial(iTrial)) - preSaccRejMS : end) = nan;
        o2RightTargEEG(iTrial, rightTargAlign + trialData.rt(rGoTargTrial(iTrial)) - preSaccRejMS : end) = nan;
    end
    
    % Reject trials with saturated signals
    if sum(diff(o1RightTargEEG(iTrial, rightTargAlign + rejectWindow)) == 0) > satThreshold || ...
            sum(diff(o2RightTargEEG(iTrial, rightTargAlign + rejectWindow)) == 0) > satThreshold
        rSatTrial(iTrial) = 1;
    end
    
    
    % Baseline-shift the signals
    iBaseO1 = mean(o1RightTargEEG(iTrial, rightTargAlign + baseWindow));
    o1RightTargEEG(iTrial, :)  = o1RightTargEEG(iTrial, :) - iBaseO1;
    iBaseO2 = mean(o2RightTargEEG(iTrial, rightTargAlign + baseWindow));
    o2RightTargEEG(iTrial, :)  = o2RightTargEEG(iTrial, :) - iBaseO2;
    
end %  iTrial = 1 : nRightTrial

% Delete rejected trials
rtR(find(rSatTrial),:)          = nan;
o1RightTargEEG(find(rSatTrial),:) = nan;
o2RightTargEEG(find(rSatTrial),:) = nan;


% Loop through left-target trials
% ----------------------------------------
lSatTrial = zeros(nLeftTrial, 1);
for iTrial = 1 : nLeftTrial
    
    % Truncate signal just before saccade
    if strcmp(alignEvent, 'checkerOn')
        o2LeftTargEEG(iTrial, leftTargAlign + trialData.rt(lGoTargTrial(iTrial)) - preSaccRejMS : end) = nan;
        o1LeftTargEEG(iTrial, leftTargAlign + trialData.rt(lGoTargTrial(iTrial)) - preSaccRejMS : end) = nan;
    end
    
    % Reject trials with saturated signals
    if sum(diff(o2LeftTargEEG(iTrial, leftTargAlign + rejectWindow)) == 0) > satThreshold || ...
            sum(diff(o1LeftTargEEG(iTrial, leftTargAlign + rejectWindow)) == 0) > satThreshold
        lSatTrial(iTrial) = 1;
    end
    
    % Baseline-shift the signals
    iBaseO1 = mean(o2LeftTargEEG(iTrial, leftTargAlign + baseWindow));
    o2LeftTargEEG(iTrial, :)  = o2LeftTargEEG(iTrial, :) - iBaseO1;
    iBaseO2 = mean(o1LeftTargEEG(iTrial, leftTargAlign + baseWindow));
    o1LeftTargEEG(iTrial, :)  = o1LeftTargEEG(iTrial, :) - iBaseO2;
    
end %  iTrial = 1 : nLeftTrial

% Delete rejected trials
rtL(find(lSatTrial))          = nan;
o2LeftTargEEG(find(lSatTrial),:) = nan;
o1LeftTargEEG(find(lSatTrial),:) = nan;





%%

plotTrial = false;
if plotTrial
    clf
    hold all;
    for iTrial = 1 :nLeftTrial
        iTrial
        cla
        xlim([0 800]);
        plot(o1RightTargEEG(iTrial, rightTargAlign - 400 : rightTargAlign + 400), 'b')
        plot(o2RightTargEEG(iTrial, rightTargAlign - 400 : rightTargAlign + 400), 'r')
        pause
    end
    for iTrial = 1 : nRightTrial
        clf
        plot(o2LeftTargEEG(iTrial, leftTargAlign - 400 : leftTargAlign + 400), 'b')
        plot(o1LeftTargEEG(iTrial, leftTargAlign - 400 : leftTargAlign + 400), '--b')
        pause
    end
end



% Extract a segment of the signal to display (and use to make a difference
% waveform)
o1RightTargEEG = o1RightTargEEG(:, rightTargAlign + displayRange(1) : rightTargAlign + displayRange(2));
o2RightTargEEG = o2RightTargEEG(:, rightTargAlign + displayRange(1) : rightTargAlign + displayRange(2));
o2LeftTargEEG = o2LeftTargEEG(:, leftTargAlign + displayRange(1) : leftTargAlign + displayRange(2));
o1LeftTargEEG = o1LeftTargEEG(:, leftTargAlign + displayRange(1) : leftTargAlign + displayRange(2));

% o1RTargMean = o1RTargMean(rightTargAlign + displayRange(1) : rightTargAlign + displayRange(2));
% o2RTargMean = o2RTargMean(rightTargAlign + displayRange(1) : rightTargAlign + displayRange(2));
% o2LTargMean = o2LTargMean(rightTargAlign + displayRange(1) : rightTargAlign + displayRange(2));
% o1LTargMean = o1LTargMean(rightTargAlign + displayRange(1) : rightTargAlign + displayRange(2));



% ****************************************
%   Band-pass the data?
if filterData
    
    o1RTargMean = lowpass(nanmean(o1RightTargEEG, 1)', stopHz);
    o2RTargMean = lowpass(nanmean(o2RightTargEEG, 1)', stopHz);
    o2LTargMean = lowpass(nanmean(o2LeftTargEEG, 1)', stopHz);
    o1LTargMean = lowpass(nanmean(o1LeftTargEEG, 1)', stopHz);
    
else
    o1RTargMean = nanmean(o1RightTargEEG, 1);
    o2RTargMean = nanmean(o2RightTargEEG, 1);
    o2LTargMean = nanmean(o2LeftTargEEG, 1);
    o1LTargMean = nanmean(o1LeftTargEEG, 1);
end
% ****************************************


% Difference waveforms:
o1RightTargIpsiDiff = o1RTargMean - o1LTargMean;
o2LeftTargIpsiDiff  = o2LTargMean - o2RTargMean;
lrDiff              = o1RTargMean - o2RTargMean;
rlDiff              = o2LTargMean - o1LTargMean;

%%
% xMax = 1.1 * max(rtMeanL, rtMeanR);
xMax        = 500;
if plotFlag
    
    plot(ax(lE), o1RTargMean, 'color', rightColor, 'linewidth', 2)
    %     plot(ax(lE), nanmean(o2RightTargEEGAlign, 1), 'b', 'linewidth', 2)
    %     plot(ax(lE), lrDiff, diffColor, 'linewidth', 2)
    plot(ax(lE), o1LTargMean, 'color', leftColor, 'linewidth', 2)
    plot(ax(lE), o1RightTargIpsiDiff, 'color', diffColor, 'linewidth', 2)
    ylim(ax(lE), yAxLim)
    plot(ax(lE), [-displayRange(1) -displayRange(1)], ylim(ax(lE)))
    if strcmp(alignEvent, 'checkerOn')
        plot(ax(lE), [rtMeanR + -displayRange(1)  rtMeanR + -displayRange(1)], ylim(ax(lE)), '--', 'color', rightColor)
        plot(ax(lE), [rtMeanL + -displayRange(1)  rtMeanL + -displayRange(1)], ylim(ax(lE)), '--', 'color', leftColor)
    end
    xlim(ax(lE), [0 xMax])
    set(ax(lE), 'xtick', (0 : 100 : displayRange(2) - displayRange(1)))
    set(ax(lE), 'XTickLabel', [displayRange(1) : 100 : displayRange(2)])
    legend(ax(lE), 'o1 Right (Contra)', 'o1 Left (Ipsi)', 'Contra - Ipsi')
    
    plot(ax(rE), o2LTargMean, 'color', leftColor, 'linewidth', 2)
    %     plot(ax(rE), nanmean(o1LeftTargEEGAlign, 1), 'b', 'linewidth', 2)
    %     plot(ax(rE), rlDiff, diffColor, 'linewidth', 2)
    plot(ax(rE), o2RTargMean, 'color', rightColor, 'linewidth', 2)
    plot(ax(rE), o2LeftTargIpsiDiff, 'color', diffColor, 'linewidth', 2)
    ylim(ax(rE), yAxLim)
    plot(ax(rE), [-displayRange(1) -displayRange(1)], ylim(ax(rE)))
    if strcmp(alignEvent, 'checkerOn')
        plot(ax(rE), [rtMeanR + -displayRange(1)  rtMeanR + -displayRange(1)], ylim(ax(rE)), '--', 'color', rightColor)
        plot(ax(rE), [rtMeanL + -displayRange(1)  rtMeanL + -displayRange(1)], ylim(ax(rE)), '--', 'color', leftColor)
    end
    xlim(ax(rE), [0 xMax])
    set(ax(rE), 'xtick', (0 : 100 : displayRange(2) - displayRange(1)))
    set(ax(rE), 'XTickLabel', [displayRange(1) : 100 : displayRange(2)])
    legend(ax(rE), 'o2 Left (Contra)', 'o2 Right (Ipsi)', 'Contra - Ipsi')
end



data.o1RightTargEEG	= o1RightTargEEG;
data.o2RightTargEEG	= o2RightTargEEG;
data.o2LeftTargEEG  = o2LeftTargEEG;
data.o1LeftTargEEG 	= o1LeftTargEEG;
data.rtR            = rtR;
data.rtL            = rtL;

data.alignTime      = -displayRange(1);
data.displayRange      = displayRange;



end
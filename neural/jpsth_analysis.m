function [jpsthData, unitNames] = jpsth_analysis(monkey, brainArea, sessionNumber, alignmentEvent, timeWindow, outcome, target, plotFlag)

% JPSTH_ANALYSYS

sessionNumberString = num2str(sessionNumber, '%03.0f');

coinWidth = 10;
binWidth = 1;

% dataFileType = 'dataset';
% switch dataFileType
%     case 'dataset'
%         % dataset converted files:
%         dataPath = ['~/Volumes/teba/SchallLab/data/', regexprep(monkey, '(\<[a-z])','${upper($1)}'), '/datasetFiles/'];
%         fileName = [monkey, '_', brainArea, '_', sessionNumberString];
%
%     case 'translated'
%         % original translated files:
%         translatedFlag = 1;
%         recordingSystem = 'merged';
%         dataPath = get_monkey_data_path(monkey, brainArea, translatedFlag, recordingSystem);
%         directory = dir(dataPath);
%
%         foundSession = 0;
%         iSession = 1;
%         while ~foundSession
%             %     iSession
%             %         directory(iSession).name
%             if regexp(directory(iSession).name, sessionNumberString)
%                 foundSession = 1;
%                 fileName = directory(iSession).name;
%             end
%             iSession = iSession + 1;
%         end
% end
%
%
%
% load([dataPath, fileName]);

fileName = [monkey, '_', brainArea, '_', sessionNumberString];
load(fileName);

% Get t
nUnit = size(trialData.trialOutcome.singleUnit, 2);

if nUnit < 2
    fprintf('%s has only one or no single units recorded\n', fileName)
    return
end


% 1. Make a list of the possible 2-neuron combinations to analyze, filtering
% out neurons with little or no spikes
unitCompare = nchoosek(1 : nUnit, 2);
unitNames = cell(size(unitCompare));



% 2. Loop through and get the jpsth data, assigning each pairwise data
% structure a unique name.
for iCompare = 1: size(unitCompare, 1)
    
    [trialList, alignTimeList] = trial_selection(monkey, brainArea, sessionNumber, outcome, alignmentEvent, target);
    nTarget = length(trialList);
    for jTarget = 1 : nTarget
        jTrialList = trialList{jTarget};
        jAlignTimeList = alignTimeList{jTarget};
        % [alignEventVector trialList] = asdf(trialData, alignmentEvent);
        
        singleUnit1 = trialData.singleUnit(:, unitCompare(iCompare, 1));
        singleUnit2 = trialData.singleUnit(:, unitCompare(iCompare, 2));
        nTrial = length(jTrialList);
        unit1Spike = singleUnit1(jTrialList, :);
        unit2Spike = singleUnit2(jTrialList, :);
        maxSpike1 = max(cell2mat(cellfun(@length, unit1Spike, 'uniformoutput', false)));
        maxSpike2 = max(cell2mat(cellfun(@length, unit2Spike, 'uniformoutput', false)));
        spikeData1 = zeros(nTrial, maxSpike1);
        spikeData2 = zeros(nTrial, maxSpike2);
        for iTrial = 1 : nTrial
            spikeData1(iTrial, 1 : length(unit1Spike{iTrial,:})) = unit1Spike{iTrial,:};
            spikeData2(iTrial, 1 : length(unit2Spike{iTrial,:})) = unit2Spike{iTrial,:};
        end
        
        
        alignedSpikeData1 = alignTimeStamps(spikeData1, jAlignTimeList);
        alignedSpikeData2 = alignTimeStamps(spikeData2, jAlignTimeList);
        timeStamps1 = trimTimeStamps(alignedSpikeData1, timeWindow);
        timeStamps2 = trimTimeStamps(alignedSpikeData2, timeWindow);
        spikes1 = spikeCounts(timeStamps1, timeWindow, binWidth);
        spikes2 = spikeCounts(timeStamps2, timeWindow, binWidth);
        
%         sameNueronTest = intersect(
        jpsthData(iCompare, jTarget) = jpsth(spikes1, spikes2, coinWidth);
        unitNames(iCompare, 1:2) = sessionData.collectedData.neurophysiology.singleUnit.names(unitCompare(iCompare, 1:2));
    end
end
% 3. Use Z-scores in jpsth data to determine significance testing




% 4. Display the jpsth-- create separate function
if plotFlag
    jpsth_display(jpsthData, unitNames)
end



% function  [alignTimeList trialList] = asdf(trialData, alignmentEvent)
%
% switch alignmentEvent
%     case 'fixation'
%         allTrialTimes = cell2mat(trialData.fixationWindowEntered);
%     case 'target'
%         allTrialTimes = cell2mat(trialData.targetOnset);
%     case 'stopSignal'
%         allTrialTimes = cell2mat(trialData.stopOnset);
%     case 'saccade'
%         % list of trials with a reponse to target
%         list = ~isnan(cell2mat(trialData.responseToTargetIndex));
%         % remove NaNs from the cell arrays
%         index = trialData.responseToTargetIndex(list);
%         onset = trialData.responseOnset(list);
%         % extract the saccade times desired from the lists of all responses
%         saccadeTimes = cell2mat(cellfun(@(x,y) x(y), onset, index, 'uniformoutput', false));
%         % initialize allTrialTimes as all NaNs
%         allTrialTimes = nan(length(list), 1);
%         % fill in the trials that have responses with time of response
%         allTrialTimes(list) = saccadeTimes;
%     otherwise
%         fprintf('jpsth_analysis --> function asdf:  Don''t have code yet for aligning on %s ', alignmentEvent);
% end
% trialList = find(~isnan(allTrialTimes));
% alignTimeList = allTrialTimes(trialList);
% end






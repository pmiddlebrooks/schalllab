function  [trialList] = gng_trial_selection(subjectID, sessionID, outcomeArray, goCheckerPercentRange)
% function  [trialList, alignTimeList] = ccm_trial_selection(monkey, sessionID, outcomeArray, alignmentEvent, target)

% outcome: 'goCorrect', 'goIncorrect', 'stopCorrect', 'stopIncorrect', 'all'
%
%
%
% target: 'all' to collapse across targets, 'each' to anaylyze them
% separately
%
%
%
%
%
%


% Load the data
[dataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessionID);
% If the file hasn't already been copied to a local directory, do it now
if exist(localDataFile, 'file') ~= 2
    copyfile(dataFile, localDataPath)
end
load(localDataFile);







% Convert cells to doubles if necessary
trialData = cell_to_mat(trialData);





% Get list(s) of trials w.r.t. the outcome
if strcmp(outcomeArray, 'all')
    outcomeList = ones(nTrial, 1);
else
    outcomeList = zeros(nTrial, 1);
    for iOutcomeIndex = 1 : length(outcomeArray)
        iOutcome = outcomeArray{iOutcomeIndex};
        
        outcomeList = outcomeList + strcmp(trialData.trialOutcome, iOutcome);
    end
end
trialLogical = trialLogical & outcomeList;
% sum(trialLogical)
% 



% Get list(s) of trials w.r.t. the checkerboard proportion of target1
% checkers
if strcmp(goCheckerPercentRange, 'all')
    goProportion = (0 : 100) / 100;
else
    goProportion = goCheckerPercentRange / 100;
end
trialLogical = trialLogical & ismember(cell2mat(trialData.goCheckerProportion), goProportion);
sum(trialLogical)




trialList = find(trialLogical);

return








% Get list(s) of trials w.r.t. the targets (targetTrialList) and a list of
% the target locations (targetList) if necessary
switch target
    case 'all'
        nList = 1;
        targetTrialList = {allTrialList(~cellfun(@isnan, trialData.targOn))};
    case 'each'
        targetLoc = unique(cell2mat(trialData.targetLocation), 'rows');
        targetLocations = targetLoc(~isnan(targetLoc(:,1)), :);
        nTarget = size(targetLocations, 1);
        nList = nTarget;
        for iTarget = 1 : nTarget
            targetList{iTarget} = targetLocations(iTarget, :);
            targetTrialList{iTarget} = allTrialList(ismember(cell2mat(trialData.targetLocation), targetLocations(iTarget,:),'rows'));
        end
    otherwise
        disp('trial_selection.m:  "target" input needs to be either "all" or "each"')
        return
end




% Get list(s) of trials w.r.t. the targets
switch alignmentEvent
    case 'responseOnset'
        % list of trials with a reponse to target
        alignTrialList = allTrialList(~isnan(cell2mat(trialData.saccToTargIndex)));
        % remove NaNs from the cell arrays
        index = trialData.saccToTargIndex(alignTrialList);
        onset = trialData.responseOnset(alignTrialList);
        % extract the saccade times desired from the alignTrialLists of all responses
        saccadeTimes = cell2mat(cellfun(@(x,y) x(y), onset, index, 'uniformoutput', false));
        % initialize allTrialTimes as all NaNs
        allTrialTimes = nan(length(alignTrialList), 1);
        % fill in the trials that have responses with time of response
        allTrialTimes(alignTrialList) = saccadeTimes;
    case 'rewardOnset'
        trialData.rewardOnset = cell2mat(trialData.rewardOnset);
        alignTrialList = allTrialList(~isnan(trialData.rewardOnset(:, 1)));
        allTrialTimes = trialData.rewardOnset(:, 1);
    otherwise
        allTrialTimes = cell2mat(trialData.(alignmentEvent));
        alignTrialList = allTrialList(~isnan(cell2mat(trialData.(alignmentEvent))));
end


for iList = 1 : nList
%     trialList{iList} = allTrialList(outcomeList & targetList{iList} & alignTrialList);
%     alignTimeList{iList} = allTrialTimes(outcomeList & targetList{iList} & alignTrialList);
    trialList{iList} = intersect(outcomeTrialList, intersect(targetTrialList{iList}, alignTrialList))';
    alignTimeList{iList} = allTrialTimes(trialList{iList});
    
    
    % trialList = find(~isnan(allTrialTimes));
    % alignTimeList = allTrialTimes(trialList);
end



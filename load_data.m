function [trialData, SessionData, ExtraVariable] = load_data(subjectID, sessionID, variables, multiUnit)
% function [trialData, SessionData] = load_data(subjectID, sessionID)
%
% Loads a data file and does some minimal processing common to lots of
% analyses
%

if nargin < 4
    multiUnit = false;
end


ExtraVariable = struct();

if ismember(lower(subjectID), {'joule', 'broca', 'xena', 'chase', 'hoagie', 'norm', 'andy', 'nebby', 'shuffles'})
    monkeyOrHuman = 'monkey';
else
    monkeyOrHuman = 'human';
end


% Load the data
if strcmp(sessionID(end-3:end), '.mat')
    sessionID(end-3 : end) = [];
end

[dataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessionID, monkeyOrHuman);

% If the file hasn't already been copied to a local directory, do it now
if exist(localDataFile, 'file') ~= 2
    copyfile(dataFile, localDataPath)
    disp(sessionID)
end

if nargin > 2
trialData = load(localDataFile, variables{:});
else
    trialData = load(localDataFile);
end
SessionData = trialData.SessionData;
trialData = rmfield(trialData, 'SessionData');

if multiUnit
    if ismember('spikeData', variables)
    [SessionData.spikeUnitArray, trialData.spikeData] = convert_to_multiunit(SessionData.spikeUnitArray, trialData.spikeData);
    else
    [SessionData.spikeUnitArray, ~] = convert_to_multiunit(SessionData.spikeUnitArray);
    end
end

if isfield(SessionData, 'taskID')
    task = SessionData.taskID;
elseif isfield(SessionData.task, 'taskID')
    task = SessionData.task.taskID;
end


% Convert cells to doubles if necessary
if ~strcmp(task, 'maskbet')
    trialData = cell_to_mat(trialData);
    trialData.iTrial = (1 : size(trialData.trialOutcome,1))';
end








if strcmp(task, 'ccm')
    %    fieldnames(trialData)';
    pSignalArray = unique(trialData.targ1CheckerProp);
    pSignalArray(isnan(pSignalArray)) = [];
    ExtraVariable.pSignalArray = pSignalArray;
end




if strcmp(task, 'ccm') || strcmp(task, 'cmd')
    % Need to do a little SSD value adjusting, due to ms difference and 1-frame
    % differences in SSD values
    ssd = trialData.stopSignalOn - trialData.responseCueOn;
    
%     %    Old method
%            ssdArray = unique(trialData.stopSignalOn - trialData.responseCueOn);
%            ssdArray(isnan(ssdArray)) = [];
%            if ~isempty(ssdArray)
%                diffExist = 1;
%                while diffExist
%                    a = diff(ssdArray);
%                    addOne = ssdArray(a == 1);
%                    [d,i] = ismember(ssd, addOne);
%                    if sum(d) == 0
%                        diffExist = 0;
%                    else
%                        ssd(d) = ssd(d) + 1;
%                        % ssdArray(a == 1) = ssdArray(a == 1) + 1;
%                        ssdArray = unique(ssd(~isnan(ssd)));
%                    end
%                end
%            end
%            trialData.ssd = ssd;
%            ExtraVariable.ssdArray = ssdArray;
    
    % New method (3/18/14)
    trialData.ssd = ssd_session_adjust(ssd);
    ExtraVariable.ssdArray = unique(trialData.ssd(~isnan(trialData.ssd)));
end


% Some of Xena's early sessions, the TEMPO code was outputing the wrong
% targ1CheckerProp if it was supposed to be .53 (it was outputing .52 due
% to how I set up the TEMPO code to calculate it). Fix that here
switch subjectID
    case 'xena'
        if ismember(.52, trialData.targ1CheckerProp) && ismember(.47, trialData.targ1CheckerProp)
            trialData.targ1CheckerProp(trialData.targ1CheckerProp == .52) = .53;
        end
    otherwise
end


if ~strcmp(task, 'maskbet')
    trialData.rt = trialData.responseOnset - trialData.responseCueOn;
    
    if strcmp(task, 'ccm')
        % If there isn't a distractor angle variable, assume distractor is
        % 180 degrees from target
        if ~ismember('distAngle', fieldnames(trialData))
            trialData.distAngle = trialData.targAngle + 180;
        end
        angleMat = unique([trialData.targAngle trialData.distAngle], 'rows');
        ExtraVariable.targAngleArray = angleMat(:,1);
        ExtraVariable.distAngleArray = angleMat(:,2);
    else
        ExtraVariable.targAngleArray = unique(trialData.targAngle(~isnan(trialData.targAngle)));
    end
    
    % Want to get rid of trials that were mistakenly recored with targets at
    % the wrong angles (happens sometimes at the beginning of a task session
    % recording when the angles were set wrong). For now, sue the criteria that
    % a target must have at ewast 7 trials to considered legitimate
    lowNTarg = zeros(size(trialData.trialOutcome, 1), 1);
    for i = 1 : length(ExtraVariable.targAngleArray)
        iTrial = trialData.targAngle == ExtraVariable.targAngleArray(i);
        if sum(iTrial) <= 7
            lowNTarg(iTrial) = 1;
        end
    end
    if sum(lowNTarg)
    trialData = structfun(@(x) x(~logical(lowNTarg), :), trialData, 'uni', false);
    end
else
    %     trialData.decRT = trialData.decResponseOnset - trialData.decResponseCueOn;
    %     trialData.betRT = trialData.betResponseOnset - trialData.betResponseCueOn;
    ExtraVariable.soaArray = unique(trialData.soa(~isnan(trialData.soa)));
end
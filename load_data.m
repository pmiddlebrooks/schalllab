function [trialData, SessionData, ExtraVariable] = load_data(subjectID, sessionID, variables, multiUnit)
% function [trialData, SessionData] = load_data(subjectID, sessionID)
%
% Loads a data file and does some minimal processing common to lots of
% analyses
%
% To load all spike units, input "spikeData" as part of the variables list.

warning('off','all')

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
    if ~isempty(dataFile)
        copyfile(dataFile, localDataPath)
        disp(sessionID)
    else
        error(sprintf('You don"t seem to have access to that data file from your local computer\n'))
    end
end


% Load the SessionData information to use it if necessary
load(localDataFile, 'SessionData');



% Do we want to load all the data, or specific variables?
if nargin < 3
    trialData = load(localDataFile);
else
    % If asking for all the lfp Data, add all the lfp channels to the list
    % of variables
    if ismember('lfpData', variables)
        lfpChannel = SessionData.lfpChannel;
        lfpVars = {};
        for i = 1 : length(lfpChannel)
            lfpVars = [lfpVars, sprintf('lfp%s', num2str(lfpChannel(i), '%02i'))];
        end
        variables = [variables, lfpVars];
    end
    
    % If user is not asking for spikeUnits/spikeData or lfp data, simply load whatever
    % variables requested. Otherwise some special processing may be
    % required depending on nature of spike data requested
    if ~cell2mat(regexp(variables, 'spike.*'))
        trialData = load(localDataFile, variables{:});
    else
        if cell2mat(regexp(variables, 'spike.*'))
            % If user doesn't request multiUnit, load all or individual units
            % requested.
            if ~multiUnit
                if ismember('spikeData', variables)
                    spikeUnitArray = SessionData.spikeUnitArray;
                    variables = [variables, spikeUnitArray];
                end
                trialData = load(localDataFile, variables{:});
            else % if user requests multiUnit...
                % If user requests all spikeData, load list of all units and
                % figure out unique set of names without individual units
                % (spikeUnit01, 02, etc, without 01a, 01b, etc).
                if ismember('spikeData', variables)
                    unitList = unique(cellfun(@(x) x(1:11), SessionData.spikeUnitArray, 'uni', false));
                else
                    unitList = variables(strncmp(variables, 'spikeUnit', 8));
                end
                
                % Once we have list of units, load each set of individual units
                % separately and create multi-units out of each set, adding
                % each multiUnit to the trialData struct.
                % remove the units from the variables list and load all
                % non-units:
                variables(strncmp(variables, 'spikeUnit', 8)) = [];
                trialData = load(localDataFile, variables{:});
                
                % Load each set of units separately, then concatenate them into
                % multiunit, adding each multiUnit to trialData.
                for j = 1 : length(unitList)
                    spikeData = load(localDataFile, '-regexp', sprintf('%s.*', unitList{j}));
                    jNames = fieldnames(spikeData);
                    tempCell = cell(size(spikeData.(jNames{1}), 1), 1);
                    for k = 1 : length(jNames)
                        spikeData.(jNames{k}) = cellfun(@(x) reshape(x, [], length(x)), spikeData.(jNames{k}), 'uni', false);
                        
                        tempCell = cellfun(@(x,y) [x,y], tempCell, spikeData.(jNames{k}), 'uni', false);
                        %                 trialData.(unitList{j}) = arrayfun(@(IDX) [spikeData{IDX,:}], 1:size(spikeData,1), 'Uniform', 0);
                    end
                    trialData.(unitList{j}) = tempCell;
                end
            end
            trialData.SessionData.spikeUnitArray = unique(cellfun(@(x) x(1:11), SessionData.spikeUnitArray, 'uni', false));
        end
    end
end
SessionData = trialData.SessionData;
trialData = rmfield(trialData, 'SessionData');



% if multiUnit
%     if ismember('spikeData', variables)
%         [SessionData.spikeUnitArray, trialData.spikeData] = convert_to_multiunit(SessionData.spikeUnitArray, trialData.spikeData);
%     else
%         [SessionData.spikeUnitArray, ~] = convert_to_multiunit(SessionData.spikeUnitArray);
%     end
% end

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
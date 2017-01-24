function concatenate_human_sessions(sessionID, concatID, eyeOrKey)

% dataFolder = 'data/';
humanDataPath = '/Volumes/middlepg/HumanData/ChoiceStopTask/';

% humanDataPath = '~/schalllab/';
switch eyeOrKey
    case 'eye'
        effectorTag = 'em';
        effector = 'saccade';
        %         effector = 'saccadeMod';
    case 'key'
        effectorTag = 'kp';
        effector = 'keypress';
end
[tebDataFile, localDataPath, localDataFile] = data_file_path(sessionID(1:2), [sessionID(3:end), effector]);

% Load single session data
sessionDataFile = [humanDataPath, sessionID, effector, '.mat'];
load(sessionDataFile);
clear SessionData
SessionData.taskID = 'ccm';

concatData = dataset(...
    {trialData.trialOutcome, 'trialOutcome'}, ...
    {trialData.responseOnset, 'responseOnset'},...
    {trialData.targOn, 'targOn'},...
    {trialData.targ1CheckerProp, 'targ1CheckerProp'},...
    {trialData.stopSignalOn, 'stopSignalOn'},...
    {trialData.targAngle, 'targAngle'},...
    {trialData.saccAngle, 'saccAngle'},...
    {trialData.responseCueOn, 'responseCueOn'});

% If this is the first file we're adding to a concatenated file (within a single subject), just use
% the sessionID data
if isempty(concatID)
    concatID = [sessionID(1:2), 'All'];
    concatDataFile = [humanDataPath, concatID, effector, '.mat'];
    %     concatDataFile = [humanDataPath, concatID, 'saccade', '.mat'];
    sessionIDArray = {[sessionID, effector]};
    trialData = concatData;
    % Otherwise add the sessionID data to the concatenated file
    save(concatDataFile, 'trialData', 'SessionData', 'sessionIDArray')
else
    if strcmp(concatID, 'hu')
        % If it's the first across-subjects concatenation
        concatDataFile = [humanDataPath, concatID, 'All', effector, '.mat'];
        
        trialData = concatData;
        sessionIDArray = [sessionID, effector];
        save(concatDataFile, 'trialData', 'SessionData', 'sessionIDArray')
        
    else
        
        % Load concatenated data
        concatDataFile = [humanDataPath, concatID, effector, '.mat'];
        load(concatDataFile)
        
        trialData = [trialData; concatData];
        sessionIDArray = [sessionIDArray; [sessionID, effector]];
        save(concatDataFile, 'trialData', 'SessionData', 'sessionIDArray', '-append')
    end
end


copyfile(concatDataFile, localDataPath)

return

%% MANUALLY CONCATENATE SESSIONS

% Saccade

sArray = {'bz0907saccade','bz0924saccade','bb0924saccade'};
sArray = {'pg0928saccade','pg1001saccade','pm1002saccade'};

sArray = {'dm0726saccade','dm0731saccade'};
sArray = {'ts0425saccade','tn0729saccade'};
sArray = {'oe0712saccade','oe0717saccade','og0717saccade'};
sArray = {'xb0723saccade', 'xb0724saccade'};
sArray = {'kf0226saccade','kf0301saccade'};
sArray = {'cb1001saccade','cb1002saccade'};


nSession = length(sArray);

humanDataPath = '/Volumes/middlepg/HumanData/ChoiceStopTask/';

trialData = cell2table({});

for i = 1 : nSession
    iSession = sArray{i}
    iTD = load([humanDataPath,iSession]);
    td = iTD.trialData;
    td.rt = td.responseOnset - td.responseCueOn;
    
    MIN_RT = 120;
    MAX_RT = 1200;
    nSTD   = 3;
[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTD);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);
    td = dataset2table(td);
trialData = [trialData; td];

end
SessionData = iTD.SessionData;
SessionData.taskID = 'ccm';
save(['~/schalllab/local_data/human/', sArray{2}(1:2), 'Allsaccade'], 'trialData', 'SessionData')





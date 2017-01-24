function ccm_collapse_sessions(sessionIDArray, subjectID)

% dataFolder = 'data/';
tebaDataPath = '/Volumes/SchallLab/data/';
humanDataPath = '/Volumes/middlepg/HumanData/ChoiceStopTask/';
localDataPath = '/Users/paulmiddlebrooks/matlab/local_data/';
% humanDataPath = '~/schalllab/';

tempTrialData = [];
switch subjectID
    case 'Broca'
        for iSession = 1 : length(sessionIDArray)
            sessionID = sessionIDArray{iSession};
            
            MatlabData = dir([tebaDataPath, sessionID, '.mat']);
            
            matlabFile = [tebaDataPath, MatlabData.name];
            load(matlabFile)
            
            if ~strcmp(SessionData.taskID, 'ccm')
                fprintf('%s: Not a choice countermanding session, try again', sessionID)
                return
            end
            
            tempTrialData = [tempTrialData; trialData];
        end
        
        trialData = tempTrialData;
        
        collapsedFileName = [humanDataPath, subjectID, 'All', sessionID(6:end)];
        save(collapsedFileName, 'trialData', 'SessionData', 'sessionIDArray');
        disp(collapsedFileName)
    otherwise
        for iSession = 1 : length(sessionIDArray)
            sessionIDArray
            sessionID = sessionIDArray{iSession}

            MatlabData = dir([humanDataPath, sessionID, '.mat']);

            matlabFile = [humanDataPath, MatlabData.name]
            load(matlabFile, 'trialData', 'SessionData')
            
            if ~strcmp(SessionData.taskID, 'ccm')
                fprintf('%s: Not a choice countermanding session, try again', sessionID)
                return
            end
            
            tempTrialData = [tempTrialData; trialData];
        end
        
        trialData = tempTrialData;
        
        collapsedFileName = [humanDataPath, subjectID, 'All', sessionID(6:end)]
        save(collapsedFileName, 'trialData', 'SessionData', 'sessionIDArray');
        disp(collapsedFileName)
end


% 
% tempTrialData = [];
% for iSession = 1 : length(sessionIDArray)
%     sessionID = sessionIDArray{iSession};
%     
%     MatlabData = dir([humanDataPath, sessionID, '.mat']);
%     
%     matlabFile = [humanDataPath, MatlabData.name];
%     load(matlabFile)
%     
%     if ~strcmp(SessionData.taskID, 'ccm')
%         fprintf('%s: Not a choice countermanding session, try again', sessionID)
%         return
%     end
%     
%     tempTrialData = [tempTrialData; trialData];
% end
% 
% trialData = tempTrialData;
% 
% collapsedFileName = [humanDataPath, subjectID, 'All', sessionID(7:end)];
% save(collapsedFileName, 'trialData', 'SessionData', 'sessionIDArray');
% disp(collapsedFileName)
% 

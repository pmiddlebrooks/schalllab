function [tebaDataPath, localDataPath] = subject_data_path(subjectID)


% location = 'home';
location = 'work';


tebaPath = '/Volumes/SchallLab/data/';
humanPath = '/Volumes/middlepg/HumanData/ChoiceStopTask/';

switch location
    case 'work'
        localDataPath = '/Users/paulmiddlebrooks/matlab/local_data/';
        
        switch lower(subjectID)
            case 'broca'
                animal = 'monkey';
                tebaDataPath = [tebaPath, 'Broca/'];
            case 'xena'
                animal = 'monkey';
                tebaDataPath = [tebaPath, 'Xena/Plexon/'];
            otherwise
                animal = 'human';
                %                 fileName = [MatlabData.name];
                tebaDataPath = humanPath;
        end
        
    case 'home'
        localDataPath = '/Users/paulmiddlebrooks/matlab/local_data/';
        switch upper(subjectID)
            case 'broca'
                animal = 'monkey';
                fileName = [sessionID, '.mat'];
                dataFile = [];
                localDataFile = [localDataPath, fileName];
            case 'xena'
                animal = 'monkey';
                fileName = [sessionID, '.mat'];
                dataFile = [];
                localDataFile = [localDataPath, fileName];
            otherwise
                animal = 'human';
                %                 MatlabData = dir([humanDataPath, subjectID, sessionID, '.mat']);
                fileName = [subjectID, sessionID];
                dataFile = [];
                localDataFile = [localDataPath, fileName, '.mat'];
        end
end
function sessionList = find_eeg_sessions(subjectID, taskID)

%%

sessionList = {};


[tebaDataPath, localDataPath] = subject_data_path(subjectID);
% tebaDataPath = '/Volumes/SchallLab/data/';
humanDataPath = '/Volumes/middlepg/HumanData/ChoiceStopTask/';
monkeyDir = [tebaDataPath];
directory = dir(monkeyDir);


i = 34;
firstIndex = i;



if firstIndex < length(directory)
   exampleFile = 'bp046n01.mat';
   nameLength = length(exampleFile);
   
   
   for i = firstIndex : length(directory)
      fprintf('%s: %d\n', directory(i).name, i)
      % If it's a valid session to add to the list, do it here
      if length(directory(i).name) == nameLength && ...           % if it's a valid length file name (adheres to my naming principle)
            strncmpi(directory(i).name, [directory(i).name(1), 'p'], 2) && ...    % if the name of the file is xpxxx...
            str2double(directory(i).name(3:5)) < 900 && ...         % if the part number is less than 900 (some test files have parts 999, e.g.
            strcmp(directory(i).name(end-3:end), '.mat')            % if it's a matlab file
         
         
         % Load the data
         [trialData, SessionData, ExtraVar] = load_data(subjectID, directory(i).name(1:end-4));
         
         
         iTaskID = SessionData.taskID;
         if strcmpi(taskID, iTaskID)
            
            
            % Check for eegData
            
            
            [a, b] = ismember('eegData', fieldnames(trialData));
            if a
               fprintf('Session %s added to list\n', directory(i).name)
               disp(directory(i).name)
               sessionList = [sessionList; directory(i).name];
            end
            
         end
         
      end
   end
   
   
   
   %    save(sessionsFile, 'sessions', '-append')
end



%
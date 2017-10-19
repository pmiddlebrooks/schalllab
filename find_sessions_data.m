function sessionList = find_sessions_data(subjectID, taskID, dataIncludeArray, dataExcludeArray)

%%

% dataIncludeArray: e.g. {'lfp', 'spike'} would only return sessions that
% have at least lfps and neurons recorded during the session
% dataExcludeArray: e.g. {'eeg'} would only return sessions that don't
% have eegs recorded during the session

if nargin < 4
   dataExcludeArray = {};
end

sessionList = {};


[tebaDataPath, localDataPath] = subject_data_path(subjectID);
% tebaDataPath = '/Volumes/SchallLab/data/';
humanDataPath = '/Volumes/middlepg/HumanData/ChoiceStopTask/';
monkeyDir = [tebaDataPath];
directory = dir(monkeyDir);


% Figure out the best first index in the directory to begin looking
i = [];
switch lower(subjectID)
   case 'broca'
      if sum(ismember('eeg', dataIncludeArray))
         i = [i 398];  % start at a reasonable eeg collection, bp036n01
      end
      if sum(ismember('lfp', dataIncludeArray))
         i = [i 712];  % start at first lfp collection session, bp085n02
      end
      if sum(ismember('spike', dataIncludeArray))
%          i = [i 685];  % start at first single neuron collection session, bp075n01
         i = [i 750];  % start at first single neuron collection session, bp075n01
      end
   case 'xena'
      disp('Don"t have xena sessions added yet')
end



firstIndex = min(i);



if firstIndex < length(directory)
   exampleFile = 'bp046n01.mat';
   nameLength = length(exampleFile);
   
   
   for i = firstIndex : length(directory)
      %       fprintf('%s: %d\n', directory(i).name, i)
          fprintf('index %d, session %s\n', i, (directory(i).name(1:end-4)))
     
      % Check whether it's a valid file to consider,
      if (length(directory(i).name) == nameLength ||length(directory(i).name) == nameLength+3)  && ...           % if it's a valid length file name (adheres to my naming principle)
            strncmpi(directory(i).name, [directory(i).name(1), 'p'], 2) && ...    % if the name of the file is xpxxx...
            str2double(directory(i).name(3:5)) < 970 && ...         % if the part number is less than 900 (some test files have parts 999, e.g.
            strcmp(directory(i).name(end-3:end), '.mat')            % if it's a matlab file
         
         % Load the data
         [trialData, SessionData] = load_data(subjectID, directory(i).name(1:end-4));
         
         % Is it the right task?
         iTaskID = SessionData.taskID;
         if strcmpi(taskID, iTaskID)
            
            % Reset flags each time
            eegFlag = 0;
            lfpFlag = 0;
            spikeFlag = 0;
            
            % Are we looking for eeg data?
            if sum(ismember('eeg', dataIncludeArray)) || sum(ismember('eeg', dataExcludeArray))
               % Check for eegData
               [a, b] = ismember('eegData', fieldnames(trialData));
               if a
                  eegFlag = 1;
               end
            end
            
            % Are we looking for lfp data?
            if sum(ismember('lfp', dataIncludeArray)) || sum(ismember('lfp', dataExcludeArray))
               % Check for eegData
               [a, b] = ismember('lfpData', fieldnames(trialData));
               if a
                  lfpFlag = 1;
               end
            end
            
            % Are we looking for spike data?
            if sum(ismember('spike', dataIncludeArray)) || sum(ismember('spike', dataExcludeArray))
               % Check for eegData
               if sum(strncmp(fieldnames(trialData), 'spikeUnit', 9))
                  spikeFlag = 1;
               end
            end
            
            
            
            
            if eegFlag + lfpFlag + spikeFlag == length(dataIncludeArray)
               % If it meets all criteria, add it to the list
               fprintf('Session %s added to list\n', directory(i).name)
               %                disp(directory(i).name)
               sessionList = [sessionList; directory(i).name(1:end-4)];
            end
         end % iTaskID = SessionData.taskID;
      end  % Is it a valid file to consider?
   end
   
   
   
   %    save(sessionsFile, 'sessions', '-append')
end



%
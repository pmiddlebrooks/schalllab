function find_task_sessions(subjectID, taskID)


%%
% taskID = 'vis'
% subjectID = 'Broca'
% location = 'home';
location = 'work';

% % Make sure we're looking for a task that exists
validTasks = {'amp', 'vis', 'del', 'mem', 'cmd', 'ccm', 'gng', 'mcm'};
% validTasks = {'ccm', 'gng', 'mcm'};


% mArray = {'Xena', 'Broca'};
% for m = 1 : length(mArray)
%     subjectID = mArray{m}
%   taskID = 'ccm';


%     for t = 1 : length(validTasks)
%         taskID = validTasks{t}



if ~ismember(taskID, validTasks)
   fprintf('%s is not a valid task ID: choose one of the following:\n\n', taskID)
   disp(validTasks)
   return
end



% Open the current sessions file, or create one if it doesn't exist
sessionsFile = [subjectID, '_sessions.mat'];
if exist(sessionsFile, 'file') ~= 2
   sessions = struct();
   save(sessionsFile, 'sessions', '-mat')
end
load(sessionsFile);




[tebaDataPath, localDataPath] = subject_data_path(subjectID);
% tebaDataPath = '/Volumes/SchallLab/data/';
humanDataPath = '/Volumes/middlepg/HumanData/ChoiceStopTask/';
monkeyDir = [tebaDataPath];
directory = dir(monkeyDir);


% If there's already a session list for that task, we'll add to it. Else
% make a new one from scratch
i = 1;
if isfield(sessions, taskID)
   if ~isempty(sessions.(taskID))
      
      found = 0;
      lastStopFound = 0;
      lastNoStopFound = 0;
      while ~found && i <= length(directory)
         if ~strcmp(taskID, 'ccm')
            if strcmp(sessions.(taskID)(end).name, directory(i).name(1:end-4))
               found = 1;
            end
         elseif strcmp(taskID, 'ccm')
            if strcmp(sessions.(taskID).stop(end).name, directory(i).name(1:end-4))
               lastStopFound = 1;
            end
            if strcmp(sessions.(taskID).noStop(end).name, directory(i).name(1:end-4))
               lastNoStopFound = 1;
            end
            if lastStopFound && lastNoStopFound
               found = 1;
            end
         end
         i = i + 1;
      end
   end
   
end

firstIndex = i;



if firstIndex < length(directory)
   exampleFile = 'bp046n01.mat';
   nameLength = length(exampleFile);
   iCt = size(sessions.(taskID), 2) + 1;
   if isfield(sessions.ccm, 'stop')
      iStop = size(sessions.ccm.stop, 2) + 1;
      iNoStop = size(sessions.ccm.noStop, 2) + 1;
   else
      iStop = 1;
      iNoStop = 1;
   end
   for i = firstIndex : length(sessions)
      
      
      
      
      %         % Go through the list to see if any session names match the
      %         % current directory file. If it's not on the list, consider below
      %         % whether to add it. If it's already on the list, skip it and move
      %         % on
      %         onList = false;
      %         j = i;
      %         while ~onList && j <= length(directory)
      %             if ~strcmp(taskID, 'ccm')
      %                 if strcmp(sessions.(taskID)(j).name, directory(i).name(1:end-4))
      %                     onList = false;
      %                 end
      %             elseif strcmp(taskID, 'ccm')
      %                 if strcmp(sessions.(taskID).stop(j).name, directory(i).name(1:end-4)) || ...
      %                         strcmp(sessions.(taskID).noStop(end).name, directory(i).name(1:end-4))
      %                     onList = 1;
      %                 end
      %             end
      %             j = j + 1;
      %         end
      %
      %         if onList
      %             continue
      %         end
      
      
      
      % If it's a valid session to add to the list, do it here
      if length(directory(i).name) == nameLength && ...           % if it's a valid length file name (adheres to my naming principle)
            strncmpi(directory(i).name, [directory(i).name(1), 'p'], 2) && ...    % if the name of the file is xpxxx...
            str2double(directory(i).name(3:5)) < 900 && ...         % if the part number is less than 900 (some test files have parts 999, e.g.
            strcmp(directory(i).name(end-3:end), '.mat')            % if it's a matlab file
         
         
         % Load the data
         [dataFile, localDataPath, localDataFile] = data_file_path(subjectID, directory(i).name(1:end-4));
         % If the file hasn't already been copied to a local directory, do it now
         if exist(localDataFile, 'file') ~= 2
            copyfile(dataFile, localDataPath)
         end
         load(localDataFile);
         
         %             trialData = cell_to_mat(trialData);
         iTaskID = SessionData.taskID;
         if strcmpi(taskID, iTaskID)
            if ~strcmp(taskID, 'ccm')
               sessions.(taskID)(iCt).name = directory(i).name(1:end-4);
               iCt = iCt + 1;
               % If it's choice countermanding, also distinguish
               % between sessions with and without stopping:
            elseif strcmp(taskID, 'ccm')
               directory(i).name
               
               if length(unique(trialData.targ1CheckerProp)) > 1
                  if sum(~isnan(trialData.stopSignalOn))
                     sessions.(taskID).stop(iStop).name = directory(i).name(1:end-4);
                     directory(i).name(1:end-4)
                     iStop = iStop + 1;
                  else
                     sessions.(taskID).noStop(iNoStop).name = directory(i).name(1:end-4);
                     directory(i).name(1:end-4)
                     iNoStop = iNoStop + 1;
                  end
               end
            end
            
         end
         
      end
   end
   
   
   
   save(sessionsFile, 'sessions', '-append')
end



%     end
% end


sessions_info(subjectID, taskID)

return







function sessions_info(subjectID, taskID)

sessionsFile = [subjectID, '_sessions.mat'];
load(sessionsFile)

switch taskID
   case 'ccm'
      fprintf('Choice Countermanding Sessions and info:\n\n')
      
      fprintf('\n\nSessions WITHOUT Stopping\n')
      fprintf('Session \t Trials \tPsych Fn\n')
      for i = 1 : length(sessions.ccm.noStop)
         % If session info has not been generated for this session, do
         % that and save it to the sessions.mat file. Otherwise skip
         % this and load the info that has been previously saved:
         if ~isfield(sessions.ccm.noStop, 'nTrial')
            sessions.ccm.noStop(1).nTrial = [];
         end
         if isempty(sessions.ccm.noStop(i).nTrial)
            % Load the data
            [dataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessions.ccm.noStop(i).name);
            load(localDataFile);
            
            
            selectOpt = ccm_trial_selection;  % Default trial selection option values
            selectOpt.rightCheckerPct = 'collapse';
            selectOpt.ssd = 'none';
            selectOpt.outcome     = {'goCorrectTarget', 'goCorrectDistractor'};
            
            goTrial = ccm_trial_selection(trialData, selectOpt);
            nTrial = length(goTrial);
            data = ccm_psychometric(subjectID, sessions.ccm.noStop(i).name, 'plotFlag', 0);
            goPsychFn = data.goPsychFn;
            %                 [goPsychFn(1),goPsychFn(end)]
            
            sessions.ccm.noStop(i).nTrial     = nTrial;
            sessions.ccm.noStop(i).goPsychFn  = goPsychFn;
            
         else
            nTrial      = sessions.ccm.noStop(i).nTrial;
            goPsychFn   = sessions.ccm.noStop(i).goPsychFn;
            %                 [goPsychFn(1),goPsychFn(end)]
         end
         if nTrial > 80
            criteria = nTrial >= 100 && goPsychFn(1) <= .2 && goPsychFn(end) >= .8;
            if criteria
               fprintf('%s \t%d \t \t%.2f - %.2f \t*****\n', sessions.ccm.noStop(i).name, nTrial, goPsychFn(1), goPsychFn(end))
            else
               fprintf('%s \t%d \t \t%.2f - %.2f\n', sessions.ccm.noStop(i).name, nTrial, goPsychFn(1), goPsychFn(end))
            end
         end
      end
      
      
      fprintf('Sessions WITH Stopping\n')
      fprintf('Session \tTrials \tp(Stoptrials) \tPsych Fn \tInh Fn\n')
      if ~isfield(sessions.ccm.stop, 'nTrial')
         sessions.ccm.stop(1).nTrial = [];
      end
      for i = 1 : length(sessions.ccm.stop)
         
         % If session info has not been generated for this session, do
         % that and save it to the sessions.mat file. Otherwise skip
         % this and load the info that has been previously saved:
         if isempty(sessions.ccm.stop(i).nTrial)
            % Load the data
            %                 [dataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessions.ccm.stop(i).name);
            %                 load(localDataFile);
            [trialData, SessionData, ExtraVar] = load_data(subjectID, sessions.ccm.stop(i).name);
            
            selectOpt = ccm_trial_selection;  % Default trial selection option values
            selectOpt.rightCheckerPct = 'collapse';
            selectOpt.ssd = 'none';
            selectOpt.outcome     = {'goCorrectTarget', 'goCorrectDistractor'};
            goTrial = ccm_trial_selection(trialData, selectOpt);
            
            selectOpt.ssd = 'collapse';
            selectOpt.outcome     = {'stopCorrect', 'stopIncorrectTarget', 'stopIncorrectDistractor', ...
               'targetHoldAbort', 'stopIncorrectPreSSDTarget', ...
               'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
            stopTrial = ccm_trial_selection(trialData, selectOpt);
            
            nTrial      = length(stopTrial) + length(goTrial);
            nStop       = length(stopTrial);
            stopPct     = length(stopTrial) / nTrial;
            data        = ccm_psychometric(subjectID, sessions.ccm.stop(i).name, 'plotFlag', 0);
            goPsychFn   = data.goPsychFn;
            data        = ccm_inhibition(subjectID, sessions.ccm.stop(i).name, 'plotFlag', 0);
            inhFn       = data.inhibitionFnGrand;
            
            sessions.ccm.stop(i).nTrial     = nTrial;
            sessions.ccm.stop(i).nStop      = nStop;
            sessions.ccm.stop(i).stopPct    = stopPct;
            sessions.ccm.stop(i).goPsychFn  = goPsychFn;
            sessions.ccm.stop(i).inhFn      = inhFn;
         else
            nTrial      = sessions.ccm.stop(i).nTrial;
            nStop       = sessions.ccm.stop(i).nStop;
            stopPct     = sessions.ccm.stop(i).stopPct;
            goPsychFn   = sessions.ccm.stop(i).goPsychFn;
            inhFn       = sessions.ccm.stop(i).inhFn;
            
         end
         if nTrial > 400
            criteria = nStop >= 280 && goPsychFn(1) <= .2 && goPsychFn(end) >= .8 && inhFn(1) <= .1 && inhFn(end) >= .9;
            if criteria
               fprintf('%s \t%d \t%.2f \t\t%.2f - %.2f \t%.2f - %.2f \t*******\n', sessions.ccm.stop(i).name, nTrial, stopPct, goPsychFn(1), goPsychFn(end), inhFn(1), inhFn(end))
            else
               fprintf('%s \t%d \t%.2f \t\t%.2f - %.2f \t%.2f - %.2f\n', sessions.ccm.stop(i).name, nTrial, stopPct, goPsychFn(1), goPsychFn(end), inhFn(1), inhFn(end))
            end
         end
      end
      
   otherwise
      fprintf('Don"t have task info code yet for %s\n', taskID)
end

save(sessionsFile, 'sessions', '-append')

return

%%




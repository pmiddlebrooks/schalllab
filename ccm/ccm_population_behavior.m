function populationData = ccm_population_behavior

%%
subjectID = ...
    {'bz', ...
    'pm', ...
    'cb'};
sessionArray = ...
    {'Allsaccade', ...
    'Allsaccade', ...
    'Allsaccade'};


iPlotFlag = 1;
nSession = length(sessionArray);

populationData = [];
for iSession = 1 : nSession
    
    iSessionID = sessionArray{iSession};
    iSubjectID = subjectID{iSession};
    
        iSessionData = ccm_session_behavior(iSubjectID, iSessionID, iPlotFlag);
  populationData = [populationData; iSessionData];
end % for iSession = 1 : length(sessionArray)

save('ccm_population_human.mat', 'populationData', '-mat');



%%
iSubjectID = 'Broca';
task = 'ccm';
sessionArray = task_session_array(iSubjectID, task,'behavior2');

iPlotFlag = 0;
nSession = length(sessionArray);
populationData = [];
for iSession = 1 : nSession
    
    iSessionID = sessionArray{iSession};

            iSessionData = ccm_session_behavior(iSubjectID, iSessionID,iPlotFlag);
  populationData = [populationData; iSessionData];
end % for iSession = 1 : length(sessionArray)
save('ccm_population_broca.mat', 'populationData', '-mat');


end % function

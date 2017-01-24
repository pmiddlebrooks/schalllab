function merge_translated_files(subjectID, sessionList)


% Merges multiple files (translated from plexon to matlab), and save result
% as a new file with a new (+1) number appended.
% e.g. sessionList = {'bp063n01', bp063n02'} results in a merged file
% called 'bp063n03';

tebaDataPath    = ['/Volumes/SchallLab/data/',subjectID,'/'];
trialData       = [];
for i = 1 : length(sessionList)
   
   [trialDataX, SessionData] = load_data(subjectID, sessionList{i});
   trialDataX = struct2table(trialDataX);
   trialData = [trialData; trialDataX];
   clear trialDataX
end

trialData = table2struct(trialData, 'ToScalar', true);
trialData.SessionData = SessionData;

sessionListAll = dir([tebaDataPath,sessionList{1}(1:5),'*.plx']);
lastN = str2num(sessionListAll(end).name(7:8));
newN = num2str(1 + lastN, '%02d');
mergedSessionName = [sessionList{1}(1:6), newN];

disp(mergedSessionName)

save([tebaDataPath, mergedSessionName, '.mat'], '-struct', 'trialData','-v7.3')
saveLocalName = fullfile(local_data_path, lower(subjectID), [mergedSessionName]);
save(saveLocalName, '-struct', 'trialData','-v7.3')




% % Create local files of each individual unit data
% if isfield(SessionData, 'spikeUnitArray')
%     nUnit = length(SessionData.spikeUnitArray);
%     
%     for j = 1 : nUnit
%         jUnitName = SessionData.spikeUnitArray{j};
%         saveFileName = [mergedSessionName, '_', jUnitName];
%         
%         spikeData = trialData.spikeData(:, j);
%         save(fullfile(local_data_path, subjectID, saveFileName), 'spikeData')
%     end
% end
% 
% 
% 
% % Save a local version of behavior-only (might be duplicate,but that's
% % ok- if the session is only behavior to begin with, it won't take up
% % much memory)
% variables = fieldnames(trialData);
% removeVar = {'spikeData', 'lfpData', 'eegData'};
% 
% % If it has pyshiology data, remove it and save the file locally
% physData = ismember(variables, removeVar);
% if sum(physData)
%     trialData(:, physData) = [];
%     save([saveLocalName,'_behavior'], 'trialData', 'SessionData','-v7.3')
%     
% end

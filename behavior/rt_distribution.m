function rt = rt_distribution(subjectID, sessionID, trialList)


% Load the data
[dataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessionID);
% If the file hasn't already been copied to a local directory, do it now
if exist(localDataFile, 'file') ~= 2
    copyfile(dataFile, localDataPath)
end
load(localDataFile, 'trialData');

responseOnset   = cell2mat(trialData.responseOnset(trialList));
responseCueOn = cell2mat(trialData.responseCueOn(trialList));

rt = responseOnset - responseCueOn;

    
function localDataPath = local_data_path



homeDataPath = '/Users/paulmiddlebrooks/Dropbox/local_data';
% homeDataPath = '/Users/paulmiddlebrooks/schalllab/local_data/';
tebaDataPath = '/Volumes/SchallLab/data/';
if isdir(tebaDataPath)
    localDataPath = '/Users/paulmiddlebrooks/Dropbox/local_data/';
elseif isdir (homeDataPath)
    localDataPath = '/Users/paulmiddlebrooks/Dropbox/local_data/';
else
    disp('If you''re at work you may need to connect to teba')
    return
end



%%

plexonFile = 'File1.plx';
tebaPath = '/Volumes/SchallLab/data/Broca';
tebaFile = fullfile(tebaPath,plexonFile);
plx  = readPLXFileC(tebaFile, 'all')


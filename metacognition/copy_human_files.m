function success = copy_human_files(eyelinkFileName, matlabFileName)

edfFile = [eyelinkFileName];
% ascFile = [eyelinkFileName(1:8), '.asc'];
matFile = [matlabFileName, '.mat'];

% copy the .edf file to the usb drive
copyfile(edfFile, 'e:ChoiceStopTask/data/')
copyfile(matFile, 'e:ChoiceStopTask/data/')

% convert the edf file to an asc file
success = system(['edf2asc e:/ChoiceStopTask/data/', edfFile]);
%%

function plexon_translate_datafile_batch(monkey, partNumber)
% partNumber: e.g. for bp016n01, partNumber = 16


partNumberString = num2str(partNumber, '%03.0f')

if regexp('broca', monkey, 'ignorecase')
        tebaDataPath = ['b:/'];
elseif regexp('xena', monkey, 'ignorecase')
        tebaDataPath = ['x:/'];
else
        disp('I think you typed the monkey name wrong: see lines 10-20 in plexon_translate_datafile.m')
        return
end

dataStructure = dir([tebaDataPath,'*', partNumberString, '*.plx'])
size(dataStructure, 1)
for i = 1 : size(dataStructure, 1)
    if isempty(regexp(dataStructure(i).name, '_legacy')) &&  isempty(regexp(dataStructure(i).name, '2012'))
    dataStructure(i).name(1:end-4)
 plexon_translate_datafile(monkey, dataStructure(i).name(1:end-4))
    end
end

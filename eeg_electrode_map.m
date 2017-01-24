function [electrodeArray] = eeg_electrode_map(subjectID)

% names of electrodes, from right to left on elctrode pin in implant

switch lower(subjectID)
    case 'broca'
        electrodeArray = {'o2', 'o1', 'cz', 'fcz', 'fz'};
%         electrodeArray = {'fz', 'fcz', 'cz', 'o1', 'o2'};
    case 'xena'
        electrodeArray = {'fz', 'f1', 'f2', 'cz'};
    case 'human'
        electrodeArray = {'o2', 'o1', 'cz', 'fcz', 'fz'};
end
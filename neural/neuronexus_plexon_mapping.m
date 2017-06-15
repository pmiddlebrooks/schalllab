function [unitIndex, unitArrayNew] = neuronexus_plexon_mapping(unitArray, nChannel)
%
% Arranges spike unit array from neuronexus probe in order relative to
% plexon's recording assignments, since the 1:nChannel arrangement of the
% probe gets recorded out of order in plexon.
%
% Returns array from deepest to most shallow lead

plexonMap = 1:nChannel;

unitIndex = [];
unitArrayNew = [];

switch nChannel
    case 32
%         neuronexusMap = fliplr([9:16,25:32,17:24,1:8]);
        neuronexusMap = [9:16,25:32,17:24,1:8];
        for i = neuronexusMap
            iInd = find(cellfun(@(x) ~isempty(regexp(x, sprintf('.*%02d', i))),unitArray));
            unitIndex = [unitIndex, iInd];
            unitArrayNew = [unitArrayNew, unitArray(iInd)];
        end
        
end

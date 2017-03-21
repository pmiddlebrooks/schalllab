function [spikeUnitArrayNew, spikeDataNew] = convert_to_multiunit(spikeUnitArray, spikeData)
% function [spikeDataNew, spikeUnitArrayNew] = convert_to_multiunit(spikeData, spikeUnitArray)
%
% collapses a cell array of single unit data into a cell array of
% multi-units, combining all units on each electrode
%
% INPUT:
%
% spikeData: cell array of spike timings. One row per trial, one column per
% unit
% spikeUnitArray: cell array of the names of the units.
%
% OUTPUT:
%
% spikeDataNew: cell array of multi-units. One row per trial, one column per electrode.
%
% spikeUnitArraynew: cell array of the names of the electrodes



%%
% Testing functionality
% [td s] = load_data('joule', 'jp125n04', [ccm_min_vars,'spikeData']);

%%
% How many units on each electrode were there?
unitIndList = cellfun(@(x) str2double(regexp(x,'\d*','Match')), spikeUnitArray);

% Initialize empty cell to build name array for spike units
spikeUnitArrayNew = {};
spikeDataNew = [];

% If asking for the spike data too, process it
if nargin > 1
    % Reshape the spike data cells to concatenate them later
    spikeData = cellfun(@(x) reshape(x, [], length(x)), spikeData, 'uni', false);
    
    % Initialize a cell to fill for spikeData
    spikeDataNew = cell(size(spikeData, 1), length(unique(unitIndList)));
end

unitArray = unique(unitIndList);
for i = 1 : length(unitArray)
    iInd = unitArray(i);
    
    iUnitInd = unitIndList == iInd;
    
    spikeUnitArrayNew = [spikeUnitArrayNew, sprintf('spikeUnit%s', num2str(iInd, '%02i'))];
    if nargin > 1
        spikeDataNew(:,i) = arrayfun(@(IDX) [spikeData{IDX,iUnitInd}], 1:size(spikeData,1), 'Uniform', 0);
    end
end

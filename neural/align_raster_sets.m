function [alignedRasters, alignmentIndex] = align_raster_sets(rasterSet, alignmentIndexList)

% Input a cell array of sets of aligned rasters, each set with an alignment
% index, and align the set of rasters, computing the new alignment index

if ~iscell(alignmentIndexList)
    alignmentIndexList = num2cell(alignmentIndexList);
end

nRasterSet = length(rasterSet);

% Determine which raster set has the most time from beginning to alignemnt
% index, and use that value to pad other raster sets' begninnings
maxBegin = max(cell2mat(alignmentIndexList));
alignmentIndex = maxBegin;


ends = cell2mat(cellfun(@(x, y) size(x, 2) - y, rasterSet, alignmentIndexList, 'uniformoutput', false));
maxEnd = max(ends);

alignedRasters = [];
for i = 1 : nRasterSet
    iRaster = rasterSet{i};
    nTrial = size(iRaster, 1);
    
    iBeginPad = nan(nTrial, maxBegin - alignmentIndexList{i});
    iEndPad = nan(nTrial, maxEnd - (size(iRaster, 2) - alignmentIndexList{i}));
    
    iRaster = [iBeginPad, iRaster, iEndPad];
    
    alignedRasters = [alignedRasters; iRaster];
end
    
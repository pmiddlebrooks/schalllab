function [alignedRasters, alignmentIndex] = spike_to_raster(spikeTrainMatrix, alignmentTimeList, debugMode)

% [alignedRasters, alignmentIndex] = spike_to_raster(spikeTrainMatrix, alignmentTimeList)
%
% Converts a spike train matrix or vector into rasters. If a
% single trial is input, a single raster is returned. If a matrix or cell array of trials
% is input, an matrix of rasters is returned aligned on alignmentIndex
%
% spikeTrainMatrix: A vector/matrix of trial/trials containing spike times.
% Each row is a trial, each column is one ms
%
% alignmentTimeList (optional): A vector list of times to align each respective trial to form an aligned matrix of rows of rasters
%
%
%
% alignedRasters: The raster or aligned rasters
%
% alignmentIndex: The index that the rasters are aligned to if a spike
% train matrix was input.
%
%
% Two main parts:
% 1. Change the spikeTrain matrix to rasters
% 2. Align the rasters on alignmentTimeList
%
% spikeTrainMatrix:         [trials X spike trains], each nonzero number in
% a row is the time of a spike on that trial
%
% alignmentTimeList:        [trials X 1], each number is the time to align
% the rasters for each trial.
%
% debugMode: 0 (default) to supress warning messages (like, "no spikes in spikeTrainMatrix, etc), 1 to print them to screen

if nargin < 3
    debugMode = 0;
end

if isempty(spikeTrainMatrix)
%     disp('In spike_to_raster.m, apparently there are no spikes in the spike train')
    alignedRasters = [];
    alignmentIndex = [];
    return
end

nTrial = [];

% Determine how many trials we're dealing with
% and
% If it's a cell, transform it into a matrix with NaNs padding
if iscell(spikeTrainMatrix)
nTrial = length(spikeTrainMatrix);
    % reshape cell contents if necessary
    if sum(cellfun(@(x) size(x, 1)>1, spikeTrainMatrix))
        spikeTrainMatrix = cellfun(@(x) x', spikeTrainMatrix, 'uniformoutput', false);
    end
    maxSpikes = round(max(cellfun(@length, spikeTrainMatrix)));
    spikeTrainCell = cellfun(@(x) [x, nan(1, maxSpikes-length(x))], spikeTrainMatrix, 'uniformoutput', false);
    spikeTrainMatrix = cell2mat(spikeTrainCell);
%     if size(spikeTrainMatrix, 2) == 1
%     alignedRasters = [];
%     alignmentIndex = [];
%     return
%     end 
else
% If it's a vector, there's only one trial.
if min(size(spikeTrainMatrix)) == 1 && ~sum(isnan(spikeTrainMatrix))
        nTrial = 1;
        if size(spikeTrainMatrix, 2) == 1
            spikeTrainMatrix = spikeTrainMatrix';
        end
else
    nTrial = size(spikeTrainMatrix, 1);
end
end
    
%         % If it's a matrix, trials might be included as rows or columns.
%         % Figure that out and 
% elseif ~sum(isnan(spikeTrainMatrix(:, 1)))
%         nTrial = size(spikeTrainMatrix(:,1));
% elseif ~sum(isnan(spikeTrainMatrix(1, :)))
%         nTrial = size(spikeTrainMatrix(1,:));
% end



if nargin < 2
   alignmentTimeList = ones(nTrial, 1);
end






% Return empty matricies if alignTimeList is all NaNs (e.g. if the few
% trials we wanted to analyze were discarded for not making criteria as
% valid saccades to target, etc)
if nargin > 1 && sum(~isnan(alignmentTimeList)) == 0
    alignedRasters = [];
    alignmentIndex = [];
    return
end










% Make sure there are no spikes at negative times- sesms that some are at
% -1, -2, etc
spikeTrainMatrix(spikeTrainMatrix < 0) = 1;
spikeTrainMatrix(spikeTrainMatrix > 1e5) = nan;

% Figure out the latest time a spike occurs among all the trials. This will
% be used to help initialize the size of rasters matrices
lastSpikeTime = round(max(spikeTrainMatrix(:)));
% If there were no spikes, let user know, and just make a bit raster of
% zeros big enough, with align times in the middle. This is for the rare
% case when a neuron did not fire at all over many trials (in a certain
% condition, or for a very low firing neuron, etc).
if isempty(lastSpikeTime)
    if debugMode
    disp('In spike_to_raster.m, apparently there are no spikes in the spike train')
    end
    alignedRasters = int8(zeros(nTrial, 6000));
    alignmentIndex = 3000;
    return
end


% if size(spikeTrainMatrix, 2) == 1 && size(spikeTrainMatrix, 1) > 1 && nTrial == 1
%     spikeTrainMatrix = spikeTrainMatrix';
% end

    




% % EXTRA_PAD ms will be added to the end of alignedRasters once
% % alignedRasters is assigned. This is done because we want the full trial-
% % otherwise the raster would extend only to the last spike.
% if nargin > 2 && ~isempty(trialDuration)
%     EXTRA_PAD = zeros(nTrial, trialDuration - lastSpikeTime);
% end
% 



% % Make sure the user inputs a trialDuration in the case that the user input
% % a single trial of spike data
% if nTrial == 1 && (nargin < 3 || isempty(trialDuration))
%     disp('You need to enter a trialDuration for a single trial so the raster continues beyond the last spike time')
%     return
% end






% ////////////////////////////////////////////////////////////////////////
% 1. Change the spikeTrain matrix to rasters


% Initialize the alignedRasters all to zeros, nTrial trials long
unAlignedRasters = zeros(nTrial, lastSpikeTime);
% Loop through and add rasters to each trial index when spikes occured, and
% add nans on each trial after the last spike on that trial until the last
% spike from all the trials
for iTrial = 1 : nTrial
    %     iLastSpikeTime = max(spikeTrainMatrix(iTrial, :));
    %     if iLastSpikeTime > 0
    %         unAlignedRasters(iTrial, :) = zeros(1, lastSpikeTime);
    %     end
    %     unAlignedRasters(iTrial, iLastSpikeTime+1 : lastSpikeTime) = nan;
    spikeIndices = round(nonzeros(spikeTrainMatrix(iTrial, :)));
    unAlignedRasters(iTrial, spikeIndices(spikeIndices > 0)) = 1;
end


% size(unAlignedRasters)



% ////////////////////////////////////////////////////////////////////////
% 2. Align the rasters on alignmentTimeList
% function [aligned] = align_rows_on_indices(rows, alignedindexlist)

% If user input a single raster trial or the rasters are already
% aligned, we're done here. Else, we need to align stuff
if isempty(alignmentTimeList) || length(alignmentTimeList) == 1
    alignedRasters = unAlignedRasters;
else
    % Make sure the trial number in the spikeTrainMatrix and the alignmentTime list match:
    if nTrial ~= length(alignmentTimeList)
        fprintf( 'In spike_to_raster.m, the number of trial between spikeTrainMatrix (%d) and alignmentTimeList (%d) does not match\n', ...
            nTrial, length(alignmentTimeList));
        return;
    end
    
    minimumIndex = min(alignmentTimeList);
    maximumIndex = max(alignmentTimeList);
    
    % if maximumIndex > lastSpikeTime
    %     disp( 'In spike_to_raster.m, at least one of the alignmentTimeList indices is later than the latest spike time in spikeTrainMatrix');
    %     return;
    % end
    
    % maximumShift is the amount we have to increase the matrix size.
    maximumShift = maximumIndex - minimumIndex;
    if maximumShift == 0 && nTrial > 1
        if debugMode
        disp( 'In spike_to_raster.m, check your rsater to make sure: seems the rasters were already aligned when input' );
        end
    end
    
    % Initializes the alignedRasters to be as long as the required shift,
    % considering the aligmentTimeList, and add the extra padding (see
    % EXTRA_PAD comments above)
    alignedRasters = zeros(nTrial, (lastSpikeTime + maximumShift));
    % Loop through each trial
    for iTrial = 1 : nTrial
        iShift = maximumIndex - alignmentTimeList(iTrial);
        iFirstIndex = iShift + 1;
        iLastIndex = iFirstIndex + lastSpikeTime - 1;
        if isnan(iFirstIndex)
            iFirstIndex = 1;
            iLastIndex = lastSpikeTime;
        end
        alignedRasters(iTrial, iFirstIndex : iLastIndex) = unAlignedRasters(iTrial, :);
        alignedRasters(iTrial, iLastIndex+1 : end) = nan;
        
        
    end
end
% sdf = spike_density_function(alignedRasters, 'gaussian', 10);
% plot(nanmean(sdf, 1))
if nargin > 1 && ~isempty(alignmentTimeList)
    alignmentIndex = max(alignmentTimeList);
else
    alignmentIndex = [];
end

% For low-firing neurons, there may be so few spikes (or none at all) that
% the returned rasters will be shorter than the trial actually was. Add
% some padding here if that seems to be the case
catchLimit = 600;  % ms
if size(alignedRasters, 2) - alignmentIndex < catchLimit
    addNaN = catchLimit - (size(alignedRasters, 2) - alignmentIndex);
    alignedRasters = [alignedRasters, nan(nTrial, addNaN)];
end



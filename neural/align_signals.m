function [alignedSignals, alignmentIndex] = align_signals(signalMatrix, alignmentTimeList, window)

% [alignedSignals, alignmentIndex] = align_signals(s, alignmentTimeList, trialDuration)
%
% This function aligns a signal matrix. If a
% single trial is input, a single signal vector is returned. If a matrix of trials
% is input, an matrix of rasters is returned aligned on an alignmentIndex
%
% signalMatrix: A cell array of signals containing raw signal values
% Each element is a trial, sampled at 1 kHz
%
% alignmentTimeList (optional): A vector list of times to align each respective trial to form an aligned matrix of rows of rasters
%
%
% alignedSignals: The aligned signals, in a double matrix
%
% alignmentIndex: The index that the signals are aligned to.
%
% window: e.g. (-99 : 400)  A range of time around which to cut of the aligned signal: may be
% useful to keep the data smaller (uses less memory)



% signalMatrix:         [trials X spike trains], each nonzero number in
% a row is the time of a spike on that trial

% alignmentTimeList:        [trials X 1], each number is the time to align
% the rasters for each trial

alignedSignals = [];
alignmentIndex = [];
if isempty(signalMatrix)
   %     disp('In spike_to_raster.m, apparently there are no spikes in the spike train')
   return
end
% Return empty matricies if alignTimeList is all NaNs (e.g. if the few
% trials we wanted to analyze were discarded for not making criteria as
% valid saccades to target, etc)
if nargin > 1 && sum(~isnan(alignmentTimeList)) == 0
   return
end
if nargin < 3
   window = [];
end


if ~iscell(signalMatrix)
   disp('Need the signal matrix to be input as a cell array')
   return
end
% reshape cell contents if necessary
if sum(cellfun(@(x) size(x, 1)>1, signalMatrix))
   signalMatrix = cellfun(@(x) x', signalMatrix, 'uniformoutput', false);
end

maxTrialLength = max(cellfun(@length, signalMatrix));
nTrial = length(signalMatrix);









% ////////////////////////////////////////////////////////////////////////
% 2. Align the signals on alignmentTimeList
% ////////////////////////////////////////////////////////////////////////

% If user input a single raster trial or the rasters are already
% aligned, we're done here. Else, we need to align stuff
if isempty(alignmentTimeList) || length(alignmentTimeList) == 1
   alignedSignals = cell2mat(signalMatrix);
else
   % Make sure the trial number in the signalMatrix and the alignmentTime list match:
   if nTrial ~= length(alignmentTimeList)
      fprintf( 'In align_signals.m, the number of trial between signalMatrix (%d) and alignmentTimeList (%d) does not match\n', ...
         nTrial, length(alignmentTimeList));
      return;
   end
   
   minimumIndex = min(alignmentTimeList);
   maximumIndex = max(alignmentTimeList);
   
   % maximumShift is the amount we have to increase the matrix size.
   maximumShift = maximumIndex - minimumIndex;
   if maximumShift == 0 && nTrial > 1
      disp( 'In align_signals.m, check your signal matrix to make sure: seems the signals were already aligned when input' );
   end
   
   % Initializes the alignedSignals to be as long as the required shift,
   % considering the aligmentTimeList, and add the extra padding (see
   % EXTRA_PAD comments above)
   alignedSignals = nan(nTrial, (maxTrialLength + maximumShift));
   % Loop through each trial
   for iTrial = 1 : nTrial
      trialDuration = length(signalMatrix{iTrial});
      iShift = maximumIndex - alignmentTimeList(iTrial);
      iFirstIndex = iShift + 1;
      iLastIndex = iFirstIndex + trialDuration - 1;
      if isnan(iFirstIndex)
         iFirstIndex = 1;
         iLastIndex = trialDuration;
      end
      %         iLastIndex - iFirstIndex
      %         length(signalMatrix{iTrial})
      alignedSignals(iTrial, iFirstIndex : iLastIndex) = signalMatrix{iTrial};
      
      
      
   end
end

% alignmentIndex is the maximum of alignTimeList
if nargin > 1 && ~isempty(alignmentTimeList)
   alignmentIndex = max(alignmentTimeList);
else
   alignmentIndex = [];
end



% Finally, crop a window of time of the signal to return if desired.
if ~isempty(window)
   if window(1) <= -alignmentIndex
      alignedSignals = [nan(nTrial, alignmentIndex + window(1)), alignedSignals];
      alignmentIndex = -window(1) + 1;
   end
   if window(end) > size(alignedSignals, 2) - alignmentIndex
      alignedSignals = [alignedSignals, nan(nTrial, window(end))];
   end
   alignedSignals = alignedSignals(:, window + alignmentIndex);
   if nargin > 1 && ~isempty(alignmentIndex)
      alignmentIndex = -window(1) + 1;
   end
end



function signal = signal_baseline_correct(signal, baseWindow, alignIndex)


% INPUT
%
% signal:          a matrix of trials (rows) of signal (eeg,lfp) data (time in columns) aligned.
%
% baseWindow: window of time relative to alignment index to zero (thus
%   shift whole signal by that amount).
%
% alignIndex: Index with the signal upon which the data is aligned.

if isempty(signal)
return
end

baseShift = nanmean(signal(:,alignIndex + baseWindow), 2);
baseShift = repmat(baseShift, 1, size(signal, 2));


signal = signal - baseShift;
function [satTrial] = signal_reject_saturate(signal, varargin)


% INPUT
%
% signal:          a matrix of trials (rows) of eeg  or lfp data (time in columns) aligned.
% vargargin:    
%   1. rejectWindow: can narrow the windown within which to probe for
%   rejection criterium. E.g. might want to reject if the signal satureates
%   within the 100 : 200 ms window of a signal 1000 ms long
%   2. satThreshold: How many ms will we allow a signal to be saturated
%   without rejecting it?
rejectWindow    = [];  % By default, treat whole sample as a window to reject
satThreshold    = 50; % How many ms will we allow a signal to be saturated?
alignIndex      = []; % Where does the rejectWindow sit in the signal?

for i = 1 : 2 : length(varargin)
   switch varargin{i}
      case 'rejectWindow'
         rejectWindow = varargin{i+1};
      case 'satThreshold'
         satThreshold = varargin{i+1};
      case 'alignIndex'
         alignIndex = varargin{i+1};
      otherwise
   end
end



nTrial      = size(signal, 1);
satTrial    = zeros(nTrial, 1);

for i = 1 : nTrial
      % Reject trials with saturated signals
      if sum(diff(signal(i, alignIndex + rejectWindow)) == 0) > satThreshold
         satTrial(iTrial) = 1;
      end
end
satTrial = find(satTrial);

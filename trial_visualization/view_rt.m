
%%


subjectID = 'xena';
sessionID = 'xp036n02';

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
trialData = struct2dataset(trialData);
trialData(isnan(trialData.rt), :) = [];

nTrial = size(trialData.trialOutcome, 1);


figure(33);
clf
for i = 1 : nTrial
   
   cla
   hold on
   
%    plot(trialData.eyeX{i}, 'r');
%    plot(trialData.eyeY{i}, 'b');
%    plot(trialData.responseOnset(i), trialData.eyeX{i}(trialData.responseOnset(i)), 'ok')
%    plot(trialData.responseOnset(i), trialData.eyeY{i}(trialData.responseOnset(i)), 'ok')
%    xlim([trialData.responseOnset(i) - 70 trialData.responseOnset(i) + 70])
   
   plot(trialData.eyeX{i}, trialData.eyeY{i});
   plot(trialData.eyeX{i}(trialData.saccBegin{i}(trialData.saccToTargIndex(i))), trialData.eyeY{i}(trialData.saccBegin{i}(trialData.saccToTargIndex(i))), 'og');
   plot(trialData.eyeX{i}(trialData.saccBegin{i}(trialData.saccToTargIndex(i)) + trialData.saccDuration{i}(trialData.saccToTargIndex(i))), trialData.eyeY{i}(trialData.saccBegin{i}(trialData.saccToTargIndex(i)) + trialData.saccDuration{i}(trialData.saccToTargIndex(i))), 'or');
   
   % % EEG data
   %    subplot(1,2,1)
   %    cla
   %    plot(trialData.eegData{i,1})
   %    subplot(1,2,2)
   %    cla
   %    plot(trialData.eegData{i,3})
   
   disp(trialData(i,:))
   
   pause
end
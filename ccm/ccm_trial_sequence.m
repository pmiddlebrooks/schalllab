function [trialList] = ccm_trial_sequence(trialData, seqOpt)

% %% TESTING
% 
% [trialData, S, E] = load_data('broca','bp174n02');
% 
% %%
% % Get default options structure:
% selectOpt = ccm_trial_selection;
% 
% seqOpt = selectOpt;
% seqOpt(1).outcome = {'goCorrectTarget'};
% seqOpt(1).ssd = 'none';
% 
% seqOpt(2) = selectOpt;
% seqOpt(2).outcome = {'stopStop'};
% 
% seqOpt(3) = selectOpt;
% seqOpt(3).outcome = {'goCorrectTarget'};
% seqOpt(3).ssd = 'none';


% % Required selectOpt structure fields:
%    selectOpt.outcome            = 'valid';
%    selectOpt.rightCheckerPct   = 'collapse';
%    selectOpt.ssd               = 'collapse';
%    selectOpt.targDir           = 'collapse';
%    selectOpt.responseDir       = 'collapse';

nSequence = length(seqOpt);
nTrial = size(trialData.trialOutcome, 1);
trialList = 1:nTrial;


for i = 1 : nSequence
    
    iTrialList = ccm_trial_selection(trialData, seqOpt(i));
    
    trialList = intersect(trialList, iTrialList-i+1);
end
function Data = ccm_interspike_interval(subject, session, unit, Opt)
%
% function Data = ccm_session_data(subject, sessionID, plotFlag, Opt.unitArray)
%
% Creates a processed data struct for other functions to use, and
% optionally plots session neurophysiology data
%
%
% Single neuron analyses for choice countermanding task. Only plots the
% sdfs. To see rasters, use ccm_single_neuron_rasters, which displays all
% conditions in a given epoch
%
% input:
%   subject: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%
%   Opt: A structure with various ways to select/organize data:
%   The following fields of Opt are relevant for ccm_session_data, with
%   possible values (default listed first):
%
%    Opt.trialData: if Options structure contains trialData and
%    SessionData, don't need to load the data
%   Opt.SessionData: ditto
%    Opt.dataType = 'neuron', 'lfp', 'erp';
%
%    Opt.figureHandle   = 1000;
%    Opt.printPlot      = false, true;
%    Opt.plotFlag       = true, false;
%    Opt.hemisphere = which hemsiphere were the data recorded from? left or
%    right?
%
%
%
% Returns Data structure with fields:
%
%   Data.(epoch).colorCoherence(x).(outcome).ssd(x)
%
%   outcome can be:  goTarg, goDist, stopTarg, stopDist, stopStop
%   ssd(x):  only applies for stop trials, else the field is absent
%   epoch name: fixOn, targOn, checkerOn, etc.

%%
if nargin < 4
    Opt = ccm_options;
    Opt.trialData = [];
end

% User may have input the data already, so you don't need to load it
if nargin < 4
    Opt = ccm_options;
    Opt.multiUnit = false;
    Opt.printFlag = 0;
    Opt.plotFlag = 0;
    Opt.rf = [];
    Opt.hemisphere = [];
end

% Get RF if it wasn't input:
rfOpt                 = ccm_options;
rfOpt.collapseSignal  = true;
rfOpt.collapseTarget  = true;
rfOpt.doStops         = false;
rfOpt.plotFlag        = false;
rfOpt.printPlot       = false;
rfOpt.multiUnit       = Opt.multiUnit;
    iUnit               = [{session}, {unit}];
    rfData               = ccm_session_data(subject, iUnit, rfOpt);
Opt.rf = ccm_find_saccade_rf(rfData);



% CONSTANTS
MIN_RT          = 120;
MAX_RT          = 1200;
STD_MULTIPLE    = 3;
saccadeEpoch = 100; % How many ms previous from saccade onset do we start the analysis?

% Load and process data
% ============================================================
variables = [ccm_min_vars, 'spikeData'];
[trialData, SessionData, ExtraVar] = load_data(subject, session, variables, Opt.multiUnit);
pSignalArray    = ExtraVar.pSignalArray;


if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a chioce countermanding saccade session, try again\n')
    return
end


% Get rid of trials with outlying RTs
[allRT, rtOutlierTrial] = truncate_rt(trialData.rt, MIN_RT, MAX_RT, STD_MULTIPLE);
trialData.rt(rtOutlierTrial) = nan;

trialData = ccm_delete_nan_rt(trialData);


% Which spike data index do we need to process?
spikeDataInd = ismember(SessionData.spikeUnitArray, unit);
% Make sure user input a dataType that was recorded during the session
if ~sum(spikeDataInd)
    fprintf('Session %s apparently does not contain %s data \n', session, unit)
    return
end






% Determine which color coherence to use:
% ============================================================
switch Opt.rf
    case 'left'
        coherence = pSignalArray(1);
    case 'right'
        coherence = pSignalArray(end);
    otherwise
        switch Opt.hemisphere
            case 'left'
                coherence = pSignalArray(end);
            case 'right'
                coherence = pSignalArray(1);
        end
end







% Get trials to analyze
% ============================================================
% Get default trial selection Opt
selectOpt = ccm_options;
selectOpt.rightCheckerPct = coherence * 100;
selectOpt.rightCheckerPct = coherence * 100;
selectOpt.outcome     = {'goCorrectTarget'};

trial = ccm_trial_selection(trialData, selectOpt);





% Process checker onset data
% ============================================================
checkerWinBegin = num2cell(trialData.checkerOn(trial));
checkerWinEnd = num2cell(trialData.checkerOn(trial) + trialData.rt(trial));
Data.checkerSpikeInterval = cellfun(@(in1,in2,in3) diff(sort(in1(in1 > in2 & in1 < in3))), trialData.spikeData(trial,spikeDataInd), checkerWinBegin, checkerWinEnd, 'uni', false);





% Process saccade onset data
% ============================================================
saccadeWinBegin = num2cell(trialData.responseOnset(trial) - saccadeEpoch);
saccadeWinEnd = num2cell(trialData.responseOnset(trial));
spikeData = cellfun(@(in1, in2) [in1 in2], trialData.spikeData(trial,spikeDataInd), num2cell(trialData.responseOnset(trial)), 'uni', false); 
% Data.saccadeSpikeInterval = cellfun(@(in1,in2,in3) diff(sort(in1(in1 > in2 & in1 < in3))), trialData.spikeData(trial,spikeDataInd), saccadeWinBegin, saccadeWinEnd, 'uni', false);
Data.saccadeSpikeInterval = cellfun(@(in1,in2,in3) diff(sort(in1(in1 > in2 & in1 < in3))), spikeData, saccadeWinBegin, saccadeWinEnd, 'uni', false);

figure()
hold all
cellfun(@(x) plot(1:length(x), flipud(x(:))), Data.saccadeSpikeInterval)
return



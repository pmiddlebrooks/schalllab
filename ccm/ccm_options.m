function Opt = ccm_options
% function Opt = ccm_options
%
% Returns default Opt structure for use as input in choice
% countermanding analyses functions. A structure with various ways to
% select/organize data:
%
% Opt has the following fields with
%   possible values (default listed first):
%
%    Opt.dataType = 'neuron', 'lfp', 'erp';
%
%    Opt.figureHandle   = 1000;
%    Opt.printPlot      = false, true;
%    Opt.plotFlag       = true, false;
%    Opt.collapseSignal = false, true;
%    Opt.collapseTarg         = false, true;
%    Opt.doStops        = true, false;
%    Opt.filterData 	= false, true;
%    Opt.stopHz         = 50, <any number, above which signal is filtered;
%    Opt.normalize      = false, true;
%    Opt.howProcess      = how to step through the list of units
%                                 'each' to plot all,
%                                 'step' (default): step through to see one
%                                 plot at a time, pausing between
%                                 'print' (default): step through each plot
%                                 individually, printing each to file
%    Opt.unitArray      = units you want to analyze (default gets filled
%                                   with all available).
%                                 {'spikeUnit17a'}, input a cell of units for a list of individaul units
%   Opt.baselineCorrect = false, true; Baseline correct analog signals?
%
%     Opt.outcome  = array of strings indicating the outcomes to include:
%           {'collapse',
%           'valid'
%           'goCorrectTarget', 'goCorrectDistractor',
%           'goIncorrect',
%           'stopCorrect',
%           'stopIncorrectTarget', 'stopIncorrectDistractor',
%           'targetHoldAbort', 'distractorHoldAbort',
%           'fixationAbort', 'saccadeAbort', 'checkerStimulusAbort'}
%
%           valid = any non-aborts
%
%     Opt.choiceAccuracy  = default is collapse across all choices. Segments
%           with respect to weather a choice was correct or error (does not count
%           if no choice was made.
%           Opt: 'collapse' (default), 'correct', 'error'.
%
%     Opt.rightCheckerPct  = range of checkerboard percentage of right target checkers:
%           {'collapse', 'right', 'left'
%           a double array containing the values, e.g. [40 50 60]
%     Opt.ssd    = range of SSDs to include in the trial list:
%           {'collapse', 'any', 'none', or
%           a double array containing the values, e.g. [43 86 129]
%     Opt.allowRtPreSsd    = whether to allow noncanceled stop trials with RTs before SSDs
%           true (default) or false
%     Opt.targDir  = the angle of the CORRECT TARGET
%           {'collapse','left','right','up','down','upLeft','upRight',downLeft','downRight'<list of possible angels, e.g. [45 -45]
%     Opt.responseDir  = the angle of target to which a response was made
%           {'collapse','none' (for stopStop trials),'left','right','up','down','upLeft','upRight',downLeft','downRight'<list of possible angels, e.g. [45 -45]

Opt.dataType = 'neuron';

Opt.trialData           = [];
Opt.howProcess        = 'each';%'print';
Opt.unitArray        = [];

Opt.figureHandle     = 1000;
Opt.printPlot        = true;
Opt.plotFlag         = true;

Opt.collapseSignal   = false;
Opt.collapseTarg     = true;
Opt.doStops          = true;

Opt.filterData       = false;
Opt.stopHz           = 50;
Opt.normalize        = false;
Opt.baselineCorrect  = true;
Opt.epochName       = 'checkerOn';
Opt.eventMarkName       = 'responseOnset';
Opt.epochWindow     = -299:300;


Opt.outcome             = 'valid';
Opt.choiceAccuracy      = 'collapse';
Opt.rightCheckerPct     = 'collapse';
Opt.ssd                 = 'any';
Opt.allowRtPreSsd     = true;
Opt.targDir             = 'collapse';
Opt.responseDir    = {'collapse'}; % {'left', 'right'};
Opt.ssrt                 = [];


Opt.include50           = false; % include 50% color coherence conditions?
Opt.deleteAborts     = true;
Opt.USE_PRE_SSD   	= true; % Use noncanceled stop trials when RT precedes SSD?
Opt.USE_TWO_COLORS         = false; % Only use one easy and one hard color for both sides

Opt.latencyMatchMethod 	= 'ssrt';
Opt.minTrialPerCond     = 15;
Opt.cellType            =      'presacc';

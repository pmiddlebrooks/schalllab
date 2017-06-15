function Opt = mem_options
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
%           'goCorrectTarget', 'targetHoldAbort',
%           'fixationAbort', 'saccadeAbort'}
%
%           valid = any non-aborts
%
%     Opt.targDir  = the angle of the CORRECT TARGET
%           {'collapse','left','right','up','down','upLeft','upRight',downLeft','downRight'<list of possible angels, e.g. [45 -45]
%     Opt.responseDir  = the angle of target to which a response was made
%           {'collapse','none' (for stopStop trials),'left','right','up','down','upLeft','upRight',downLeft','downRight'<list of possible angels, e.g. [45 -45]

Opt.task                = 'mem';  % CAn use 'del' to analyzed delayed saccade task with mem codes.
Opt.dataType            = 'neuron';
Opt.multiUnit           = false;
Opt.trialData           = [];
Opt.howProcess        = 'each';%'print';
Opt.unitArray        = [];

Opt.figureHandle     = 1000;
Opt.printPlot        = true;
Opt.plotFlag         = true;

Opt.collapseSignal   = false;
Opt.collapseTarg     = true;

Opt.filterData       = false;
Opt.stopHz           = 50;
Opt.normalize        = false;
Opt.baselineCorrect  = true;
Opt.epochName       = 'targOn';
Opt.eventMarkName       = 'responseOnset';
Opt.epochWindow     = -299:300;


Opt.outcome             = 'goCorrectTarget';
Opt.targDir             = 'collapse';
Opt.responseDir    = {'collapse'}; % {'left', 'right'};


Opt.deleteAborts     = true;

Opt.cellType            =      'presacc';

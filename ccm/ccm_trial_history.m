function data = ccm_trial_history(subjectID, sessionID, plotFlag, figureHandle)

if nargin < 3
   plotFlag = 1;
end
if nargin < 4
   figureHandle = 4445;
end

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
ssdArray = ExtraVar.ssdArray;
pSignalArray = ExtraVar.pSignalArray;

% Constants
DO_STOPS = 1;
MIN_RT = 120;
MAX_RT = 1200;
nSTD   = 3;


if ~strcmp(SessionData.taskID, 'ccm')
   fprintf('Not a choice countermanding session, try again\n')
   return
end






if plotFlag
   % axes names
   axXGoTarg = 1;
   axXGoTarg2 = 2;
   axXGo = 3;
   axXGo2 = 4;
   axGoDistX = 5;
   axGoDistX2 = 6;
   axStopP = 7;
   axTriplet = 8;
   
   nRow = 3;
   nColumn = 3;
   screenOrSave = 'screen';
   [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nRow, nColumn, 'portrait', figureHandle);
   clf
   choicePlotXMargin = .03;
   ssdMargin = 20;
   ylimArray = [];
   
   
   ax(axXGoTarg) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
   hold(ax(axXGoTarg), 'on')
   ax(axXGoTarg2) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
   hold(ax(axXGoTarg2), 'on')
   ax(axXGo) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 1) yAxesPosition(2, 1) axisWidth axisHeight]);
   hold(ax(axXGo), 'on')
   ax(axXGo2) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2) yAxesPosition(2, 2) axisWidth axisHeight]);
   hold(ax(axXGo2), 'on')
   ax(axGoDistX) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 1) yAxesPosition(3, 1) axisWidth axisHeight]);
   hold(ax(axGoDistX), 'on')
   ax(axGoDistX2) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 2) yAxesPosition(3, 2) axisWidth axisHeight]);
   hold(ax(axGoDistX2), 'on')
   ax(axStopP) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 3) yAxesPosition(2, 3) axisWidth axisHeight]);
   hold(ax(axStopP), 'on')
   ax(axTriplet) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 3) yAxesPosition(1, 3) axisWidth axisHeight]);
   hold(ax(axTriplet), 'on')
end




% Get rid of trials without valid conclusions: i.e.
deleteTrial = strcmp(trialData.trialOutcome, 'noFixation') | ...
   strcmp(trialData.trialOutcome, 'fixationAbort') | ...
   strcmp(trialData.trialOutcome, 'choiceStimulusAbort') | ...
   strcmp(trialData.trialOutcome, 'saccadeAbort');

trialData(deleteTrial, :) = [];


% Truncate RTs
[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTD);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);





nTrial = size(trialData.trialOutcome, 1);










% ********************************************************************
% xxxxxxxxxxxxxxxxx
% ********************************************************************

% Go Target RT as a function of previous trial outcome
% ********************************************************************
goTargGoTargTrial       = cell(1, length(pSignalArray));
goDistGoTargTrial       = cell(1, length(pSignalArray));
goGoTargTrial           = cell(1, length(pSignalArray));
stopTargGoTargTrial     = cell(1, length(pSignalArray));
stopDistGoTargTrial     = cell(1, length(pSignalArray));
stopIncorrectGoTargTrial  = cell(1, length(pSignalArray));
stopStopGoTargTrial	= cell(1, length(pSignalArray));

goTargGoTargGoTargTrial       = cell(1, length(pSignalArray));
goTargGoDistGoTargTrial	= cell(1, length(pSignalArray));
goTargStopIncorrectGoTargTrial  = cell(1, length(pSignalArray));
goTargStopCorrectGoTargTrial	= cell(1, length(pSignalArray));

goTargGoTargRT       = cell(1, length(pSignalArray));  % rt to the target after a go target trial
goTargRTGoTarg       = cell(1, length(pSignalArray));  % rt to the target before a go target trial
goDistGoTargRT       = cell(1, length(pSignalArray));  % rt to the target after a go distractor trial
goDistRTGoTarg       = cell(1, length(pSignalArray));  % rt to the distractor before a go target trial
goGoTargRT           = cell(1, length(pSignalArray));
stopTargGoTargRT     = cell(1, length(pSignalArray));
stopDistGoTargRT     = cell(1, length(pSignalArray));
stopIncorrectGoTargRT  = cell(1, length(pSignalArray));
stopStopGoTargRT   = cell(1, length(pSignalArray));

goTargGoTargGoTargRT       = cell(1, length(pSignalArray));  % rt to the target after a go target trial
goTargRTGoTargGoTarg       = cell(1, length(pSignalArray));  % rt to the target before a go target trial
goTargStopIncorrectGoTargRT  = cell(1, length(pSignalArray)); % rt to last goTarget in trial sequence goTarget --> stopIncorrect --> goTarget
goTargStopCorrectGoTargRT   = cell(1, length(pSignalArray)); % rt to last goTarget in trial sequence goTarget --> stopStop --> goTarget
goTargGoDistGoTargRT   = cell(1, length(pSignalArray)); % rt to last goTarget in trial sequence goTarget --> stopStop --> goTarget
goTargRTStopIncorrectGoTarg  = cell(1, length(pSignalArray));  % rt to first goTarget in trial sequence goTarget --> stopInorrect --> goTarget
goTargRTStopCorrectGoTarg   = cell(1, length(pSignalArray)); % rt to first goTarget in trial sequence goTarget --> stopStop --> goTarget
goTargRTGoDistGoTarg   = cell(1, length(pSignalArray)); % rt to first goTarget in trial sequence goTarget --> stopStop --> goTarget


% Go (to target or distractor) RT as a function of previous trial outcome
% ********************************************************************
goTargGoTrial       = cell(1, length(pSignalArray));
goDistGoTrial       = cell(1, length(pSignalArray));
goGoTrial           = cell(1, length(pSignalArray));
stopTargGoTrial     = cell(1, length(pSignalArray));
stopDistGoTrial     = cell(1, length(pSignalArray));
stopIncorrectGoTrial  = cell(1, length(pSignalArray));
stopStopGoTrial	= cell(1, length(pSignalArray));

goTargGoRT       = cell(1, length(pSignalArray));
goDistGoRT       = cell(1, length(pSignalArray));
goGoRT           = cell(1, length(pSignalArray));
stopTargGoRT     = cell(1, length(pSignalArray));
stopDistGoRT     = cell(1, length(pSignalArray));
stopIncorrectGoRT  = cell(1, length(pSignalArray));
stopStopGoRT   = cell(1, length(pSignalArray));



% Trial outcomes as a function of previous Go to distractor trial
% ********************************************************************
% goDistGoTargTrial       = cell(1, length(pSignalArray));
goDistGoDistTrial       = cell(1, length(pSignalArray));
% goDistGoTrial           = cell(1, length(pSignalArray));
goDistStopTargTrial     = cell(1, length(pSignalArray));
goDistStopDistTrial     = cell(1, length(pSignalArray));
goDistStopIncorrectTrial  = cell(1, length(pSignalArray));
% goDistStopCorrectTrial	= cell(1, length(pSignalArray));

% goDistGoTargRT       = cell(1, length(pSignalArray));
goDistGoDistRT       = cell(1, length(pSignalArray));
% goDistGoRT           = cell(1, length(pSignalArray));
goDistStopTargRT     = cell(1, length(pSignalArray));
goDistStopDistRT     = cell(1, length(pSignalArray));
goDistStopIncorrectRT  = cell(1, length(pSignalArray));
% goDistStopCorrectRT   = cell(1, length(pSignalArray));


% Stop probability in this signal strength as a function of previous
% trial outcome in any signal strength
% ********************************************************************
pGoDistStopCorrect       = cell(1, length(pSignalArray));
pGoTargStopCorrect       = cell(1, length(pSignalArray));
pStopCorrectStopCorrect 	= cell(1, length(pSignalArray));
pStopIncorrectStopCorrect  = cell(1, length(pSignalArray));







% Get default trial selection options
selectOpt = ccm_trial_selection;

% Get trials without respect to signal strength:
selectOpt.ssd = 'none';
selectOpt.outcome     = {'goCorrectTarget'};
goTargTrial = ccm_trial_selection(trialData, selectOpt);
back1GoTargTrial = goTargTrial + 1;
back2GoTargTrial = goTargTrial + 2;
selectOpt.outcome     = {'goCorrectDistractor'};
goDistTrial = ccm_trial_selection(trialData, selectOpt);
back1GoDistTrial = goDistTrial + 1;
goTrial = sort([goTargTrial; goDistTrial]);
back1GoTrial = goTrial + 1;


selectOpt.ssd = 'collapse';
selectOpt.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
stopTargTrial       = ccm_trial_selection(trialData, selectOpt);
back1StopTargTrial  = stopTargTrial + 1;
selectOpt.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
stopDistTrial       = ccm_trial_selection(trialData, selectOpt);
back1StopDistTrial  = stopDistTrial + 1;
stopIncorrectTrial = sort([stopTargTrial; stopDistTrial]);
back1StopIncorrectTrial = stopIncorrectTrial + 1;
selectOpt.outcome       = {'stopCorrect'};
stopStopTrial    = ccm_trial_selection(trialData, selectOpt);
back1StopCorrectTrial = stopStopTrial + 1;










for iPropIndex = 1 : length(pSignalArray);
   iPct = pSignalArray(iPropIndex) * 100;
   selectOpt.rightCheckerPct = iPct;
   
   
   
   % Get trial with respect to signal strength:
   selectOpt.ssd = 'none';
   selectOpt.outcome     = {'goCorrectTarget'};
   iGoTargTrial = ccm_trial_selection(trialData, selectOpt);
   iBack1GoTargTrial = iGoTargTrial + 1;
   selectOpt.outcome     = {'goCorrectDistractor'};
   iGoDistTrial = ccm_trial_selection(trialData, selectOpt);
   iBack1GoDistTrial = iGoDistTrial + 1;
   iGoTrial = sort([iGoTargTrial; iGoDistTrial]);
   iBack1GoTrial = iGoTrial + 1;
   
   selectOpt.ssd = 'collapse';
   selectOpt.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
   iStopTargTrial       = ccm_trial_selection(trialData, selectOpt);
   iBack1StopTargTrial  = iStopTargTrial + 1;
   selectOpt.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
   iStopDistTrial       = ccm_trial_selection(trialData, selectOpt);
   iBack1StopDistTrial  = iStopDistTrial + 1;
   iStopIncorrectTrial = sort([iStopTargTrial; iStopDistTrial]);
   iBack1StopIncorrectTrial = iStopIncorrectTrial + 1;
   selectOpt.outcome       = {'stopCorrect'};
   iStopCorrectTrial    = ccm_trial_selection(trialData, selectOpt);
   iBack1StopCorrectTrial = iStopCorrectTrial + 1;
   
   
   
   
   
   
   % Go Target RT in this signal strength as a function of previous trial
   % outcome in any signal strength
   % ********************************************************************
   % Last trial was Go trial to target
   goTargGoTargTrial{iPropIndex}       = intersect(back1GoTargTrial, iGoTargTrial);
   goTargGoTargRT{iPropIndex}      = allRT(goTargGoTargTrial{iPropIndex});
   goTargRTGoTarg{iPropIndex}      = allRT(goTargGoTargTrial{iPropIndex} - 1);
   % Last trial was Go trial to distractor
   goDistGoTargTrial{iPropIndex}       = intersect(back1GoDistTrial, iGoTargTrial);
   goDistGoTargRT{iPropIndex}      = allRT(goDistGoTargTrial{iPropIndex});
   goDistRTGoTarg{iPropIndex}      = allRT(goDistGoTargTrial{iPropIndex} - 1);
   % Last trial was any Go trial (to target or distractor)
   goGoTargTrial{iPropIndex}       = intersect(back1GoTrial, iGoTargTrial);
   goGoTargRT{iPropIndex}      = allRT(goGoTargTrial{iPropIndex});
   % Last trial was Stop Incorrect, to target
   stopTargGoTargTrial{iPropIndex}     = intersect(back1StopTargTrial, iGoTargTrial);
   stopTargGoTargRT{iPropIndex}    = allRT(stopTargGoTargTrial{iPropIndex});
   % Last trial was Stop Incorrect, to distractor
   stopDistGoTargTrial{iPropIndex}     = intersect(back1StopDistTrial, iGoTargTrial);
   stopDistGoTargRT{iPropIndex}    = allRT(stopDistGoTargTrial{iPropIndex});
   % Last trial was Stop Incorrect, to target or distractor
   stopIncorrectGoTargTrial{iPropIndex}     = intersect(back1StopIncorrectTrial, iGoTargTrial);
   stopIncorrectGoTargRT{iPropIndex}    = allRT(stopIncorrectGoTargTrial{iPropIndex});
   % Last trial was Stop Correct
   stopStopGoTargTrial{iPropIndex}  = intersect(back1StopCorrectTrial, iGoTargTrial);
   stopStopGoTargRT{iPropIndex} = allRT(stopStopGoTargTrial{iPropIndex});
   
   
   
   % Go (to target or distractor) RT i this signal strength as a function
   % of previous trial outcome in any signal strength
   % ********************************************************************
   % Last trial was Go trial to target
   goTargGoTrial{iPropIndex}       = intersect(back1GoTargTrial, iGoTrial);
   goTargGoRT{iPropIndex}      = allRT(goTargGoTrial{iPropIndex});
   % Last trial was Go trial to distractor
   goDistGoTrial{iPropIndex}       = intersect(back1GoDistTrial, iGoTrial);
   goDistGoRT{iPropIndex}      = allRT(goDistGoTrial{iPropIndex});
   % Last trial was any Go trial (to target or distractor)
   goGoTrial{iPropIndex}       = intersect(back1GoTrial, iGoTrial);
   goGoRT{iPropIndex}      = allRT(goGoTrial{iPropIndex});
   % Last trial was Stop Incorrect, to target
   stopTargGoTrial{iPropIndex}     = intersect(back1StopTargTrial, iGoTrial);
   stopTargGoRT{iPropIndex}    = allRT(stopTargGoTrial{iPropIndex});
   % Last trial was Stop Incorrect, to distractor
   stopDistGoTrial{iPropIndex}     = intersect(back1StopDistTrial, iGoTrial);
   stopDistGoRT{iPropIndex}    = allRT(stopDistGoTrial{iPropIndex});
   % Last trial was Stop Incorrect, to target or distractor
   stopIncorrectGoTrial{iPropIndex}     = intersect(back1StopIncorrectTrial, iGoTrial);
   stopIncorrectGoRT{iPropIndex}    = allRT(stopIncorrectGoTrial{iPropIndex});
   % Last trial was Stop Correct
   stopStopGoTrial{iPropIndex}  = intersect(back1StopCorrectTrial, iGoTrial);
   stopStopGoRT{iPropIndex} = allRT(stopStopGoTrial{iPropIndex});
   
   
   
   
   
   
   
   % Trial outcome RTs in this signal strength as a function of previous Go
   % to distractor trial in any signal strength
   % ********************************************************************
   % Next trial was Go trial to target
   %     goDistGoTargTrial{iPropIndex}       = intersect(back1GoDistTrial, goTargTrial);
   %     goDistGoTargRT{iPropIndex}      = allRT(goTargGoTargTrial{iPropIndex});
   % Next trial was Go trial to distractor
   goDistGoDistTrial{iPropIndex}       = intersect(back1GoDistTrial, iGoDistTrial);
   goDistGoDistRT{iPropIndex}      = allRT(goDistGoDistTrial{iPropIndex});
   % Next trial was any Go trial (to target or distractor)
   %     goDistGoTrial{iPropIndex}       = intersect(back1GoDistTrial, goTrial);
   %     goDistGoRT{iPropIndex}      = allRT(goDistGoTrial{iPropIndex});
   % Next trial was Stop Incorrect, to target
   goDistStopTargTrial{iPropIndex}     = intersect(back1GoDistTrial, iStopTargTrial);
   goDistStopTargRT{iPropIndex}    = allRT(goDistStopTargTrial{iPropIndex});
   % Next trial was Stop Incorrect, to distractor
   goDistStopDistTrial{iPropIndex}     = intersect(back1GoDistTrial, iStopDistTrial);
   goDistStopDistRT{iPropIndex}    = allRT(goDistStopDistTrial{iPropIndex});
   % Next trial was Stop Incorrect, to target or distractor
   goDistStopIncorrectTrial{iPropIndex}     = intersect(back1GoDistTrial, iStopIncorrectTrial);
   goDistStopIncorrectRT{iPropIndex}    = allRT(goDistStopIncorrectTrial{iPropIndex});
   %     % Next trial was Stop Correct
   goDistStopCorrectTrial{iPropIndex}  = intersect(back1GoDistTrial, iStopCorrectTrial);
   %     goDistStopCorrectRT{iPropIndex} = allRT(goDistStopCorrectTrial{iPropIndex});
   
   
   
   
   
   
   % Stop probability in this signal strength as a function of previous
   % trial outcome in any signal strength
   % ********************************************************************
   nGoDistStopCorrect      = length(intersect(back1GoDistTrial, iStopCorrectTrial));
   nGoDistStop            = length(intersect(back1GoDistTrial, [iStopIncorrectTrial; iStopCorrectTrial]));
   pGoDistStopCorrect{iPropIndex} = nGoDistStopCorrect / nGoDistStop;
   nGoTargStopCorrect      = length(intersect(back1GoTargTrial, iStopCorrectTrial));
   nGoTargStop             = length(intersect(back1GoTargTrial, [iStopIncorrectTrial; iStopCorrectTrial]));
   pGoTargStopCorrect{iPropIndex} = nGoTargStopCorrect / nGoTargStop;
   nStopCorrectStopCorrect      = length(intersect(back1StopCorrectTrial, iStopCorrectTrial));
   nStopCorrectStop            = length(intersect(back1StopCorrectTrial, [iStopIncorrectTrial; iStopCorrectTrial]));
   pStopCorrectStopCorrect{iPropIndex} = nStopCorrectStopCorrect / nStopCorrectStop;
   nStopIncorrectStopCorrect      = length(intersect(back1StopIncorrectTrial, iStopCorrectTrial));
   nStopIncorrectStop            = length(intersect(back1StopIncorrectTrial, [iStopIncorrectTrial; iStopCorrectTrial]));
   pStopIncorrectStopCorrect{iPropIndex} = nStopIncorrectStopCorrect / nStopIncorrectStop;
   
   
   % ********************************************************************
   % TRIPLETS
   goTargGoTargGoTargTrial{iPropIndex} = intersect(goTargGoTargTrial{iPropIndex}, back2GoTargTrial);
   goTargGoTargGoTargRT{iPropIndex}    = allRT(goTargGoTargGoTargTrial{iPropIndex});
   goTargRTGoTargGoTarg{iPropIndex}    = allRT(goTargGoTargGoTargTrial{iPropIndex} - 2);
   
   goTargGoDistGoTargTrial{iPropIndex} = intersect(goDistGoTargTrial{iPropIndex}, back2GoTargTrial);
   goTargGoDistGoTargRT{iPropIndex}    = allRT(goTargGoDistGoTargTrial{iPropIndex});
   goTargRTGoDistGoTarg{iPropIndex}    = allRT(goTargGoDistGoTargTrial{iPropIndex} - 2);
   
   goTargStopIncorrectGoTargTrial{iPropIndex} = intersect(stopIncorrectGoTargTrial{iPropIndex}, back2GoTargTrial);
   goTargStopIncorrectGoTargRT{iPropIndex}    = allRT(goTargStopIncorrectGoTargTrial{iPropIndex});
   goTargRTStopIncorrectGoTarg{iPropIndex}    = allRT(goTargStopIncorrectGoTargTrial{iPropIndex} - 2);
   
   goTargStopCorrectGoTargTrial{iPropIndex} = intersect(stopStopGoTargTrial{iPropIndex}, back2GoTargTrial);
   goTargStopCorrectGoTargRT{iPropIndex}    = allRT(goTargStopCorrectGoTargTrial{iPropIndex});
   goTargRTStopCorrectGoTarg{iPropIndex}    = allRT(goTargStopCorrectGoTargTrial{iPropIndex} - 2);
   
end




% Go Target RT as a function of previous trial outcome
% ********************************************************************
goTargGoTargRTMean = cellfun(@nanmean, goTargGoTargRT);
goDistGoTargRTMean = cellfun(@nanmean, goDistGoTargRT);
goGoTargRTMean = cellfun(@nanmean, goGoTargRT);
stopTargGoTargRTMean = cellfun(@nanmean, stopTargGoTargRT);
stopDistGoTargRTMean = cellfun(@nanmean, stopDistGoTargRT);
stopIncorrectGoTargRTMean = cellfun(@nanmean, stopIncorrectGoTargRT);
stopStopGoTargRTMean = cellfun(@nanmean, stopStopGoTargRT);

% Go (to target or distractor) RT as a function of previous trial outcome
% ********************************************************************
goTargGoRTMean = cellfun(@nanmean, goTargGoRT);
goDistGoRTMean = cellfun(@nanmean, goDistGoRT);
goGoRTMean = cellfun(@nanmean, goGoRT);
stopTargGoRTMean = cellfun(@nanmean, stopTargGoRT);
stopDistGoRTMean = cellfun(@nanmean, stopDistGoRT);
stopIncorrectGoRTMean = cellfun(@nanmean, stopIncorrectGoRT);
stopStopGoRTMean = cellfun(@nanmean, stopStopGoRT);


% Trial outcomes as a function of previous Go to distractor trial
% ********************************************************************
% goDistGoTargRTMean = cellfun(@nanmean, goDistGoTargRT);
goDistGoDistRTMean = cellfun(@nanmean, goDistGoDistRT);
% goDistGoRTMean = cellfun(@nanmean, goDistGoRT);
goDistStopTargRTMean = cellfun(@nanmean, goDistStopTargRT);
goDistStopDistRTMean = cellfun(@nanmean, goDistStopDistRT);
goDistStopIncorrectRTMean = cellfun(@nanmean, goDistStopIncorrectRT);
% goDistStopCorrectRTMean = cellfun(@nanmean, goDistStopCorrectRT);








% ********************************************************************
% TRIPLETS
goTargGoTargRT1 = cellfun(@(x) x', goTargGoTargRT, 'uniformOutput', false);
goTargRTGoTarg1 = cellfun(@(x) x', goTargRTGoTarg, 'uniformOutput', false);
goTargGoTargGoTargRT = cellfun(@(x) x', goTargGoTargGoTargRT, 'uniformOutput', false);
goTargRTGoTargGoTarg = cellfun(@(x) x', goTargRTGoTargGoTarg, 'uniformOutput', false);
goTargGoDistGoTargRT = cellfun(@(x) x', goTargGoDistGoTargRT, 'uniformOutput', false);
goTargRTGoDistGoTarg = cellfun(@(x) x', goTargRTGoDistGoTarg, 'uniformOutput', false);
goTargStopIncorrectGoTargRT = cellfun(@(x) x', goTargStopIncorrectGoTargRT, 'uniformOutput', false);
goTargRTStopIncorrectGoTarg = cellfun(@(x) x', goTargRTStopIncorrectGoTarg, 'uniformOutput', false);
goTargStopCorrectGoTargRT = cellfun(@(x) x', goTargStopCorrectGoTargRT, 'uniformOutput', false);
goTargRTStopCorrectGoTarg = cellfun(@(x) x', goTargRTStopCorrectGoTarg, 'uniformOutput', false);



% goTargGoTargRTMean = cellfun(@mean, goTargGoTargRT);
% goTargRTGoTargMean = cellfun(@mean, goTargRTGoTarg);
% goDistGoTargRTMean = cellfun(@mean, goDistGoTargRT);
% goDistRTGoTargMean = cellfun(@mean, goDistRTGoTarg);
%


goTargGoTargDiff1 = cellfun(@(x,y) x - y, goTargGoTargRT1, goTargRTGoTarg1, 'uniformOutput', false);
goTargGoTargDiff = cellfun(@(x,y) x - y, goTargGoTargGoTargRT, goTargRTGoTargGoTarg, 'uniformOutput', false);
goDistGoTargDiff = cellfun(@(x,y) x - y, goTargGoDistGoTargRT, goTargRTGoDistGoTarg, 'uniformOutput', false);
stopIncorrectGoTargDiff = cellfun(@(x,y) x - y, goTargStopIncorrectGoTargRT, goTargRTStopIncorrectGoTarg, 'uniformOutput', false);
stopStopGoTargDiff = cellfun(@(x,y) x - y, goTargStopCorrectGoTargRT, goTargRTStopCorrectGoTarg, 'uniformOutput', false);



goTargGoTargDiffMean1 = cellfun(@nanmean, goTargGoTargDiff1);
goTargGoTargDiffMean = cellfun(@nanmean, goTargGoTargDiff);
goDistGoTargDiffMean = cellfun(@nanmean, goDistGoTargDiff);
stopIncorrectGoTargDiffMean = cellfun(@nanmean, stopIncorrectGoTargDiff);
stopStopGoTargDiffMean = cellfun(@nanmean, stopStopGoTargDiff);











% ********************************************************************
%            DO ANALYSIS COLLAPSED AcROSS SIGNAL STRENGTHS
% ********************************************************************

% Peform the same analyses without respect for signal strength
% ********************************************************************

% Trials before Go (to targ or distractor) trials
goTargGoTargRTMeanCollapse = nanmean(cell2mat(goTargGoTargRT'));
goDistGoTargRTMeanCollapse = nanmean(cell2mat(goDistGoTargRT'));
stopTargGoTargRTMeanCollapse = nanmean(cell2mat(stopTargGoTargRT'));
stopDistGoTargRTMeanCollapse = nanmean(cell2mat(stopDistGoTargRT'));
stopStopGoTargRTMeanCollapse = nanmean(cell2mat(stopStopGoTargRT'));

goGoTargRTMeanCollapse          = nanmean(cell2mat(goGoTargRT'));
stopIncorrectGoTargRTMeanCollapse = nanmean(cell2mat(stopIncorrectGoTargRT'));

% Trials before Go (to targ or distractor) trials
goTargGoRTMeanCollapse = nanmean(cell2mat(goTargGoRT'));
goDistGoRTMeanCollapse = nanmean(cell2mat(goDistGoRT'));
stopTargGoRTMeanCollapse = nanmean(cell2mat(stopTargGoRT'));
stopDistGoRTMeanCollapse = nanmean(cell2mat(stopDistGoRT'));
stopStopGoRTMeanCollapse = nanmean(cell2mat(stopStopGoRT'));

goGoRTMeanCollapse          = nanmean(cell2mat(goGoRT'));
stopIncorrectGoRTMeanCollapse = nanmean(cell2mat(stopIncorrectGoRT'));



% Trials after Go Distractor trials
% goDistGoTargRTMeanCollapse = nanmean(cell2mat(goDistGoTargRT));
goDistGoDistRTMeanCollapse = nanmean(cell2mat(goDistGoDistRT'));
goDistStopTargRTMeanCollapse = nanmean(cell2mat(goDistStopTargRT'));
goDistStopDistRTMeanCollapse = nanmean(cell2mat(goDistStopDistRT'));

% goDistGoRTMeanCollapse          = nanmean(cell2mat(goDistGoRT));
goDistStopIncorrectRTMeanCollapse = nanmean(cell2mat(goDistStopIncorrectRT'));




% Stop probability in this signal strength as a function of previous
% trial outcome in any signal strength
% ********************************************************************
nGoDistStopCorrect      = length(intersect(back1GoDistTrial, stopStopTrial));
nGoDistStop            = length(intersect(back1GoDistTrial, [stopIncorrectTrial; stopStopTrial]));
pGoDistStopCorrectCollapse = nGoDistStopCorrect / nGoDistStop;
nGoTargStopCorrect      = length(intersect(back1GoTargTrial, stopStopTrial));
nGoTargStop             = length(intersect(back1GoTargTrial, [stopIncorrectTrial; stopStopTrial]));
pGoTargStopCorrectCollapse = nGoTargStopCorrect / nGoTargStop;
nStopCorrectStopCorrect      = length(intersect(back1StopCorrectTrial, stopStopTrial));
nStopCorrectStop            = length(intersect(back1StopCorrectTrial, [stopIncorrectTrial; stopStopTrial]));
pStopCorrectStopCorrectCollapse = nStopCorrectStopCorrect / nStopCorrectStop;
nStopIncorrectStopCorrect      = length(intersect(back1StopIncorrectTrial, stopStopTrial));
nStopIncorrectStop            = length(intersect(back1StopIncorrectTrial, [stopIncorrectTrial; stopStopTrial]));
pStopIncorrectStopCorrectCollapse = nStopIncorrectStopCorrect / nStopIncorrectStop;









% Trial triplet analyses
goGoBeforeTrial         = intersect(goTargTrial+1, goTargTrial-1) - 1;
goGoAfterTrial          = intersect(goTargTrial+1, goTargTrial-1) + 1;
goDistGoBeforeTrial     = intersect(goDistTrial, intersect(goTargTrial-1, goTargTrial+1)) - 1;
goDistGoAfterTrial      = intersect(goDistTrial, intersect(goTargTrial-1, goTargTrial+1)) + 1;
goStopCorrGoBeforeTrial = intersect(stopStopTrial, intersect(goTargTrial-1, goTargTrial+1)) - 1;
goStopCorrGoAfterTrial  = intersect(stopStopTrial, intersect(goTargTrial-1, goTargTrial+1)) + 1;
goStopIncorrGoBeforeTrial = intersect(stopIncorrectTrial, intersect(goTargTrial-1, goTargTrial+1)) - 1;
goStopIncorrGoAfterTrial = intersect(stopIncorrectTrial, intersect(goTargTrial-1, goTargTrial+1)) + 1;

goGoBeforeRT            = allRT(goGoBeforeTrial);
goGoAfterRT             = allRT(goGoAfterTrial);
goDistGoBeforeRT        = allRT(goDistGoBeforeTrial);
goDistGoAfterRT         = allRT(goDistGoAfterTrial);
goStopCorrGoBeforeRT    = allRT(goStopCorrGoBeforeTrial);
goStopCorrGoAfterRT     = allRT(goStopCorrGoAfterTrial);
goStopIncorrGoBeforeRT  = allRT(goStopIncorrGoBeforeTrial);
goStopIncorrGoAfterRT   = allRT(goStopIncorrGoAfterTrial);






% ********************************************************************
% TRIPLETS
% goTargGoTargDiffCollapse1 = cell(nSession, 1);
% goTargGoTargDiffCollapse = cell(nSession, 1);
% goDistGoTargDiffCollapse = cell(nSession, 1);
% stopIncorrectGoTargDiffCollapse = cell(nSession, 1);
% stopStopGoTargDiffCollapse = cell(nSession, 1);
% for i = 1 : nSession
%     goTargGoTargDiffCollapse1{i} = cell2mat(goTargGoTargDiff1(i,:));
%     goTargGoTargDiffCollapse{i} = cell2mat(goTargGoTargDiff(i,:));
%     goDistGoTargDiffCollapse{i} = cell2mat(goDistGoTargDiff(i,:));
%     stopIncorrectGoTargDiffCollapse{i} = cell2mat(stopIncorrectGoTargDiff(i,:));
%     stopStopGoTargDiffCollapse{i} = cell2mat(stopStopGoTargDiff(i,:));
% end
%
goTargGoTargDiffCollapse1 = nanmean(cell2mat(goTargGoTargDiff1));
goTargGoTargDiffCollapse = nanmean(cell2mat(goTargGoTargDiff));
goDistGoTargDiffCollapse = nanmean(cell2mat(goDistGoTargDiff));
stopIncorrectGoTargDiffCollapse = nanmean(cell2mat(stopIncorrectGoTargDiff));
stopStopGoTargDiffCollapse = nanmean(cell2mat(stopStopGoTargDiff));










% Output the data as a structure:

% Go Target RT as a function of previous trial outcome
% ********************************************************************
data.goTargGoTargRT         = goTargGoTargRT;
data.goTargRTGoTarg         = goTargRTGoTarg;

data.goTargGoTargGoTargRT         = goTargGoTargGoTargRT;
data.goTargRTGoTargGoTarg         = goTargRTGoTargGoTarg;
data.goTargGoDistGoTargRT  = goTargGoDistGoTargRT;
data.goTargRTGoDistGoTarg    = goTargRTGoDistGoTarg;
data.goTargStopIncorrectGoTargRT  = goTargStopIncorrectGoTargRT;
data.goTargStopCorrectGoTargRT    = goTargStopCorrectGoTargRT;
data.goTargRTStopIncorrectGoTarg  = goTargRTStopIncorrectGoTarg;
data.goTargRTStopCorrectGoTarg    = goTargRTStopCorrectGoTarg;

data.goDistGoTargRT         = goDistGoTargRT;
data.goDistRTGoTarg         = goDistRTGoTarg;
data.goGoTargRT             = goGoTargRT;
data.stopTargGoTargRT       = stopTargGoTargRT;
data.stopDistGoTargRT       = stopDistGoTargRT;
data.stopIncorrectGoTargRT  = stopIncorrectGoTargRT;
data.stopStopGoTargRT    = stopStopGoTargRT;

% Go (to target or distractor) RT as a function of previous trial outcome
% ********************************************************************
data.goTargGoRT             = goTargGoRT;
data.goDistGoRT             = goDistGoRT;
data.goGoRT                 = goGoRT;
data.stopTargGoRT           = stopTargGoRT;
data.stopDistGoRT           = stopDistGoRT;
data.stopIncorrectGoRT      = stopIncorrectGoRT;
data.stopStopGoRT        = stopStopGoRT;

% Trial outcomes as a function of previous Go to distractor trial
% ********************************************************************
data.goDistGoDistRT         = goDistGoDistRT;
data.goDistStopTargRT       = goDistStopTargRT;
data.goDistStopDistRT       = goDistStopDistRT;
data.goDistStopIncorrectRT  = goDistStopIncorrectRT;

% Trial outcomes as a function of previous Go to distractor trial
% ********************************************************************
data.pGoDistStopCorrect   	= pGoDistStopCorrect;
data.pGoTargStopCorrect   	= pGoTargStopCorrect;
data.pStopCorrectStopCorrect   	= pStopCorrectStopCorrect;
data.pStopIncorrectStopCorrect   	= pStopIncorrectStopCorrect;
data.pGoDistStopCorrectCollapse   	= pGoDistStopCorrectCollapse;
data.pGoTargStopCorrectCollapse   	= pGoTargStopCorrectCollapse;
data.pStopCorrectStopCorrectCollapse   	= pStopCorrectStopCorrectCollapse;
data.pStopIncorrectStopCorrectCollapse   	= pStopIncorrectStopCorrectCollapse;











if plotFlag
   
   goColor = [0 .7 0];
   goTargColor = [0 .7 0];
   %     goDistColor = [0 .7 0];
   goDistColor = goTargColor ./ 2;
   stopStopColor = [.5 .5 .5];
   stopIncorrectColor = [1 0 0];
   stopTargColor = [1 0 0];
   %     stopDistColor = [1 0 0];
   stopDistColor = stopTargColor ./ 2;
   ylimRange = [200 700];
   
   % *****************************************************************
   % Go Target preceded by various outcomes, broken down by
   % target/distractor choice
   plot(ax(axXGoTarg), pSignalArray, goTargGoTargRTMean, '-o', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axXGoTarg), pSignalArray, goDistGoTargRTMean, '--o', 'color', goDistColor, 'linewidth', 1, 'markerfacecolor', goDistColor, 'markeredgecolor', goDistColor)
   plot(ax(axXGoTarg), pSignalArray, stopTargGoTargRTMean, '-o', 'color', stopTargColor, 'linewidth', 1, 'markerfacecolor', stopTargColor, 'markeredgecolor', stopTargColor)
   plot(ax(axXGoTarg), pSignalArray, stopDistGoTargRTMean, '--o', 'color', stopDistColor, 'linewidth', 1, 'markerfacecolor', stopDistColor, 'markeredgecolor', stopDistColor)
   plot(ax(axXGoTarg), pSignalArray, stopStopGoTargRTMean, '-o', 'color', stopStopColor, 'linewidth', 1, 'markerfacecolor', stopStopColor, 'markeredgecolor', stopStopColor)
   legend(ax(axXGoTarg), {'goTargGoTarg', 'goDistGoTarg', 'stopTargGoTarg', 'stopDistGoTarg', 'stopStopGoTarg'})
   
   plot(ax(axXGoTarg), pSignalArray(1) - choicePlotXMargin+.01, goTargGoTargRTMeanCollapse, '-d', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axXGoTarg), pSignalArray(1) - choicePlotXMargin+.01, goDistGoTargRTMeanCollapse, '--d', 'color', goDistColor, 'linewidth', 1, 'markerfacecolor', goDistColor, 'markeredgecolor', goDistColor)
   plot(ax(axXGoTarg), pSignalArray(1) - choicePlotXMargin+.01, stopTargGoTargRTMeanCollapse, '-d', 'color', stopTargColor, 'linewidth', 1, 'markerfacecolor', stopTargColor, 'markeredgecolor', stopTargColor)
   plot(ax(axXGoTarg), pSignalArray(1) - choicePlotXMargin+.01, stopDistGoTargRTMeanCollapse, '--d', 'color', stopDistColor, 'linewidth', 1, 'markerfacecolor', stopDistColor, 'markeredgecolor', stopDistColor)
   plot(ax(axXGoTarg), pSignalArray(1) - choicePlotXMargin+.01, stopStopGoTargRTMeanCollapse, '-d', 'color', stopStopColor, 'linewidth', 1, 'markerfacecolor', stopStopColor, 'markeredgecolor', stopStopColor)
   % errorbar(ax(axXGoTarg), signalStrengthLeft ,goLeftToTargMean, goLeftToTargStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
   
   set(ax(axXGoTarg), 'xtick', pSignalArray)
   set(ax(axXGoTarg), 'xtickLabel', pSignalArray*100)
   set(get(ax(axXGoTarg), 'ylabel'), 'String', 'RT')
   set(ax(axXGoTarg),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
   set(ax(axXGoTarg),'YLim',ylimRange)
   plot(ax(axXGoTarg), [.5 .5], ylim, '--k')
   title(ax(axXGoTarg), 'trial then GoTarget')
   
   
   
   % Go Target preceded by various outcomes, collapsed target/distractor choice
   plot(ax(axXGoTarg2), pSignalArray, goGoTargRTMean, '-o', 'color', goColor, 'linewidth', 1, 'markerfacecolor', goColor, 'markeredgecolor', goColor)
   plot(ax(axXGoTarg2), pSignalArray, stopIncorrectGoTargRTMean, '-o', 'color', stopIncorrectColor, 'linewidth', 1, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
   plot(ax(axXGoTarg2), pSignalArray, stopStopGoTargRTMean, '-o', 'color', stopStopColor, 'linewidth', 1, 'markerfacecolor', stopStopColor, 'markeredgecolor', stopStopColor)
   legend(ax(axXGoTarg2), {'goGoTarg', 'stopIncorrectGoTarg', 'stopStopGoTarg'})
   
   plot(ax(axXGoTarg2), pSignalArray(1) - choicePlotXMargin+.01, goGoTargRTMeanCollapse, '-d', 'color', goColor, 'linewidth', 1, 'markerfacecolor', goColor, 'markeredgecolor', goColor)
   plot(ax(axXGoTarg2), pSignalArray(1) - choicePlotXMargin+.01, stopIncorrectGoTargRTMeanCollapse, '-d', 'color', stopIncorrectColor, 'linewidth', 1, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
   plot(ax(axXGoTarg2), pSignalArray(1) - choicePlotXMargin+.01, stopStopGoTargRTMeanCollapse, '-d', 'color', stopStopColor, 'linewidth', 1, 'markerfacecolor', stopStopColor, 'markeredgecolor', stopStopColor)
   
   % errorbar(ax(axXGoTarg2), signalStrengthLeft ,goLeftToTargMean, goLeftToTargStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
   
   set(ax(axXGoTarg2), 'xtick', pSignalArray)
   set(ax(axXGoTarg2), 'xtickLabel', pSignalArray*100)
   set(get(ax(axXGoTarg2), 'ylabel'), 'String', 'RT')
   set(ax(axXGoTarg2),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
   set(ax(axXGoTarg2),'YLim',ylimRange)
   plot(ax(axXGoTarg2), [.5 .5], ylim, '--k')
   title(ax(axXGoTarg2), 'trial then GoTarget')
   
   
   
   
   
   
   % *****************************************************************
   % Go (to target or distractor) preceded by various outcomes, broken down by
   % target/distractor choice
   plot(ax(axXGo), pSignalArray, goTargGoRTMean, '-o', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axXGo), pSignalArray, goDistGoRTMean, '--o', 'color', goDistColor, 'linewidth', 1, 'markerfacecolor', goDistColor, 'markeredgecolor', goDistColor)
   plot(ax(axXGo), pSignalArray, stopTargGoRTMean, '-o', 'color', stopTargColor, 'linewidth', 1, 'markerfacecolor', stopTargColor, 'markeredgecolor', stopTargColor)
   plot(ax(axXGo), pSignalArray, stopDistGoRTMean, '--o', 'color', stopDistColor, 'linewidth', 1, 'markerfacecolor', stopDistColor, 'markeredgecolor', stopDistColor)
   plot(ax(axXGo), pSignalArray, stopStopGoRTMean, '-o', 'color', stopStopColor, 'linewidth', 1, 'markerfacecolor', stopStopColor, 'markeredgecolor', stopStopColor)
   legend(ax(axXGo), {'goTargGo', 'goDistGo', 'stopTargGo', 'stopDistGo', 'stopStopGo'})
   
   plot(ax(axXGo), pSignalArray(1) - choicePlotXMargin+.01, goTargGoRTMeanCollapse, '-d', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axXGo), pSignalArray(1) - choicePlotXMargin+.01, goDistGoRTMeanCollapse, '--d', 'color', goDistColor, 'linewidth', 1, 'markerfacecolor', goDistColor, 'markeredgecolor', goDistColor)
   plot(ax(axXGo), pSignalArray(1) - choicePlotXMargin+.01, stopTargGoRTMeanCollapse, '-d', 'color', stopTargColor, 'linewidth', 1, 'markerfacecolor', stopTargColor, 'markeredgecolor', stopTargColor)
   plot(ax(axXGo), pSignalArray(1) - choicePlotXMargin+.01, stopDistGoRTMeanCollapse, '--d', 'color', stopDistColor, 'linewidth', 1, 'markerfacecolor', stopDistColor, 'markeredgecolor', stopDistColor)
   plot(ax(axXGo), pSignalArray(1) - choicePlotXMargin+.01, stopStopGoRTMeanCollapse, '-d', 'color', stopStopColor, 'linewidth', 1, 'markerfacecolor', stopStopColor, 'markeredgecolor', stopStopColor)
   % errorbar(ax(axXGo), signalStrengthLeft ,goLeftToTargMean, goLeftToTargStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
   
   set(ax(axXGo), 'xtick', pSignalArray)
   set(ax(axXGo), 'xtickLabel', pSignalArray*100)
   set(get(ax(axXGo), 'ylabel'), 'String', 'RT')
   set(ax(axXGo),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
   %         set(ax(axXGo),'YLim',[min([goTargGoRTMean; stopTargGoRTMean; stopStopGoRTMean]) - 25, max([goTargGoRTMean; stopTargGoRTMean; stopStopGoRTMean]) + 25])
   set(ax(axXGo),'YLim',ylimRange)
   plot(ax(axXGo), [.5 .5], ylim, '--k')
   title(ax(axXGo), 'trial then Go')
   
   
   
   % Go (to target or distractor) preceded by various outcomes, collapsed target/distractor choice
   plot(ax(axXGo2), pSignalArray, goGoRTMean, '-o', 'color', goColor, 'linewidth', 1, 'markerfacecolor', goColor, 'markeredgecolor', goColor)
   plot(ax(axXGo2), pSignalArray, stopIncorrectGoRTMean, '-o', 'color', stopIncorrectColor, 'linewidth', 1, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
   plot(ax(axXGo2), pSignalArray, stopStopGoRTMean, '-o', 'color', stopStopColor, 'linewidth', 1, 'markerfacecolor', stopStopColor, 'markeredgecolor', stopStopColor)
   legend(ax(axXGo2), {'goGo', 'stopIncorrectGo', 'stopStopGo'})
   
   plot(ax(axXGo2), pSignalArray(1) - choicePlotXMargin+.01, goGoRTMeanCollapse, '-d', 'color', goColor, 'linewidth', 1, 'markerfacecolor', goColor, 'markeredgecolor', goColor)
   plot(ax(axXGo2), pSignalArray(1) - choicePlotXMargin+.01, stopIncorrectGoRTMeanCollapse, '-d', 'color', stopIncorrectColor, 'linewidth', 1, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
   plot(ax(axXGo2), pSignalArray(1) - choicePlotXMargin+.01, stopStopGoRTMeanCollapse, '-d', 'color', stopStopColor, 'linewidth', 1, 'markerfacecolor', stopStopColor, 'markeredgecolor', stopStopColor)
   % errorbar(ax(axXGo2), signalStrengthLeft ,goLeftToTargMean, goLeftToTargStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
   
   set(ax(axXGo2), 'xtick', pSignalArray)
   set(ax(axXGo2), 'xtickLabel', pSignalArray*100)
   set(get(ax(axXGo2), 'ylabel'), 'String', 'RT')
   set(ax(axXGo2),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
   %         set(ax(axXGo2),'YLim',[min([goTargGoRTMean; stopTargGoRTMean; stopStopGoRTMean]) - 25, max([goTargGoRTMean; stopTargGoRTMean; stopStopGoRTMean]) + 25])
   set(ax(axXGo2),'YLim',ylimRange)
   plot(ax(axXGo2), [.5 .5], ylim, '--k')
   title(ax(axXGo2), 'trial then Go')
   
   
   
   
   
   
   % *****************************************************************
   % Go distractor followed by various outcomes, broken down by
   % target/distractor choice
   plot(ax(axGoDistX), pSignalArray, goDistGoTargRTMean, '-o', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axGoDistX), pSignalArray, goDistGoDistRTMean, '--o', 'color', goDistColor, 'linewidth', 1, 'markerfacecolor', goDistColor, 'markeredgecolor', goDistColor)
   plot(ax(axGoDistX), pSignalArray, goDistStopTargRTMean, '-o', 'color', stopTargColor, 'linewidth', 1, 'markerfacecolor', stopTargColor, 'markeredgecolor', stopTargColor)
   plot(ax(axGoDistX), pSignalArray, goDistStopDistRTMean, '--o', 'color', stopDistColor, 'linewidth', 1, 'markerfacecolor', stopDistColor, 'markeredgecolor', stopDistColor)
   legend(ax(axGoDistX), {'goDistGoTarg', 'goDistgoDist', 'goDistStopTarg', 'goDistStopDist'})
   
   plot(ax(axGoDistX), pSignalArray(1) - choicePlotXMargin+.01, goDistGoTargRTMeanCollapse, '-d', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axGoDistX), pSignalArray(1) - choicePlotXMargin+.01, goDistGoDistRTMeanCollapse, '--d', 'color', goDistColor, 'linewidth', 1, 'markerfacecolor', goDistColor, 'markeredgecolor', goDistColor)
   plot(ax(axGoDistX), pSignalArray(1) - choicePlotXMargin+.01, goDistStopTargRTMeanCollapse, '-d', 'color', stopTargColor, 'linewidth', 1, 'markerfacecolor', stopTargColor, 'markeredgecolor', stopTargColor)
   plot(ax(axGoDistX), pSignalArray(1) - choicePlotXMargin+.01, goDistStopDistRTMeanCollapse, '--d', 'color', stopDistColor, 'linewidth', 1, 'markerfacecolor', stopDistColor, 'markeredgecolor', stopDistColor)
   % errorbar(ax(axGoDistX), signalStrengthLeft ,goLeftToTargMean, goLeftToTargStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
   
   set(ax(axGoDistX), 'xtick', pSignalArray)
   set(ax(axGoDistX), 'xtickLabel', pSignalArray*100)
   set(get(ax(axGoDistX), 'ylabel'), 'String', 'RT')
   set(ax(axGoDistX),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
   %         set(ax(axGoDistX),'YLim',[min([goTargGoRTMean; stopTargGoRTMean; stopStopGoRTMean]) - 25, max([goTargGoRTMean; stopTargGoRTMean; stopStopGoRTMean]) + 25])
   set(ax(axGoDistX),'YLim',ylimRange)
   plot(ax(axGoDistX), [.5 .5], ylim, '--k')
   title(ax(axGoDistX), 'Go Distractor then trial')
   
   
   
   % Go distractor followed by various outcomes, collapsed target/distractor choice
   plot(ax(axGoDistX2), pSignalArray, goDistGoRTMean, '-o', 'color', goColor, 'linewidth', 1, 'markerfacecolor', goColor, 'markeredgecolor', goColor)
   plot(ax(axGoDistX2), pSignalArray, goDistStopIncorrectRTMean, '-o', 'color', stopIncorrectColor, 'linewidth', 1, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
   legend(ax(axGoDistX2), {'goDistGo', 'goDistStopIncorrect'})
   
   plot(ax(axGoDistX2), pSignalArray(1) - choicePlotXMargin+.01, goDistGoRTMeanCollapse, '-d', 'color', goColor, 'linewidth', 1, 'markerfacecolor', goColor, 'markeredgecolor', goColor)
   plot(ax(axGoDistX2), pSignalArray(1) - choicePlotXMargin+.01, goDistStopIncorrectRTMeanCollapse, '-d', 'color', stopIncorrectColor, 'linewidth', 1, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
   % errorbar(ax(axGoDistX2), signalStrengthLeft ,goLeftToTargMean, goLeftToTargStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
   
   set(ax(axGoDistX2), 'xtick', pSignalArray)
   set(ax(axGoDistX2), 'xtickLabel', pSignalArray*100)
   set(get(ax(axGoDistX2), 'ylabel'), 'String', 'RT')
   set(ax(axGoDistX2),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
   %         set(ax(axGoDistX2),'YLim',[min([goTargGoRTMean{:}; stopTargGoRTMean{:}; stopStopGoRTMean{:}]) - 25, max([goTargGoRTMean{:}; stopTargGoRTMean{:}; stopStopGoRTMean{:}]) + 25])
   set(ax(axGoDistX2),'YLim',ylimRange)
   plot(ax(axGoDistX2), [.5 .5], ylim, '--k')
   title(ax(axGoDistX2), 'Go Distractor then trial')
   
   
   
   
   % *****************************************************************
   % Stopping probability as a function of previous trial outcomes
   plot(ax(axStopP), pSignalArray, cell2mat(pGoTargStopCorrect), '-o', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axStopP), pSignalArray, cell2mat(pGoDistStopCorrect), '--o', 'color', goDistColor, 'linewidth', 1, 'markerfacecolor', goDistColor, 'markeredgecolor', goDistColor)
   plot(ax(axStopP), pSignalArray, cell2mat(pStopCorrectStopCorrect), '-o', 'color', stopStopColor, 'linewidth', 1, 'markerfacecolor', stopStopColor, 'markeredgecolor', stopStopColor)
   plot(ax(axStopP), pSignalArray, cell2mat(pStopIncorrectStopCorrect), '--o', 'color', stopIncorrectColor, 'linewidth', 1, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
   legend(ax(axStopP), {'pGoTargStopCorrect', 'pGoDistStopCorrect', 'pStopCorrectStopCorrect', 'pStopIncorrectStopCorrect'})
   
   plot(ax(axStopP), pSignalArray(1) - choicePlotXMargin+.01, pGoTargStopCorrectCollapse, '-d', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axStopP), pSignalArray(1) - choicePlotXMargin+.01, pGoDistStopCorrectCollapse, '--d', 'color', goDistColor, 'linewidth', 1, 'markerfacecolor', goDistColor, 'markeredgecolor', goDistColor)
   plot(ax(axStopP), pSignalArray(1) - choicePlotXMargin+.01, pStopCorrectStopCorrectCollapse, '-d', 'color', stopStopColor, 'linewidth', 1, 'markerfacecolor', stopStopColor, 'markeredgecolor', stopStopColor)
   plot(ax(axStopP), pSignalArray(1) - choicePlotXMargin+.01, pStopIncorrectStopCorrectCollapse, '--d', 'color', stopIncorrectColor, 'linewidth', 1, 'markerfacecolor', stopIncorrectColor, 'markeredgecolor', stopIncorrectColor)
   % errorbar(ax(axXGoTarg), signalStrengthLeft ,goLeftToTargMean, goLeftToTargStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
   
   set(ax(axStopP), 'xtick', pSignalArray)
   set(ax(axStopP), 'xtickLabel', pSignalArray*100)
   set(get(ax(axStopP), 'ylabel'), 'String', 'RT')
   set(ax(axStopP),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
   set(ax(axStopP),'YLim',[0 1])
   plot(ax(axStopP), [.5 .5], ylim, '--k')
   title(ax(axStopP), 'Stop Probability')
   
   
   
   
   
   
   % *****************************************************************
   % TRIPLETS
   plot(ax(axTriplet), pSignalArray, goTargGoTargDiffMean1, '-*', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axTriplet), pSignalArray, goTargGoTargDiffMean, '-o', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axTriplet), pSignalArray, goDistGoTargDiffMean, '--o', 'color', goDistColor, 'linewidth', 1, 'markerfacecolor', goDistColor, 'markeredgecolor', goDistColor)
   plot(ax(axTriplet), pSignalArray, stopIncorrectGoTargDiffMean, '-o', 'color', stopTargColor, 'linewidth', 1, 'markerfacecolor', stopTargColor, 'markeredgecolor', stopTargColor)
   plot(ax(axTriplet), pSignalArray, stopStopGoTargDiffMean, '--o', 'color', stopDistColor, 'linewidth', 1, 'markerfacecolor', stopDistColor, 'markeredgecolor', stopDistColor)
   legend(ax(axTriplet), {'goTarg2', 'goTargGoTarg', 'goDistGoTarg', 'stopTargGoTarg', 'stopCorrGoTarg'})
   
   plot(ax(axTriplet), pSignalArray(1) - choicePlotXMargin+.01, goTargGoTargDiffCollapse1, '-*', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axTriplet), pSignalArray(1) - choicePlotXMargin+.01, goTargGoTargDiffCollapse, '-d', 'color', goTargColor, 'linewidth', 1, 'markerfacecolor', goTargColor, 'markeredgecolor', goTargColor)
   plot(ax(axTriplet), pSignalArray(1) - choicePlotXMargin+.01, goDistGoTargDiffCollapse, '--d', 'color', goDistColor, 'linewidth', 1, 'markerfacecolor', goDistColor, 'markeredgecolor', goDistColor)
   plot(ax(axTriplet), pSignalArray(1) - choicePlotXMargin+.01, stopIncorrectGoTargDiffCollapse, '-d', 'color', stopTargColor, 'linewidth', 1, 'markerfacecolor', stopTargColor, 'markeredgecolor', stopTargColor)
   plot(ax(axTriplet), pSignalArray(1) - choicePlotXMargin+.01, stopStopGoTargDiffCollapse, '--d', 'color', stopDistColor, 'linewidth', 1, 'markerfacecolor', stopDistColor, 'markeredgecolor', stopDistColor)
   % errorbar(ax(axTriplet), signalStrengthLeft ,goLeftToTargMean, goLeftToTargStd, '.' , 'linestyle' , 'none', 'color', goColor, 'linewidth' , 1)
   
   set(ax(axTriplet), 'xtick', pSignalArray)
   set(ax(axTriplet), 'xtickLabel', pSignalArray*100)
   set(get(ax(axTriplet), 'ylabel'), 'String', 'RT')
   set(ax(axTriplet),'XLim',[pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
   %     set(ax(axTriplet),'YLim',[min([goTargGoRTMean; stopTargGoRTMean; stopStopGoRTMean]) - 25, max([goTargGoRTMean; stopTargGoRTMean; stopStopGoRTMean]) + 25])
   set(ax(axTriplet),'YLim',[-150 150])
   plot(ax(axTriplet), [.5 .5], ylim, '--k')
   title(ax(axTriplet), 'TRIPLETS')
   
   
end

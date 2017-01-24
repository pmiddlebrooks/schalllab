
%%
subjectID = 'Broca';
sessionID = '_concat';
% sessionID = 'bp093n02';

% subjectID = 'Xena';
% sessionID = '_concat';

% subjectID = 'hu';
% sessionID = 'Allsaccade';
%

% % Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
rSignalArray = ExtraVar.pSignalArray;
localFigurePath = local_figure_path;
accuracyOrResponse = 'response';
% accuracyOrResponse = 'accuracy';

% Which signal levels (difficulty conditions) to model?
% ---------------------------------------
signalModel = 'all';
switch signalModel
   case 'all'
      % Use all (right and left) trials
      rSignalArray = rSignalArray(rSignalArray ~= .5);
   case 'right'
      % Rightward-only
      rSignalArray = rSignalArray(rSignalArray > .5);
   case 'left'
      % Leftward-only
      rSignalArray = rSignalArray(rSignalArray < .5);
end

nSignal = length(rSignalArray);

%         allRT = trialData.responseOnset - trialData.responseCueOn;







aInd    = 1;  % boundary separation
TerInd  = 2;  % nondecision time
etaInd  = 3;  % inter-trial standard deviation of drift rate
zInd    = 4;  % starting point
szInd   = 5;  % inter-trial range of starting point
stInd   = 6;  % inter-trial range of non-decision time
vInd    = 7;  % drift rate
piInd   = 8;  % proportion of the data that are valid (non-outliers)
gammaInd = 9; % proportion of the data that are guesses

MIN_RT = 120;
MAX_RT = 1200;
STD_MULTIPLE = 3;
% Get rid of trials with outlying RTs
[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTD);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);

nTrial = size(trialData.trialOutcome, 1);
% Get default trial selection options
selectOpt       = ccm_trial_selection;
selectOpt.ssd       = 'none';
selectOpt.rightCheckerPct = rSignalArray .* 100;


selectOpt.outcome = {'goCorrectTarget', 'targetHoldAbort', 'goCorrectDistractor', 'distractorHoldAbort'};
selectOpt.targDir       = 'collapse';
trial = ccm_trial_selection(trialData, selectOpt);

selectOpt.outcome = {'goCorrectTarget', 'targetHoldAbort'};
selectOpt.targDir       = 'right';
trialRightCorr      = ccm_trial_selection(trialData, selectOpt);
selectOpt.outcome = {'goCorrectDistractor', 'distractorHoldAbort'};
selectOpt.targDir       = 'left';
trialLeftIncorr     = ccm_trial_selection(trialData, selectOpt);
rightResponse = zeros(nTrial, 1);
rightResponse([trialRightCorr; trialLeftIncorr]) = 1;

%         goTargTrial = ccm_trial_selection(trialData, {'goCorrectTarget'; 'targetHoldAbort'}, iPct, ssdRange, targetHemifield);
%         goDistTrial = ccm_trial_selection(trialData, {'goCorrectDistractor', 'distractorHoldAbort'}, iPct, ssdRange, targetHemifield);

% condition = nan(nTrial, 1);
% for i = 1 : length(targ1PropArray)
%     condition(trialData.targ1CheckerProp == targ1PropArray(i)) = i;
% end

% % Form a data matrix
% accuracy = ones(length(goTargTrial), 1);
% rt = trialData.responseOnset(goTargTrial) - trialData.responseCueOn(goTargTrial) ./1000;
% % targData = [trialData.targ1CheckerProp(goTargTrial), accuracy, rt];
% targData = [condition(goTargTrial), accuracy, rt];
%
%
% accuracy = zeros(length(goDistTrial), 1);
% rt = trialData.responseOnset(goDistTrial) - trialData.responseCueOn(goDistTrial) ./1000;
% % distData = [trialData.targ1CheckerProp(goDistTrial), accuracy, rt];
% distData = [condition(goDistTrial), accuracy, rt];

% For inputting accuracy instead of response direction
accuracy = strcmp(trialData.trialOutcome, 'goCorrectTarget') | strcmp(trialData.trialOutcome, 'targetHoldAbort');
switch accuracyOrResponse
   case 'accuracy'
      data = [trialData.targ1CheckerProp(trial), accuracy(trial), .001 .* trialData.rt(trial)];
   case 'response'
      % For inputting response direction instead of accuracy
      data = [trialData.targ1CheckerProp(trial), rightResponse(trial), .001 .* trialData.rt(trial)];
end






% =========================================================================
% DESIGN MATRIX
% get default
options = multiestv4;
% How many diffrent design matrices are we going to test?
if strcmp(signalModel, 'all')
   nDesign = 4;
else
   nDesign = 3;
end
% copy this design matrix, so we can make incremental changes to it
options = repmat(options, nDesign, 1);

O = ones(nSignal, 1);  % for fixed paramters
I = eye(nSignal); % for free parameters (free btwn all conditions
H = [repmat([1 0], nSignal/2, 1); repmat([0 1], nSignal/2, 1)];% for free parameters btwn response hemifields

% All fixed
designMatrixl = {O,O,O,O,O,O,O,O,O};
% Drift all free
designMatrix2 = {O,O,O,O,O,O,I,O,O};


% Non-decision time free between hemifields
designMatrix3 = {O,H,O,O,O,O,I,O,O};


% Response boundary free between hemifields
designMatrix4 = {H,H,O,O,O,O,I,O,O};


options(1).DesignMatrix = designMatrixl;
options(1).Name = 'designMatrix1';
options(2).DesignMatrix = designMatrix2;
options(2).Name = 'designMatrix2';
options(3).DesignMatrix = designMatrix3;
options(3).Name = 'designMatrix3';
options(4).DesignMatrix = designMatrix4;
options(4).Name = 'designMatrix4';


% If modeling all directions, include response bias as a parameter
if strcmp(signalModel, 'all')
   
   
   
   % Response bias free between hemifields
   designMatrix5 = {H,H,O,H,O,O,I,O,O};
   
   
   options(5).DesignMatrix = designMatrix5;
   options(5).Name = 'designMatrix5';
   
end


%
% % Ter Nondecision time free between response hemifields
% designMatrix3 = {O,H,O,O,O,O,I,O,O};
%
% % Ter Nondecision time free between response hemifields
% designMatrix4b = {H,H,O,O,O,O,I,O,O};
%
% % Response boundary all free
% designMatrix5a = {I,O,O,O,O,O,I,O,O};
% % ... with Ter Nondecision time free between response hemifields
% designMatrix5b = {I,H,O,O,O,O,I,O,O};
%
% % Response bias free between hemifields
% designMatrix6a = {O,O,O,H,O,O,I,O,O};
% designMatrix6b = {O,H,O,H,O,O,I,O,O};
% designMatrix6c = {H,O,O,H,O,O,I,O,O};
% designMatrix6d = {H,H,O,H,O,O,I,O,O};
% designMatrix6e = {I,O,O,H,O,O,I,O,O};
% designMatrix6f = {I,H,O,H,O,O,I,O,O};
%
% % Response bias all free
% designMatrix7a = {O,O,O,I,O,O,I,O,O};
% designMatrix7b = {O,H,O,I,O,O,I,O,O};
% designMatrix7c = {H,O,O,I,O,O,I,O,O};
% designMatrix7d = {H,H,O,I,O,O,I,O,O};
% designMatrix7e = {I,O,O,I,O,O,I,O,O};
% designMatrix7f = {I,H,O,I,O,O,I,O,O};
%
%
%
% options(1).DesignMatrix = designMatrixl;
% options(1).Name = 'designMatrixl';
% options(2).DesignMatrix = designMatrix2;
% options(2).Name = 'designMatrix2';
% options(3).DesignMatrix = designMatrix3;
% options(3).Name = 'designMatrix3';
% options(4).DesignMatrix = designMatrix4a;
% options(4).Name = 'designMatrix4a';
% options(5).DesignMatrix = designMatrix4b;
% options(5).Name = 'designMatrix4b';
% options(6).DesignMatrix = designMatrix5a;
% options(6).Name = 'designMatrix5a';
% options(7).DesignMatrix = designMatrix5b;
% options(7).Name = 'designMatrix5b';
% options(8).DesignMatrix = designMatrix6a;
% options(8).Name = 'designMatrix6a';
% options(9).DesignMatrix = designMatrix6b;
% options(9).Name = 'designMatrix6b';
% options(10).DesignMatrix = designMatrix6c;
% options(10).Name = 'designMatrix6c';
% options(11).DesignMatrix = designMatrix6d;
% options(11).Name = 'designMatrix6d';
% options(12).DesignMatrix = designMatrix6e;
% options(12).Name = 'designMatrix6e';
% options(13).DesignMatrix = designMatrix6f;
% options(13).Name = 'designMatrix6f';
% options(14).DesignMatrix = designMatrix7a;
% options(14).Name = 'designMatrix7a';
% options(15).DesignMatrix = designMatrix7b;
% options(15).Name = 'designMatrix7b';
% options(16).DesignMatrix = designMatrix7c;
% options(16).Name = 'designMatrix7c';
% options(17).DesignMatrix = designMatrix7d;
% options(17).Name = 'designMatrix7d';
% options(18).DesignMatrix = designMatrix7e;
% options(18).Name = 'designMatrix7e';
% options(19).DesignMatrix = designMatrix7f;
% options(19).Name = 'designMatrix7f';


output = multiestv4(data,options)

qtable(output)
tbl = qtable(output)


%%
% =========================================================================
% PLOTS AND READ-OUTS



% -------------------------------------------------------------------------
% Parameters plots
cDrift = 1;
cBnd = 3;
cTer = 2;
nColumn = 3;

nRow = size(options, 1);


figureHandle = 1;
driftMax = .9;
driftMin = -driftMax;
bndMax = .1;
biasMax = bndMax;
TerMax = .2;
xLimVec = [0 nSignal + 1];
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
clf
for hOut = 1 : size(options, 1)
   
   
   
   ax(hOut, cDrift) = axes('units', 'centimeters', 'position', [xAxesPosition(hOut, cDrift) yAxesPosition(hOut, cDrift) axisWidth axisHeight]);
   set(ax(hOut, cDrift), 'ylim', [driftMin driftMax], 'xlim', xLimVec)
   hold(ax(hOut, cDrift), 'on')
   if hOut == 1
      title('Drift (v)')
   end
   ylabel(output(hOut).Name)
   paramplot(output(hOut).Minimum(:,vInd),output(hOut).StdErr(:,vInd), ax(hOut, cDrift))
   plot([1 6], [0 0], '--k')
   
   ax(hOut, cBnd) = axes('units', 'centimeters', 'position', [xAxesPosition(hOut, cBnd) yAxesPosition(hOut, cBnd) axisWidth axisHeight]);
   set(ax(hOut, cBnd), 'ylim', [0 bndMax], 'xlim', xLimVec)
   if hOut == 1
      title('Bound and Bias: (a and z)')
   end
   hold(ax(hOut, cBnd), 'on')
   paramplot(output(hOut).Minimum(:,aInd),output(hOut).StdErr(:,aInd), ax(hOut, cBnd))
   paramplot(output(hOut).Minimum(:,zInd),output(hOut).StdErr(:,zInd), ax(hOut, cBnd))
   
   
   %     ax(hOut, cBias) = axes('units', 'centimeters', 'position', [xAxesPosition(hOut, cBias) yAxesPosition(hOut, cBias) axisWidth axisHeight]);
   %     set(ax(hOut, cBias), 'ylim', [0 biasMax], 'xlim', xLimVec)
   %     hold(ax(hOut, cBias), 'on')
   %     paramplot(output(hOut).Minimum(:,zInd),output(hOut).StdErr(:,zInd), ax(hOut, cBias))
   
   ax(hOut, cTer) = axes('units', 'centimeters', 'position', [xAxesPosition(hOut, cTer) yAxesPosition(hOut, cTer) axisWidth axisHeight]);
   set(ax(hOut, cTer), 'ylim', [0 TerMax], 'xlim', xLimVec)
   if hOut == 1
      title('Non-decision time (Ter)')
   end
   hold(ax(hOut, cTer), 'on')
   paramplot(output(hOut).Minimum(:,TerInd),output(hOut).StdErr(:,TerInd), ax(hOut, cTer))
   
   
end % hOut
print(figureHandle,[localFigurePath, subjectID, '_', sessionID, '_', signalModel, '_', accuracyOrResponse],'-dpdf', '-r300')




% CDFs PLOT
rLeft = 1;
rRight = 2;
nColumn = nSignal;
nRow = 2;

figureHandle = 10;
yLimVec = [0 1];
xLimVec = [.1 .9];

leftResp = 0;
rightResp = 1;
for hOut = 1 : size(options, 1)
   figureHandle = figureHandle + 1;
   
   [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
   % [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = model_figure(nRow, nColumn, figureHandle, 'save');
   clf
   for iSig = 1 : nSignal
      sig = rSignalArray(iSig);
      
      
      
      ax(rLeft, iSig) = axes('units', 'centimeters', 'position', [xAxesPosition(rLeft, iSig) yAxesPosition(rLeft, iSig) axisWidth axisHeight]);
      set(ax(rLeft, iSig), 'ylim', yLimVec, 'xlim', xLimVec)
      if iSig == 1
         title(output(hOut).Name)
      end
      
      %         hold(ax(rLeft, iSig), 'on')
      
      rtLeft = data(data(:,1) == sig & data(:,2) == leftResp, 3);
      edfcdf(rtLeft, sum(data(:,1) == sig), leftResp,...
         output(hOut).Minimum(iSig,:), output(hOut).Options.nQ(iSig,:), ax(rLeft, iSig))
      
      
      
      ax(rRight, iSig) = axes('units', 'centimeters', 'position', [xAxesPosition(rRight, iSig) yAxesPosition(rRight, iSig) axisWidth axisHeight]);
      set(ax(rRight, iSig), 'ylim', yLimVec, 'xlim', xLimVec)
      %         hold(ax(rRight, iSig), 'on')
      
      rtRight = data(data(:,1) == sig & data(:,2) == rightResp, 3);
      edfcdf(rtRight, sum(data(:,1) == sig), rightResp,...
         output(hOut).Minimum(iSig,:), output(hOut).Options.yQ(iSig,:), ax(rRight, iSig))
      
   end % iSig = 1 : nSignal
   
   
   print(figureHandle,[localFigurePath, subjectID, '_', sessionID, '_', signalModel, '_', output(hOut).Name, '_CDF'],'-dpdf', '-r300')
   
end % hOut

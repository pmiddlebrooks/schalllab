function [solution,LL,AIC,BIC,CDF,trialDataToFit, correct, incorrect, t] = fit_LBA_ccm(subjectID, sessionID, freeOrFixFlag, InitialParamStruct, plotFlag)


% rand('seed',5150);
% randn('seed',5150);
% fitFn = 'chi2';
% fitFn = 'LL';


% CONSTANTS
ALL_FIXED_FLAG = 0;
FREE_BTWN_TARGET_FLAG = 2;
FREE_FLAG = 1;

R_TARG_FLAG = 1;
L_TARG_FLAG = 2;
TARG_FLAG   = 1;
DIST_FLAG   = 2;
GO_FLAG     = 1;
STOP_FLAG   = 2;



MINIMIZE = 1;
collapseDifficulty = 0;
FIT_STOPS = 0;
if FIT_STOPS
   nGoStop = 2;
else
   nGoStop = 1;
end

SCREEN_OR_SAVE = 'save';


% Load the data
switch sessionID
   case 'collapse'
      % Check to see if there's a local usable version of collapsed data
      [dataFile, localDataPath, localDataFile] = data_file_path(subjectID, [sessionID]);
      if exist(localDataFile, 'file') == 2
         load(localDataFile);
         disp('local')
      else
         disp('construct')
         % Need to create a separate file with sessionArray data in
         % it specific to each subject
         sessionArray = ...
            {'bp041n03', ...
            'bp042n02', ...
            'bp043n02', ...
            'bp044n02', ...
            'bp045n02', ...
            'bp046n02', ...
            'bp047n02', ...
            'bp048n02', ...
            'bp049n02', ...
            'bp050n02', ...
            'bp051n02'};
         nSession = length(sessionArray);
         collapsedData = [];
         for iSession = 1 : nSession
            iSessionID = sessionArray{iSession};
            [dataFile, localDataPath, localDataFile] = data_file_path(subjectID, iSessionID);
            % If the file hasn't already been copied to a local directory, do it now
            if exist(localDataFile, 'file') ~= 2
               copyfile(dataFile, localDataPath)
            end
            load(localDataFile);
            if sum(strcmp('checkerArray',get(trialData,'VarNames')))
               trialData.checkerArray = [];
            end
            % Get rid of a few memory-intensive variables we don't need
            % for fitting
            trialData.eyeX = [];
            trialData.eyeY = [];
            trialData.eegData = [];
            collapsedData = [collapsedData; trialData];
         end
         trialData = collapsedData;
         % Save a version of this collapsed data set for temporary use if
         % calling this function iteratively (batching, etc)
         save([localDataPath, sessionID, '.mat'], 'trialData', 'sessionData', 'sessionArray', 'subjectID')
      end
   otherwise
      [dataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessionID);
      % If the file hasn't already been copied to a local directory, do it now
      if exist(localDataFile, 'file') ~= 2
         copyfile(dataFile, localDataPath)
      end
      load(localDataFile);
      if ~strcmp(SessionData.taskID, 'ccm')
         fprintf('Not a choice countermanding session, try again\n')
         return
      end
end

nTrial = size(trialData.trialOutcome, 1);


% Convert cells to doubles if necessary
trialData = cell_to_mat(trialData);
% fprintf('Kludged targetproportion 58 to 59: see fit_LBA_ccm.m line 34\n')
% BRAM:
% trialData.targ1CheckerProp(trialData.targ1CheckerProp == .58) = .59;
% MONKEY:
trialData.targ1CheckerProp(trialData.targ1CheckerProp == .42) = .41;





% Establish vectors for conditions:

% Proportions of checkered stimuli
signalStrength = round(trialData.targ1CheckerProp*100);
signalStrengthArray = unique(signalStrength);
if ismember(50, signalStrengthArray)
   
   signalStrengthArray = [signalStrengthArray(1:find(signalStrengthArray==50));signalStrengthArray(find(signalStrengthArray==50):end)];
end
if ismember(50, signalStrengthArray)
signalStrengthArray(signalStrengthArray == 50) = [50;50];
end
nSignalStrength = length(signalStrengthArray);
% Include option to bin the symmetric levels of difficulty, so all easy,
% harder, hardest, etc. trials will be collapsed regardless of which
% was the correct target location.
if collapseDifficulty
   %     signalStrength = round(abs(signalStrength - 50));
   signalStrength = abs(signalStrength - 50);
   signalStrengthArray = sort(unique(signalStrength));
   %     signalStrengthArray = sort(abs(signalStrengthArray(signalStrengthArray <= 50) - 50));
   %     signalStrengthArray = round(signalStrengthArray*100) / 100
   % unique(signalStrength)
   nSignalStrength = length(signalStrengthArray);
   %     signalStrengthArray = [2 5];
   %     nSignalStrength = length(signalStrengthArray);
end




% SSDs: Need to do a little SSD value adjusting, due to ms difference and 1-frame
% differences in SSD values
ssdArray = unique(trialData.stopSignalOn - trialData.responseCueOn);
ssdArray(isnan(ssdArray)) = [];
dSSD = diff(ssdArray);
ssdArray(dSSD == 1) = ssdArray(dSSD == 1) + 1;
ssdArray = unique(ssdArray);
bS = [ssdArray(1); diff(ssdArray)];
ssdArray(bS < 18) = [];
nSSD = length(ssdArray);



%SET PARAMETERS TO FIXED OR FREE, 1 FOR Fixed, fr FOR NUMBER OF CONDITIONS
%FREE TO VARY ACROSS FOR THAT PARAMETER
% FORMAT:   A  b  v T0  %%% s assumed to always be fixed
nCondition = nSignalStrength;
if nargin < 5
   plotFlag = 1;
end
if nargin < 4
   InitialParamStruct.A = 200;
   InitialParamStruct.b = 300;
   InitialParamStruct.v = .65; %found this to work.  Some models will perform worse than all constrained if this is too high.  I.e., it will get stuck in a terrible local minima.
   InitialParamStruct.T0 = 300;
   InitialParamStruct.s = .3;
end
if nargin < 3;   %set all parameter freedoms as desired
   %     A_freeOrFix = FREE_BTWN_TARGET_FLAG;
   %     b_freeOrFix = FREE_BTWN_TARGET_FLAG;
   %     v_freeOrFix = FREE_FLAG;
   %     t0_freeOrFix = FREE_BTWN_TARGET_FLAG;
   %     s_freeOrFix = ALL_FIXED_FLAG;
   
   A_freeOrFix = FREE_FLAG;
   b_freeOrFix = FREE_FLAG;
   v_freeOrFix = FREE_FLAG;
   t0_freeOrFix = FREE_FLAG;
   s_freeOrFix = ALL_FIXED_FLAG;
   
   freeOrFixFlag = [A_freeOrFix, b_freeOrFix, v_freeOrFix, t0_freeOrFix, s_freeOrFix];
end
% For free parameters, create as many paramters as there are conditions to
% vary. For fixed parameters, create only one.
freeOrFix(freeOrFixFlag == FREE_FLAG) = nCondition;
freeOrFix(freeOrFixFlag == FREE_BTWN_TARGET_FLAG) = 2;  % One fixed value for each target
freeOrFix(freeOrFixFlag == ALL_FIXED_FLAG) = 1;
A(1:nGoStop, 1:freeOrFix(1)) = InitialParamStruct.A;
b(1:nGoStop, 1:freeOrFix(2)) = InitialParamStruct.b;
v(1:nGoStop, 1:freeOrFix(3)) = InitialParamStruct.v;
T0(1:nGoStop, 1:freeOrFix(4)) = InitialParamStruct.T0;
s(1:nGoStop, 1:freeOrFix(5)) = InitialParamStruct.s;





% ****************************************************
%       GATHER AND PREPARE DATA TO FIT
% ****************************************************
% Get default trial selection options
selectOpt = ccm_trial_selection;
selectOpt.rightCheckerPct = 'collapse';

selectOpt.ssd = 'none';
selectOpt.outcome     = {'goCorrectTarget', 'targetHoldAbort'};
goTargetTrial     = ccm_trial_selection(trialData, selectOpt);
selectOpt.outcome     = {'goCorrectDistractor', 'distractorHoldAbort'};
goDistractorTrial = ccm_trial_selection(trialData, selectOpt);

selectOpt.ssd = 'collapse';
selectOpt.outcome       = {'stopIncorrectTarget', 'targetHoldAbort'};
stopTargetTrial     = ccm_trial_selection(trialData, selectOpt);
selectOpt.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort'};
stopDistractorTrial = ccm_trial_selection(trialData, selectOpt);

if FIT_STOPS
   % validTrial = [easyTargetTrial; easyDistractorTrial; hardTargetTrial; hardDistractorTrial];
   validTrial = [goTargetTrial; stopTargetTrial; goDistractorTrial; stopDistractorTrial];
else
   validTrial = [goTargetTrial; goDistractorTrial];
end

if collapseDifficulty
   targetNumber = ones(nTrial, 1);
else
   selectOpt.outcome       =  'collapse';
   selectOpt.ssd = 'any';
   selectOpt.targDir = 'right';
   rightTargetTrial     = ccm_trial_selection(trialData, selectOpt);
   selectOpt.targDir = 'left';
   leftTargetTrial     = ccm_trial_selection(trialData, selectOpt);
   targetNumber = nan(nTrial, 1);
   targetNumber(rightTargetTrial) = R_TARG_FLAG;
   targetNumber(leftTargetTrial) = L_TARG_FLAG;
end


targetDistractor = nan(nTrial, 1);
% targetDistractor([hardTargetTrial; easyTargetTrial]) = 1;
% targetDistractor([hardDistractorTrial; easyDistractorTrial]) = 0;
% if FIT_STOPS
targetDistractor([goTargetTrial; stopTargetTrial]) = TARG_FLAG;
targetDistractor([goDistractorTrial; stopDistractorTrial]) = DIST_FLAG;
% else
%     targetDistractor([goTargetTrial]) = TARG_FLAG;
%     targetDistractor([goDistractorTrial]) = DIST_FLAG;
% end
goStop = nan(nTrial, 1);
goStop([goTargetTrial; goDistractorTrial]) = GO_FLAG;
goStop([stopTargetTrial; stopDistractorTrial]) = STOP_FLAG;


responseOnset = trialData.responseOnset;
responseCueOn = trialData.responseCueOn;
RT = responseOnset - responseCueOn;



% if FIT_STOPS
trialDataToFit = dataset(...
   {targetNumber(validTrial),  	'targetNumber'},...
   {targetDistractor(validTrial), 	'targetDistractor'},...
   {goStop(validTrial),            'goStop'},...
   {RT(validTrial),                'RT'},...
   {signalStrength(validTrial),  	'signalStrength'});
% else
%     trialDataToFit = dataset(...
%         {targetNumber(validTrial),  	'targetNumber'},...
%         {targetDistractor(validTrial), 	'targetDistractor'},...
%         {RT(validTrial),                'RT'},...
%         {signalStrength(validTrial),  	'signalStrength'});
% end
%




% ****************************************************
%      DO THE FITTING
% ****************************************************

% Set upper and lower bounds on parameters. Want these kind of loose, so
% parameter space has lots of room to avoid local minima
lb.A(1:nGoStop, 1:freeOrFix(1)) = 0;
lb.b(1:nGoStop, 1:freeOrFix(2)) = 0;
lb.v(1:nGoStop, 1:freeOrFix(3)) = .5; %need this to be set to .5 minimim.  Otherwise, models with all fixed but v tend to quit with terrible fits.  This corrects the problem.
lb.T0(1:nGoStop, 1:freeOrFix(4)) = 30;

ub.A(1:nGoStop, 1:freeOrFix(1)) = 1000;
ub.b(1:nGoStop, 1:freeOrFix(2)) = 1000;
ub.v(1:nGoStop, 1:freeOrFix(3)) = 1;
ub.T0(1:nGoStop, 1:freeOrFix(4)) = 1000;

tic
if MINIMIZE == 1
   param = [A,b,v,T0];
   lower = [lb.A,lb.b,lb.v,lb.T0];
   upper = [ub.A,ub.b,ub.v,ub.T0];
   %options = optimset('MaxIter', 1000000,'MaxFunEvals', 1000000);
   options = optimset('MaxIter', 100000,'MaxFunEvals', 100000,'useparallel','always');
   %     options = optimset('MaxIter', 100,'MaxFunEvals', 100,'useparallel','always');
   [solution minval exitflag output] = fminsearchbnd(@(param) fit_LBA_logLikelihood_ccm(param,trialDataToFit),param,lower,upper,options);
   
   % options = gaoptimset('Generations',1000,'StallGenLimit',1000,...
   %    'MigrationDirection','forward','TolFun',1e-10);
   % %options.PopInitRange = [0 0 0; 100 100 100];
   % %options = gaoptimset(options,'HybridFcn', {  @fminsearch [] });
   
   %options = gaoptimset('PopulationSize',[ones(1,numel(param))*30],'useparallel','always');
   %options = gaoptimset('PopulationSize',[ones(1,numel(param))*30]);
   %[solution minval exitflag output] = ga(@(param) fitLBA_SAT_calcLL(param,trialData),numel(param),[],[],[],[],lower,upper,[],[],options);
   
   A = solution(1:nGoStop, 1:length(A));
   b = solution(1:nGoStop, length(A)+1:length(A)+length(b));
   v = solution(1:nGoStop, length(A)+length(b)+1:length(A)+length(b)+length(v));
   T0 = solution(1:nGoStop, length(A)+length(b)+length(v)+1:length(solution));
   s = repmat(s,nGoStop,max(freeOrFix));
   
else
   disp('NOT MINIMIZING: USING STARTING VALUES TO FIT')
   param = [A,b,v,T0];
   s = repmat(s,nGoStop,max(freeOrFix));
   
   minval = fit_LBA_logLikelihood_ccm(param,trialDataToFit);
   solution = param; %since we are not minimizing, take starting values as solution.
end
%Bayesian Information Criterion
n_free = numel(param);
nTrial = size(trialData.trialOutcomeToFit,1);
LL = -minval; %reset back to positive number, since we MINIMIZEd the negative to maximize the positive
[AIC BIC] = aicbic(LL,n_free,nTrial);

condTargetTrial = cell(nGoStop, nSignalStrength);
condDistractorTrial = cell(nGoStop, nSignalStrength);
CDF = cell(nGoStop, nSignalStrength);
fiftyPctIndex = 1;
for iCond = 1 : nCondition
   iPct = signalStrengthArray(iCond);
   if collapseDifficulty
      targetNumFlag = 1;  % see fit_LBA_ccm: "targetNumber"
   else
      if iPct < 50
         targetNumFlag = L_TARG_FLAG;
      elseif iPct > 50
         targetNumFlag = R_TARG_FLAG;
      elseif ismember([50 50], signalStrengthArray)
         if (iPct == 50 && fiftyPctIndex == 1)
            targetNumFlag = L_TARG_FLAG;
            fiftyPctIndex = fiftyPctIndex + 1;
         else
            targetNumFlag = R_TARG_FLAG;
         end
      end
   end
   condTargetTrial{1, iCond} = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.goStop == GO_FLAG & trialDataToFit.targetNumber == targetNumFlag);
   condDistractorTrial{1, iCond} = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.goStop == GO_FLAG & trialDataToFit.targetNumber == targetNumFlag);
   %return defective CDFs of current dataset
   CDF{1, iCond} = getDefectiveCDF(condTargetTrial{1, iCond}, condDistractorTrial{1, iCond}, trialDataToFit.RT);
   if FIT_STOPS
      condTargetTrial{2, iCond} = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.goStop == STOP_FLAG& trialDataToFit.targetNumber == targetNumFlag);
      condDistractorTrial{2, iCond} = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.goStop == STOP_FLAG& trialDataToFit.targetNumber == targetNumFlag);
      %return defective CDFs of current dataset
      CDF{2, iCond} = getDefectiveCDF(condTargetTrial{2, iCond}, condDistractorTrial{2, iCond}, trialDataToFit.RT);
   end
end





% Covert each parameter vector back to full nCondition length to make
% coding easier below in loop calls to the parameter values. It is assumed
% a parameter is free over all conditions (already has length =
% nCondition), and makes adjustments if that's not the case.
if freeOrFixFlag(1) == FREE_BTWN_TARGET_FLAG
   A = [repmat(A(1), 1, nCondition/2), repmat(A(2), 1, nCondition/2)];
elseif freeOrFixFlag(1) == ALL_FIXED_FLAG
   A = repmat(A, 1, nCondition);
end
if freeOrFixFlag(2) == FREE_BTWN_TARGET_FLAG
   b = [repmat(b(1), 1, nCondition/2), repmat(b(2), 1, nCondition/2)];
elseif freeOrFixFlag(2) == ALL_FIXED_FLAG
   b = repmat(b, 1, nCondition);
end
if freeOrFixFlag(3) == FREE_BTWN_TARGET_FLAG
   v = [repmat(v(1), 1, nCondition/2), repmat(v(2), 1, nCondition/2)];
elseif freeOrFixFlag(3) == ALL_FIXED_FLAG
   v = repmat(v, 1, nCondition);
end
if freeOrFixFlag(4) == FREE_BTWN_TARGET_FLAG
   T0 = [repmat(T0(1), 1, nCondition/2), repmat(T0(2), 1, nCondition/2)];
elseif freeOrFixFlag(4) == ALL_FIXED_FLAG
   T0 = repmat(T0, 1, nCondition);
end









correct = cell(1, nCondition);
incorrect = cell(1, nCondition);
t = cell(1, nCondition);

for jGoStop = 1 : nGoStop
   fiftyPctIndex = 1;
   for iCond = 1 : nCondition
      iPct = signalStrengthArray(iCond);
      %         correct{jGoStop, iCond} = cumsum(linearballisticPDF(1:1000,A(jGoStop, iCond),b(jGoStop, iCond),v(jGoStop, iCond),s(jGoStop, iCond)) .* (1-linearballisticCDF(1:1000,A(jGoStop, iCond),b(jGoStop, iCond),1-v(jGoStop, iCond),s(jGoStop, iCond))));
      %         incorrect{jGoStop, iCond} = cumsum(linearballisticPDF(1:1000,A(jGoStop, iCond),b(jGoStop, iCond),1-v(jGoStop, iCond),s(jGoStop, iCond)) .* (1-linearballisticCDF(1:1000,A(jGoStop, iCond),b(jGoStop, iCond),v(jGoStop, iCond),s(jGoStop, iCond))));
      %         t{jGoStop, iCond} = (1:1000) + T0(jGoStop, iCond);
      %         if collapseDifficulty
      %             targetNumFlag = 1;  % see fit_LBA_ccm: "targetNumber"
      %         else
      %             if iPct < 50
      %                 targetNumFlag = L_TARG_FLAG;
      %             elseif iPct > 50
      %                 targetNumFlag = R_TARG_FLAG;
      %             elseif ismember([50 50], signalStrengthArray)
      %                 if (iPct == 50 && fiftyPctIndex == 1)
      %                     targetNumFlag = L_TARG_FLAG;
      %                     fiftyPctIndex = fiftyPctIndex + 1;
      %                 else
      %                     targetNumFlag = R_TARG_FLAG;
      %                 end
      %             end
      %         end
      correct{jGoStop, iCond} = cumsum(linearballisticPDF(1:1000,A(jGoStop, iCond),b(jGoStop, iCond),v(jGoStop, iCond),s(jGoStop, iCond)) .* (1-linearballisticCDF(1:1000,A(jGoStop, iCond),b(jGoStop, iCond),1-v(jGoStop, iCond),s(jGoStop, iCond))));
      incorrect{jGoStop, iCond} = cumsum(linearballisticPDF(1:1000,A(jGoStop, iCond),b(jGoStop, iCond),1-v(jGoStop, iCond),s(jGoStop, iCond)) .* (1-linearballisticCDF(1:1000,A(jGoStop, iCond),b(jGoStop, iCond),v(jGoStop, iCond),s(jGoStop, iCond))));
      t{jGoStop, iCond} = (1:1000) + T0(jGoStop, iCond);
      
      %         correct{jGoStop, iCond} = cumsum(linearballisticPDF(1:1000,A(jGoStop, targetNumFlag),b(jGoStop, targetNumFlag),v(jGoStop, iCond),s(jGoStop, targetNumFlag)) .* (1-linearballisticCDF(1:1000,A(jGoStop, targetNumFlag),b(jGoStop, targetNumFlag),1-v(jGoStop, iCond),s(jGoStop, iCond))));
      %         incorrect{jGoStop, iCond} = cumsum(linearballisticPDF(1:1000,A(jGoStop, targetNumFlag),b(jGoStop, targetNumFlag),1-v(jGoStop, iCond),s(jGoStop, iCond)) .* (1-linearballisticCDF(1:1000,A(jGoStop, targetNumFlag),b(jGoStop, targetNumFlag),v(jGoStop, iCond),s(jGoStop, iCond))));
      %         t{jGoStop, iCond} = (1:1000) + T0(jGoStop, targetNumFlag);
   end
end







% %this is kludgy, but will work for plotting.  Tile parameters as if they were free, but just replicate
% %those that are fixed
if ~all(freeOrFix == 1)
   if freeOrFix(1) == 1; A = repmat(A,nGoStop,fr); end
   if freeOrFix(2) == 1; b = repmat(b,nGoStop,fr); end
   if freeOrFix(3) == 1; v = repmat(v,nGoStop,fr); end
   if freeOrFix(4) == 1; T0 = repmat(T0,nGoStop,fr); end
elseif all(freeOrFix == 1)
end


if plotFlag
   nRow = 2;
   nColumn = nCondition;
   figureHandle = 9043;
   [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = model_figure(nRow, nColumn, figureHandle, SCREEN_OR_SAVE);
   for jGoStop = 1 : nGoStop
      fiftyPctIndex = 1;
      for iCond = 1 : nCondition
         iPct = signalStrengthArray(iCond);
         %             if collapseDifficulty
         %                 targetNumFlag = 1;  % see fit_LBA_ccm: "targetNumber"
         %             else
         %                 if iPct < 50
         %                     targetNumFlag = L_TARG_FLAG;
         %                 elseif iPct > 50
         %                     targetNumFlag = R_TARG_FLAG;
         %                 elseif ismember([50 50], signalStrengthArray)
         %                     if (iPct == 50 && fiftyPctIndex == 1)
         %                         targetNumFlag = L_TARG_FLAG;
         %                         fiftyPctIndex = fiftyPctIndex + 1
         %                     else
         %                         targetNumFlag = R_TARG_FLAG;
         %                     end
         %                 end
         %             end
         ax(jGoStop, iCond) = axes('units', 'centimeters', 'position', [xAxesPosition(jGoStop, iCond) yAxesPosition(jGoStop, iCond) axisWidth axisHeight]);
         hold(ax(jGoStop, iCond), 'on')
         %             ax(jGoStop, iCond) = axes('units', 'centimeters', 'position', [xAxesPosition(jGoStop, iCond) yAxesPosition(jGoStop, iCond) axisWidth axisHeight]);
         %             cla(ax(jGoStop, iCond))
         %             hold(ax(jGoStop, iCond), 'on')
         
         ylim(ax(jGoStop, iCond), [0 1])
         xlim(ax(jGoStop, iCond), [0 1000])
         cla(ax(jGoStop, iCond))
         plot(ax(jGoStop, iCond), t{jGoStop, iCond}, correct{jGoStop, iCond}, 'k', t{jGoStop, iCond}, incorrect{jGoStop, iCond}, 'r');
         plot(ax(jGoStop, iCond), CDF{jGoStop, iCond}.correct(:,1), CDF{jGoStop, iCond}.correct(:,2), 'ok', 'markersize', 3)
         plot(ax(jGoStop, iCond),  CDF{jGoStop, iCond}.err(:,1), CDF{jGoStop, iCond}.err(:,2), 'or', 'markersize', 3)
         %             text(100,.95,['A = ' mat2str(round(A(jGoStop, iCond)*100)/100)])
         %             text(100,.90,['b = ' mat2str(round(b(jGoStop, iCond)*100)/100)])
         %             text(100,.85,['v = ' mat2str(round(v(jGoStop, iCond)*100)/100)])
         %             text(100,.80,['s = ' mat2str(round(s(jGoStop, iCond)*100)/100)])
         %             text(100,.75,['T0 = ' mat2str(round(T0(jGoStop, iCond)*100)/100)])
         text(100, -.35,['A = ' mat2str(round(A(jGoStop, iCond)*100)/100)])
         text(100, -.5,['b = ' mat2str(round(b(jGoStop, iCond)*100)/100)])
         text(100, -.65,['v = ' mat2str(round(v(jGoStop, iCond)*100)/100)])
         text(100, -.80,['s = ' mat2str(round(s(jGoStop, iCond)*100)/100)])
         text(100, -.95,['T0 = ' mat2str(round(T0(jGoStop, iCond)*100)/100)])
         %             text(100,.75,['T0 = ' mat2str(round(T0(jGoStop, targetNumFlag)*100)/100)])
         %             ttl = sprintf('Signat Pct: %d', round(signalStrengthArray(iCond)));
         %             title(ttl)
         set(ax(jGoStop, iCond), 'XTickLabel', [], 'YTickLabel', [])
         box off
      end
   end
   
   
   %[ax h] = suplabel(['LL = ' mat2str(round(minval*100)/100) ' BIC = ' mat2str(round(BIC*1000)/1000)],'t');
   %     [axs h] = suplabel(['Model: ' mat2str(freeOrFix) 'LL = ' mat2str(LL)],'t');
   %     set(h,'fontsize',12)
   %
   
   print(gcf, ['~/schalllab/tempfigures/',sessionID, '_', num2str(InitialParamStruct.A), '_', num2str(InitialParamStruct.b), '_', num2str(InitialParamStruct.v), '_', num2str(InitialParamStruct.T0), '_', num2str(InitialParamStruct.s), 'LL ', num2str(freeOrFix),'.pdf'], '-dpdf')
   
   
   
   
   %     [predictedRT] = LBA_ccm_latency_probability_function(correct, incorrect, t, signalStrengthArray, FIT_STOPS);
   %
   %     nPlotColumn = 2;
   %     nPlotRow = 3;
   %     goColor = [0 0 0];
   %     modelGoColor = [0 0 1];
   %     stopColor = [1 0 0];
   %     modelStopColor = [1 0 1];
   %     figureHandleLP = 9234;
   %     [axisWidth, axisHeight, xAxesPosition, yAxesPosition]       = standard_figure(nPlotRow, nPlotColumn, 'portrait', figureHandleLP);
   %
   %     figure(figureHandleLP)
   %     % CHRONOMETRIC: OBSERVED AND PREDICTED
   %     if ~collapseDifficulty
   %         iPlotRow = 1;
   %         iPlotColumn = 2;
   %         ax(iPlotRow, iPlotColumn) = axes('units', 'centimeters', 'position', [xAxesPosition(iPlotRow, iPlotColumn) yAxesPosition(iPlotRow, iPlotColumn) axisWidth axisHeight]);
   %         hold(ax(iPlotRow, iPlotColumn), 'on')
   %         [goTargRT, goDistRT, stopTargRT, stopDistRT] = ccm_chronometric(subjectID, sessionID, 0);
   %         %Data
   %         plot(ax(iPlotRow, iPlotColumn), signalStrengthArray, cellfun(@mean, goTargRT), '-o', 'color', goColor, 'linewidth', 1, 'markerfacecolor', goColor, 'markeredgecolor', goColor)
   %         plot(ax(iPlotRow, iPlotColumn), signalStrengthArray, cellfun(@mean, stopTargRT), '-o', 'color', stopColor, 'linewidth', 1, 'markerfacecolor', stopColor, 'markeredgecolor', stopColor)
   %         %Model
   %         plot(ax(iPlotRow, iPlotColumn), signalStrengthArray, predictedRT.correct(GO_FLAG, :), 'd', 'color', 'b', 'markerfacecolor', 'b', 'markeredgecolor', 'b')
   %         plot(ax(iPlotRow, iPlotColumn), signalStrengthArray, predictedRT.correct(STOP_FLAG, :), 'd', 'color', 'g', 'markerfacecolor', 'g', 'markeredgecolor', 'g')
   %
   %     end
   %
   %
   %
   %
   %
   %
   %     % LATENCY PROBABILITY FUNCTIONS: OBSERVED AND PREDICTED
   %
   %     targetProp = nan(2, nCondition);
   %     distractorProp = nan(2, nCondition);
   %     for jgoStop = 1 : nGoStop
   %         for iCond = 1 : nCondition
   %             targetProp(jgoStop, iCond) = length(condTargetTrial{jgoStop, iCond}) / (length(condTargetTrial{jgoStop, iCond}) + length(condDistractorTrial{jgoStop, iCond}));
   %             distractorProp(jgoStop, iCond) = 1 - targetProp(jgoStop, iCond);
   %         end
   %     end
   %
   %
   %     rtGo = [predictedRT.incorrect(1,:), fliplr(predictedRT.correct(1,:))];
   %     propGo = [distractorProp(1,:), fliplr(targetProp(1,:))];
   %     if FIT_STOPS
   %     rtStop = [predictedRT.incorrect(2,:), fliplr(predictedRT.correct(2,:))];
   %     propStop = [distractorProp(2,:), fliplr(targetProp(2,:))];
   %     end
   %
   %     % RT Distribution
   %     if collapseDifficulty
   %         iTarget = 4;
   %     elseif ~collapseDifficulty
   %         iTarget = 3;
   %     end
   %     iPlotRow = 1;
   %     iPlotColumn = 1;
   %
   %     ax(iPlotRow, iPlotColumn) = axes('units', 'centimeters', 'position', [xAxesPosition(iPlotRow, iPlotColumn) yAxesPosition(iPlotRow, iPlotColumn) axisWidth axisHeight]);
   %     hold(ax(iPlotRow, iPlotColumn), 'on')
   %     %             set(ax(iPlotRow, iPlotColumn), 'xlim', [0 1])
   %     %             set(ax(iPlotRow, iPlotColumn), 'ylim', [0 1])
   %     %             if sum(ijGoTargetTrialIndices) > 1
   %     plot(ax(iPlotRow, iPlotColumn), propGo, rtGo, 's', 'markeredgecolor', modelGoColor, 'markerfacecolor', modelGoColor)
   %         if FIT_STOPS
   %     plot(ax(iPlotRow, iPlotColumn), propStop, rtStop, 's', 'markeredgecolor', modelStopColor, 'markerfacecolor', modelStopColor)
   %         end
   %     %         plot(ax(iPlotRow, iPlotColumn), goDistractorProp(iTarget, :), goDistractorRTMean{iTarget}, 'o', 'markeredgecolor', goColor)
   %     %         plot(ax(iPlotRow, iPlotColumn), stopTargetProp(iTarget, :), stopTargetRTMean{iTarget}, 'o', 'markeredgecolor', stopColor, 'markerfacecolor', stopColor)
   %     %         plot(ax(iPlotRow, iPlotColumn), stopDistractorProp(iTarget, :), stopDistractorRTMean{iTarget}, 'o', 'markeredgecolor', stopColor)
   %
   %
   %     [goTargetProp, goDistractorProp, goTargetRT, goDistractorRT, ...
   %         stopTargetProp, stopDistractorProp, stopTargetRT, stopDistractorRT]...
   %         = ccm_latency_probability_function(subjectID, sessionID, 0);
   %
   %     goTargetRTMean{iTarget}         = cellfun(@mean, goTargetRT(iTarget, :));
   %     goDistractorRTMean{iTarget}     = cellfun(@mean, goDistractorRT(iTarget, :));
   %     stopTargetRTMean{iTarget}       = cellfun(@mean, stopTargetRT(iTarget, :));
   %     stopDistractorRTMean{iTarget}   = cellfun(@mean, stopDistractorRT(iTarget, :));
   %
   %     plot(ax(iPlotRow, iPlotColumn), goTargetProp(iTarget, :), goTargetRTMean{iTarget}, 'o', 'markeredgecolor', goColor, 'markerfacecolor', goColor)
   %     plot(ax(iPlotRow, iPlotColumn), goDistractorProp(iTarget, :), goDistractorRTMean{iTarget}, 'o', 'markeredgecolor', goColor)
   %     plot(ax(iPlotRow, iPlotColumn), stopTargetProp(iTarget, :), stopTargetRTMean{iTarget}, 'o', 'markeredgecolor', stopColor, 'markerfacecolor', stopColor)
   %     plot(ax(iPlotRow, iPlotColumn), stopDistractorProp(iTarget, :), stopDistractorRTMean{iTarget}, 'o', 'markeredgecolor', stopColor)
   %
end







disp(['Optimization ran for ' mat2str(round((toc/60)*1000)/1000) ' minutes'])

% delete(localDataFile);

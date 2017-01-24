function [solution,chi2,cdfData, cdfModel, trialDataToFit, correct, incorrect, t] = fit_LBA_chi2_ccm(subjectID, sessionID, freeOrFixFlag, InitialParamStruct, plotFlag)

chi2 = {};

rand('seed',5150);
randn('seed',5150);
% rand('state',sum(100*clock));
% randn('state',sum(100*clock));


% CONSTANTS
ALL_SHARED_FLAG = 0;
ALL_FREE_FLAG = 1;
SHARED_WITHIN_TARGET_FLAG = 2;
FIXED_FLAG = 3;

R_TARG_FLAG = 1;
L_TARG_FLAG = 2;
TARG_FLAG   = 1;
DIST_FLAG   = 2;
GO_FLAG     = 1;
STOP_FLAG   = 2;




MINIMIZE = 1;
sIsParameter = 1;
% Option to bin the symmetric levels of difficulty, so all easy,
% harder, hardest, etc. trials will be collapsed regardless of which
% was the correct target location.
collapseDifficulty = 0;
collapseAll = 0;
FIT_STOPS = 0;
if FIT_STOPS
    nGoStop = 2;
else
    nGoStop = 1;
end

SCREEN_OR_SAVE = 'screen';


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
% trialData.targ1CheckerProp(trialData.targ1CheckerProp == .42) = .41;





% Establish vectors for conditions:

% Proportions of checkered stimuli
signalStrength = round(trialData.targ1CheckerProp*100);
signalStrengthArray = unique(signalStrength);
signalStrengthArray = signalStrengthArray([6])
if ismember(50, signalStrengthArray)
    
    signalStrengthArray = [signalStrengthArray(1:find(signalStrengthArray==50));signalStrengthArray(find(signalStrengthArray==50):end)];
end
nSignalStrength = length(signalStrengthArray);
if collapseDifficulty
    signalStrength = abs(signalStrength - 50);
    signalStrengthArray = sort(unique(signalStrength));
    nSignalStrength = length(signalStrengthArray);
elseif collapseAll
    signalStrength = ones(size(trialData.trialOutcome, 1), 1);
    signalStrengthArray = 1;
    nSignalStrength = 1;
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
    InitialParamStruct.T0 = 50;
    InitialParamStruct.s = .3;
end
if nargin < 3;   %set all parameter freedoms as desired
    %     A_freeOrFix = SHARED_WITHIN_TARGET_FLAG;
    %     b_freeOrFix = SHARED_WITHIN_TARGET_FLAG;
    %     v_freeOrFix = ALL_FREE_FLAG;
    %     t0_freeOrFix = SHARED_WITHIN_TARGET_FLAG;
    %     s_freeOrFix = ALL_SHARED_FLAG;
    
    A_freeOrFix = ALL_FREE_FLAG;
    b_freeOrFix = ALL_FREE_FLAG;
    v_freeOrFix = ALL_FREE_FLAG;
    T0_freeOrFix = ALL_FREE_FLAG;
    if sIsParameter
        s_freeOrFix = ALL_SHARED_FLAG;
    else
        s_freeOrFix = [];
    end
    
    freeOrFixFlag = [A_freeOrFix, b_freeOrFix, v_freeOrFix, T0_freeOrFix, s_freeOrFix];
end
% For free parameters, create as many paramters as there are conditions to
% vary. For fixed parameters, create only one.
freeOrFix = nan(1, length(freeOrFixFlag));
freeOrFix(freeOrFixFlag == ALL_FREE_FLAG) = nCondition;
freeOrFix(freeOrFixFlag == SHARED_WITHIN_TARGET_FLAG) = 2;  % One fixed value for each target
freeOrFix(freeOrFixFlag == ALL_SHARED_FLAG) = 1;
freeOrFix(freeOrFixFlag == FIXED_FLAG) = 1;
A(1:nGoStop, 1:freeOrFix(1)) = InitialParamStruct.A;
b(1:nGoStop, 1:freeOrFix(2)) = InitialParamStruct.b;
v(1:nGoStop, 1:freeOrFix(3)) = InitialParamStruct.v;
T0(1:nGoStop, 1:freeOrFix(4)) = InitialParamStruct.T0;
if sIsParameter
    s(1:nGoStop, 1:freeOrFix(5)) = InitialParamStruct.s;
else
    s = InitialParamStruct.s;
end








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
    validTrial = [goTargetTrial; stopTargetTrial; goDistractorTrial; stopDistractorTrial];
else
    validTrial = [goTargetTrial; goDistractorTrial];
end

if collapseDifficulty
    targetNumber = ones(nTrial, 1);
else
    rightTargetTrial     = ccm_trial_selection(subjectID, sessionID, {'all'}, 'all', 'any', 'right');
    leftTargetTrial     = ccm_trial_selection(subjectID, sessionID, {'all'}, 'all', 'any', 'left');
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
meanRT = nanmean(RT)
stdRT = nanstd(RT)
iqrRT = iqr(RT)




trialDataToFit = dataset(...
    {targetNumber(validTrial),  	'targetNumber'},...
    {targetDistractor(validTrial), 	'targetDistractor'},...
    {goStop(validTrial),            'goStop'},...
    {RT(validTrial),                'RT'},...
    {signalStrength(validTrial),  	'signalStrength'});

% doNotIncludeTrial = trialDataToFit.RT < 140;
% % doNotIncludeTrial = trialDataToFit.RT < meanRT - 1.5*stdRT | trialDataToFit.RT > meanRT + 1.5*stdRT;
% % doNotIncludeTrial = trialDataToFit.RT < meanRT - 1.2*stdRT ;
% trialDataToFit.RT(doNotIncludeTrial)
% % find(doNotIncludeTrial)
% sum(doNotIncludeTrial) / length(RT)
% trialDataToFit(doNotIncludeTrial, :) = [];
%








% ****************************************************
%      DO THE FITTING
% ****************************************************

% Set upper and lower bounds on parameters. Want these kind of loose, so
% parameter space has lots of room to avoid local minima
lb.A(1:nGoStop, 1:freeOrFix(1)) = 1;
lb.b(1:nGoStop, 1:freeOrFix(2)) = 1;
lb.v(1:nGoStop, 1:freeOrFix(3)) = .52; %need this to be set to .5 minimim.  Otherwise, models with all fixed but v tend to quit with terrible fits.  This corrects the problem.
lb.T0(1:nGoStop, 1:freeOrFix(4)) = 40;
if sIsParameter
    lb.s(1:nGoStop, 1:freeOrFix(5)) = .01;
end

ub.A(1:nGoStop, 1:freeOrFix(1)) = 1000;
ub.b(1:nGoStop, 1:freeOrFix(2)) = 1000;
ub.v(1:nGoStop, 1:freeOrFix(3)) = 1;
ub.T0(1:nGoStop, 1:freeOrFix(4)) = 1000;
if sIsParameter
    ub.s(1:nGoStop, 1:freeOrFix(5)) = 1;
end

tic
if MINIMIZE == 1
    if sIsParameter
        param = [A, b, v, T0, s];
        lower = [lb.A, lb.b, lb.v, lb.T0, lb.s];
        upper = [ub.A, ub.b, ub.v, ub.T0, ub.s];
    else
        param = [A, b, v, T0];
        lower = [lb.A, lb.b, lb.v, lb.T0];
        upper = [ub.A, ub.b, ub.v, ub.T0];
    end
    
    options = optimset('MaxIter', 100000,'MaxFunEvals', 100000,'useparallel','always');
    %         options = optimset('MaxIter', 100,'MaxFunEvals', 100,'useparallel','always');
    [solution minval exitflag output] = fminsearchbnd(@(param) fit_LBA_chi2_discrepancy_ccm(param,trialDataToFit),param,lower,upper,options)
    
    % options = gaoptimset('Generations',1000,'StallGenLimit',1000,...
    %    'MigrationDirection','forward','TolFun',1e-10);
    % %options.PopInitRange = [0 0 0; 100 100 100];
    % %options = gaoptimset(options,'HybridFcn', {  @fminsearch [] });
    
    %     options = gaoptimset('PopulationSize',[ones(1,numel(param))*30],'useparallel','always');
    % %     options = gaoptimset('PopulationSize',[ones(1,numel(param))*30]);
    %     [solution minval exitflag output] = ga(@(param) fit_LBA_chi2_discrepancy_ccm(param,trialDataToFit),numel(param),[],[],[],[],lower,upper,[],[],options);
    
    A = solution(1:nGoStop, 1:length(A));
    b = solution(1:nGoStop, length(A)+1:length(A)+length(b));
    v = solution(1:nGoStop, length(A)+length(b)+1:length(A)+length(b)+length(v));
    T0 = solution(1:nGoStop, length(A)+length(b)+length(v)+1:length(A)+length(b)+length(v)+length(T0));
    if sIsParameter
        s = solution(1:nGoStop, length(A)+length(b)+length(v)+length(T0)+1:length(solution));
    else
        s = repmat(s,nGoStop,max(freeOrFix));
    end
    
    
else
    disp('NOT MINIMIZING: USING STARTING VALUES TO FIT')
    if sIsParameter
        param = [A, b, v, T0, s];
    else
        param = [A, b, v, T0];
        s = repmat(s,nGoStop,max(freeOrFix));
    end
    
    minval = fit_LBA_chi2_discrepancy_ccm(param,trialDataToFit);
    solution = param; %since we are not minimizing, take starting values as solution.
end
%Bayesian Information Criterion
n_free = numel(param);
nTrial = size(trialData.trialOutcomeToFit,1);
% [AIC BIC] = aicbic(LL,n_free,nTrial);
chi2 = minval;

condTargetTrial = cell(nGoStop, nSignalStrength);
condDistractorTrial = cell(nGoStop, nSignalStrength);
cdfData = cell(nGoStop, nSignalStrength);
fiftyPctIndex = 1;
for iCond = 1 : nCondition
    iPct = signalStrengthArray(iCond);
    if collapseDifficulty || collapseAll
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
    cdfData{1, iCond} = getDefectiveCDF(condTargetTrial{1, iCond}, condDistractorTrial{1, iCond}, trialDataToFit.RT);
    if FIT_STOPS
        condTargetTrial{2, iCond} = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.goStop == STOP_FLAG& trialDataToFit.targetNumber == targetNumFlag);
        condDistractorTrial{2, iCond} = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.goStop == STOP_FLAG& trialDataToFit.targetNumber == targetNumFlag);
        %return defective CDFs of current dataset
        cdfData{2, iCond} = getDefectiveCDF(condTargetTrial{2, iCond}, condDistractorTrial{2, iCond}, trialDataToFit.RT);
    end
end





% Covert each parameter vector back to full nCondition length to make
% coding easier below in loop calls to the parameter values. It is assumed
% a parameter is free over all conditions (already has length =
% nCondition), and makes adjustments if that's not the case.
if freeOrFixFlag(1) == SHARED_WITHIN_TARGET_FLAG
    A = [repmat(A(1), 1, nCondition/2), repmat(A(2), 1, nCondition/2)];
elseif freeOrFixFlag(1) == ALL_SHARED_FLAG
    A = repmat(A, 1, nCondition);
end
if freeOrFixFlag(2) == SHARED_WITHIN_TARGET_FLAG
    b = [repmat(b(1), 1, nCondition/2), repmat(b(2), 1, nCondition/2)];
elseif freeOrFixFlag(2) == ALL_SHARED_FLAG
    b = repmat(b, 1, nCondition);
end
if freeOrFixFlag(3) == SHARED_WITHIN_TARGET_FLAG
    v = [repmat(v(1), 1, nCondition/2), repmat(v(2), 1, nCondition/2)];
elseif freeOrFixFlag(3) == ALL_SHARED_FLAG
    v = repmat(v, 1, nCondition);
end
if freeOrFixFlag(4) == SHARED_WITHIN_TARGET_FLAG
    T0 = [repmat(T0(1), 1, nCondition/2), repmat(T0(2), 1, nCondition/2)];
elseif freeOrFixFlag(4) == ALL_SHARED_FLAG
    T0 = repmat(T0, 1, nCondition);
end

if sIsParameter
    if freeOrFixFlag(5) == SHARED_WITHIN_TARGET_FLAG
        s = [repmat(s(1), 1, nCondition/2), repmat(s(2), 1, nCondition/2)];
    elseif freeOrFixFlag(5) == ALL_SHARED_FLAG
        s = repmat(s, 1, nCondition);
    end
else
    s = repmat(s, 1, nCondition);
end






nTrialSim = 10000;

cdfModel = cell(nGoStop, nCondition);
modelTargRT = cell(nGoStop, nCondition);
modelDistRT = cell(nGoStop, nCondition);

Parameter_Mat.A = nan(nTrialSim, 1);
Parameter_Mat.b = nan(nTrialSim, 1);
Parameter_Mat.v = nan(nTrialSim, 2);
Parameter_Mat.T0 = nan(nTrialSim, 1);
Parameter_Mat.s = nan(nTrialSim, 1);

incorrect = cell(nGoStop, nCondition);
t = cell(1, nCondition);
correct = cell(nGoStop, nCondition);

propModelRT1 = cell(nGoStop, nCondition);
propModelRT2 = cell(nGoStop, nCondition);

for jGoStop = 1 : nGoStop
    fiftyPctIndex = 1;
    for iCond = 1 : nCondition
        iPct = signalStrengthArray(iCond);
        if collapseDifficulty || collapseAll
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
        
        % Get the number of target (correct) and distractor
        % (incorrect) trials
        ijTrial = trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag;
        ijTrialTarg = ijTrial & trialDataToFit.targetDistractor == TARG_FLAG;
        ijTrialDist = ijTrial & trialDataToFit.targetDistractor == DIST_FLAG;
        nDataTarg = sum(ijTrialTarg);
        nDataDist = sum(ijTrialDist);
        ijPropTarg = nDataTarg / (nDataTarg + nDataDist);
        ijPropDist = 1 - ijPropTarg;
        
        
        
        
        
        % Generate parameter matrices to simulate data
        Parameter_Mat.A = unifrnd(1, max(1, round(A(jGoStop, iCond))), nTrialSim, 1);
        Parameter_Mat.b = ones(nTrialSim, 1) * b(jGoStop, iCond);
        Parameter_Mat.T0 = ones(nTrialSim, 1) * T0(jGoStop, iCond);
        Parameter_Mat.s = ones(nTrialSim, 1) * s(jGoStop, iCond);
        
        
        % CHANGED DRIFT RATE HERE
        %         vTarg = normrnd(v(jGoStop, iCond), s(1), ijNTarg, 1);
        %         % For now, bound drift rates between 0 and 1
        %         vLB = .02;
        %         vUP = .98;
        %         vTarg(vTarg <= vLB) = vLB;
        %         vTarg(vTarg > vUP) = vUP;
        %         vDist = normrnd(v(jGoStop, iCond), s(1), ijNDist, 1);
        %         vDist(vDist <= vLB) = vLB;
        %         vDist(vDist > vUP) = vUP;
        %
        %         Parameter_Mat.v(1:ijNTarg,1) = vTarg;
        %         Parameter_Mat.v(1:ijNTarg,2) = 1-vTarg;
        %         Parameter_Mat.v(ijNTarg + 1 : nTrialSim,1) = 1-vDist;
        %         Parameter_Mat.v(ijNTarg + 1 : nTrialSim,2) = vDist;
        
        vTarg = normrnd(v(jGoStop, iCond), s(jGoStop, iCond), nTrialSim, 1);
        vDist = normrnd(1-v(jGoStop, iCond), s(jGoStop, iCond), nTrialSim, 1);
        
        % For now, bound drift rates between 0 and 1
        %         vLB = .01;
        %         vUP = .99;
        %         vTarg(vTarg <= vLB) = vLB;
        %         vTarg(vTarg > vUP) = vUP;
        %         vDist = 1 - vTarg;
        vTarg(vTarg < 0) = .0001;
        vDist(vDist < 0) = .0001;
        Parameter_Mat.v(:, 1) = vTarg;
        Parameter_Mat.v(:, 2) = vDist;
        
        modelTargRT{jGoStop, iCond} = Parameter_Mat.T0 + (Parameter_Mat.b - Parameter_Mat.A) ./ Parameter_Mat.v(:,1);
        modelDistRT{jGoStop, iCond} = Parameter_Mat.T0 + (Parameter_Mat.b - Parameter_Mat.A) ./ Parameter_Mat.v(:,2);
        ijModelTargTrial = find(modelTargRT{jGoStop, iCond} <= modelDistRT{jGoStop, iCond});
        ijModelDistTrial = find(modelTargRT{jGoStop, iCond} > modelDistRT{jGoStop, iCond});
        pModelTarg = length(ijModelTargTrial) / nTrialSim;
        
        modelRT = min([modelTargRT{jGoStop, iCond}, modelDistRT{jGoStop, iCond}], [], 2);
        
        modelTargRT{jGoStop, iCond} = modelTargRT{jGoStop, iCond}(ijModelTargTrial);
        modelDistRT{jGoStop, iCond} = modelDistRT{jGoStop, iCond}(ijModelDistTrial);
        figure(9333)
        hist(modelRT, 100)
        
        
        pQuantile = .1 : .1 : .9;
        %         pQuantile = [.1 .3 .5 .7 .9];
        %         go1RT = sort(modelTargRT{jGoStop, iCond});
        quantileGo1Model{jGoStop, iCond} = quantile(modelTargRT{jGoStop, iCond}, pQuantile);
        %         go2RT = sort(modelDistRT{jGoStop, iCond});
        quantileGo2Model{jGoStop, iCond} = quantile(modelDistRT{jGoStop, iCond}, pQuantile);
        
        ijTargRTData = trialDataToFit.RT(ijTrialTarg);
        ijDistRTData = trialDataToFit.RT(ijTrialDist);
        quantileGo1Data{jGoStop, iCond} = quantile(ijTargRTData, pQuantile);
        quantileGo2Data{jGoStop, iCond} = quantile(ijDistRTData, pQuantile);
        
        
        
        
        % Target Accumulator
        nObsTargPrev = 0;
        nExpTargPrev = 0;
        nObsTarg = nan(1, length(pQuantile)+1);
        propExp = nan(1, length(pQuantile)+1);
        nExpTarg = nan(1, length(pQuantile)+1);
        for k = 1 : length(pQuantile)
            nObsTarg(k) = sum(ijTargRTData <= quantileGo1Data{jGoStop, iCond}(k)) - nObsTargPrev;
            %             propExp(k) = sum(modelTargRT{jGoStop, iCond} <= quantileGo1Data{jGoStop, iCond}(k)) / nTrialSim;
            propExp(k) = sum(modelTargRT{jGoStop, iCond} <= quantileGo1Data{jGoStop, iCond}(k)) / length(ijModelTargTrial);
            %             propExp(k) = sum(modelRT <= quantileGo1Data{jGoStop, iCond}(k)) / nTrialSim;
            nExpTarg(k) = propExp(k) * sum(ijTrialTarg) - nExpTargPrev;
            nObsTargPrev = nObsTarg(k) + nObsTargPrev;
            nExpTargPrev = nExpTargPrev + nExpTarg(k);
        end
        nObsTarg(k+1) = sum(ijTargRTData > quantileGo1Data{jGoStop, iCond}(k));
        %         propExp(k+1) = sum(modelTargRT{jGoStop, iCond} > quantileGo1Data{jGoStop, iCond}(k)) / nTrialSim;
        propExp(k+1) = sum(modelTargRT{jGoStop, iCond} > quantileGo1Data{jGoStop, iCond}(k)) / length(ijModelTargTrial);
        %         propExp(k+1) = sum(modelRT > quantileGo1Data{jGoStop, iCond}(k)) / nTrialSim;
        nExpTarg(k+1) = propExp(k+1) * sum(ijTrialTarg);
        %         nExpTarg(k+1) = sum(ijTargRTData > quantileGo1Model{jGoStop, iCond}(k));
        
        cdfModel{jGoStop, iCond}.correct = (cumsum(nExpTarg)./nDataTarg) * pModelTarg;
        cdfModel{jGoStop, iCond}.correct
        
        nExpTarg(nExpTarg == 0) = .0001;
        chi2Accum1 = sum((nObsTarg - nExpTarg).^2 ./ nExpTarg);
        disp('TARG')
        disp([nObsTarg;nExpTarg])
        
        
        
        
        % Distractor Accumulator
        nObsDistPrev = 0;
        nExpDistPrev = 0;
        nObsDist = nan(1, length(pQuantile)+1);
        propExp = nan(1, length(pQuantile)+1);
        nExpDist = nan(1, length(pQuantile)+1);
        for k = 1 : length(pQuantile)
            nObsDist(k) = sum(ijDistRTData <= quantileGo2Data{jGoStop, iCond}(k)) - nObsDistPrev;
            %             propExp(k) = sum(modelDistRT{jGoStop, iCond} <= quantileGo2Data{jGoStop, iCond}(k)) / nTrialSim;
            propExp(k) = sum(modelDistRT{jGoStop, iCond} <= quantileGo2Data{jGoStop, iCond}(k)) / length(ijModelDistTrial);
            %             propExp(k) = sum(modelRT <= quantileGo2Data{jGoStop, iCond}(k)) / nTrialSim;
            nExpDist(k) = propExp(k) * sum(ijTrialDist) - nExpDistPrev;
            nObsDistPrev = nObsDist(k) + nObsDistPrev;
            nExpDistPrev = nExpDistPrev + nExpDist(k);
        end
        nObsDist(k+1) = sum(ijDistRTData > quantileGo2Data{jGoStop, iCond}(k));
        %         propExp(k+1) = sum(modelDistRT{jGoStop, iCond} > quantileGo2Data{jGoStop, iCond}(k)) / nTrialSim;
        propExp(k+1) = sum(modelDistRT{jGoStop, iCond} > quantileGo2Data{jGoStop, iCond}(k)) / length(ijModelDistTrial);
        %         propExp(k+1) = sum(modelRT > quantileGo2Data{jGoStop, iCond}(k)) / nTrialSim;
        nExpDist(k+1) = propExp(k+1) * sum(ijTrialDist);
        %         nExpDist(k+1) = sum(ijDistRTData > quantileGo2Model{jGoStop, iCond}(k));
        
        cdfModel{jGoStop, iCond}.err = (cumsum(nExpDist)./nDataDist) * (1-pModelTarg);
        cdfModel{jGoStop, iCond}.err
        
        nExpDist(nExpDist == 0) = .0001;
        chi2Accum2 = sum((nObsDist - nExpDist).^2 ./ nExpDist);
        
        disp('DIST')
        disp([nObsDist;nExpDist])
        
        chi2(jGoStop, iCond) = chi2Accum1 + chi2Accum2;
        
        
        
        % Get defective CDFs for the data
        cdfData{jGoStop, iCond} = getDefectiveCDF(find(ijTrialTarg), find(ijTrialDist), trialDataToFit.RT, 0, pQuantile);
    end
end







% %this is kludgy, but will work for plotting.  Tile parameters as if they were free, but just replicate
% %those that are fixed
if all(freeOrFix == 1) %for ALL FIXED condition
    %     if ~include_med
    A = repmat(A,nGoStop, nCondition);
    b = repmat(b,nGoStop, nCondition);
    v = repmat(v,nGoStop, nCondition);
    T0 = repmat(T0,nGoStop, nCondition);
    if sIsParameter
        s = repmat(s,nGoStop, nCondition);
    end
elseif ~all(freeOrFix == 1)
    if freeOrFix(1) == 1; A = repmat(A,nGoStop,max(freeOrFix)); end
    if freeOrFix(2) == 1; b = repmat(b,nGoStop,max(freeOrFix)); end
    if freeOrFix(3) == 1; v = repmat(v,nGoStop,max(freeOrFix)); end
    if freeOrFix(4) == 1; T0 = repmat(T0,nGoStop,max(freeOrFix)); end
    if sIsParameter
        if freeOrFix(5) == 1; s = repmat(s,nGoStop,max(freeOrFix)); end
    end
end


if plotFlag
    nRow = 2;
    nColumn = max(2, nCondition);
    figureHandle = 9043;
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = model_figure(nRow, nColumn, figureHandle, SCREEN_OR_SAVE);
    for jGoStop = 1 : nGoStop
        for iCond = 1 : nCondition
            
            ax(jGoStop, iCond) = axes('units', 'centimeters', 'position', [xAxesPosition(jGoStop, iCond) yAxesPosition(jGoStop, iCond) axisWidth axisHeight]);
            hold(ax(jGoStop, iCond), 'on')
            %             ax(jGoStop, iCond) = axes('units', 'centimeters', 'position', [xAxesPosition(jGoStop, iCond) yAxesPosition(jGoStop, iCond) axisWidth axisHeight]);
            %             cla(ax(jGoStop, iCond))
            %             hold(ax(jGoStop, iCond), 'on')
            
            ylim(ax(jGoStop, iCond), [0 1])
            xlim(ax(jGoStop, iCond), [0 1000])
            cla(ax(jGoStop, iCond))
            %             plot(ax(jGoStop, iCond), t{jGoStop, iCond}, correct{jGoStop, iCond}, 'k', t{jGoStop, iCond}, incorrect{jGoStop, iCond}, 'r');
            %             plot(ax(jGoStop, iCond), min(modelTargRT{jGoStop, iCond}):max(modelTargRT{jGoStop, iCond}), propModelRT1{jGoStop, iCond}, 'k');
            %             plot(ax(jGoStop, iCond), min(modelDistRT{jGoStop, iCond}):max(modelDistRT{jGoStop, iCond}), propModelRT2{jGoStop, iCond}, 'r');
            %             plot(ax(jGoStop, iCond), cdfModel{jGoStop, iCond}.correct(:,1), cdfModel{jGoStop, iCond}.correct(:,2), '*k', cdfModel{jGoStop, iCond}.err(:,1), cdfModel{jGoStop, iCond}.err(:,2), '*r')
            plot(ax(jGoStop, iCond), cdfData{jGoStop, iCond}.correct(:,1), cdfModel{jGoStop, iCond}.correct(1:end-1), '*-k', cdfData{jGoStop, iCond}.err(:,1), cdfModel{jGoStop, iCond}.err(1:end-1), '*-r')
            plot(ax(jGoStop, iCond), cdfData{jGoStop, iCond}.correct(:,1), cdfData{jGoStop, iCond}.correct(:,2), 'ok', cdfData{jGoStop, iCond}.err(:,1), cdfData{jGoStop, iCond}.err(:,2), 'or')
            text(100, -.35,['A = ' mat2str(round(A(jGoStop, iCond)*100)/100)])
            text(100, -.5,['b = ' mat2str(round(b(jGoStop, iCond)*100)/100)])
            text(100, -.65,['v = ' mat2str(round(v(jGoStop, iCond)*100)/100)])
            text(100, -.95,['T0 = ' mat2str(round(T0(jGoStop, iCond)*100)/100)])
            text(100, -.80,['s = ' mat2str(round(s(jGoStop, iCond)*100)/100)])
            %             ttl = sprintf('Signat Pct: %d', round(signalStrengthArray(iCond)));
            %             title(ttl)
            %             box off
        end
    end
    
    [axs h] = suplabel(['Model: ' mat2str(freeOrFix) '   X2 = ' mat2str(round(minval*100)/100)],'t');
    set(h,'fontsize',12)
    
    if strcmp(SCREEN_OR_SAVE, 'save')
        print(gcf, ['~/schalllab/tempfigures/',sessionID, '_', num2str(InitialParamStruct.A), '_', num2str(InitialParamStruct.b), '_', num2str(InitialParamStruct.v), '_', num2str(InitialParamStruct.T0), '_', num2str(InitialParamStruct.s), 'Chi ', num2str(freeOrFix),'.pdf'], '-dpdf')
    end
    
    
    
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

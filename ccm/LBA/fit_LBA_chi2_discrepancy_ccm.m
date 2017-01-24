function [discrepencyFn] = fit_LBA_chi2_discrepancy_ccm(param, trialDataToFit)


rand('seed',5150);
randn('seed',5150);
normrnd('seed',5150);
% rand('state',sum(100*clock));
% randn('state',sum(100*clock));

plotFlag = 1;

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

FIT_STOPS = evalin('caller','FIT_STOPS');

if FIT_STOPS
    nGoStop = 2;
else
    nGoStop = 1;
end

signalStrengthArray = evalin('caller', 'signalStrengthArray');
collapseDifficulty = evalin('caller', 'collapseDifficulty');
collapseAll = evalin('caller', 'collapseAll');

nSignalStrength = length(signalStrengthArray);
nCondition = nSignalStrength;

freeOrFix = evalin('caller','freeOrFix');
freeOrFixFlag = evalin('caller','freeOrFixFlag');
if length(freeOrFixFlag) == 5
    sIsParameter = 1;
else
    s = evalin('caller','s');
    sIsParameter = 0;
end

% s = evalin('caller','s');
% s = repmat(s,nGoStop, nCondition);

if FIT_STOPS
    %     A(1,:) = param(evalin('caller','1 : size(A, 2)'));
    %     b(1,:) = param(evalin('caller','size(A, 2)*2 + 1 : size(A, 2)*2 + size(b, 2)'));
    %     v(1,:) = param(evalin('caller','size(A, 2)*2 + size(b, 2)*2 + 1 : size(A, 2)*2 + size(b, 2)*2  + size(v, 2)'));
    %     T0(1,:) = param(evalin('caller','size(A, 2)*2 + size(b, 2)*2 + size(v, 2)*2 + 1 : size(A, 2)*2 + size(b, 2)*2 + size(v, 2)*2 + size(T0, 2)'));
    %
    %     A(2,:) = param(evalin('caller','size(A, 2) + 1 : size(A, 2)*2'));
    %     b(2,:) = param(evalin('caller','size(A, 2)*2 + size(b, 2) + 1 : size(A, 2)*2 + size(b, 2)*2'));
    %     v(2,:) = param(evalin('caller','size(A, 2)*2 + size(b, 2)*2 + size(v, 2) + 1 : size(A, 2)*2 + size(b, 2)*2  + size(v, 2)*2'));
    %     T0(2,:) = param(evalin('caller','size(A, 2)*2 + size(b, 2)*2 + size(v, 2)*2 + size(T0, 2) + 1 : size(A, 2)*2 + size(b, 2)*2 + size(v, 2)*2 + size(T0, 2)*2'));
else
    % Covert each parameter vector back to full nCondition length to make
    % coding easier below in loop calls to the parameter values. It is assumed
    % a parameter is free over all conditions (already has length =
    % nCondition), and makes adjustments if that's not the case.
    A = param(evalin('caller','1:length(A(1,:))'));
    if freeOrFixFlag(1) == SHARED_WITHIN_TARGET_FLAG
        A = [repmat(A(1), 1, nCondition/2), repmat(A(2), 1, nCondition/2)];
    elseif freeOrFixFlag(1) == ALL_SHARED_FLAG
        A = repmat(A, 1, nCondition);
    end
    b = param(evalin('caller','length(A(1,:)) + 1 : length(A(1,:))+length(b(1,:))'));
    if freeOrFixFlag(2) == SHARED_WITHIN_TARGET_FLAG
        b = [repmat(b(1), 1, nCondition/2), repmat(b(2), 1, nCondition/2)];
    elseif freeOrFixFlag(2) == ALL_SHARED_FLAG
        b = repmat(b, 1, nCondition);
    end
    v = param(evalin('caller','length(A(1,:)) + length(b(1,:))+1 : length(A(1,:))+length(b(1,:))+length(v(1,:))'));
    if freeOrFixFlag(3) == SHARED_WITHIN_TARGET_FLAG
        v = [repmat(v(1), 1, nCondition/2), repmat(v(2), 1, nCondition/2)];
    elseif freeOrFixFlag(3) == ALL_SHARED_FLAG
        v = repmat(v, 1, nCondition);
    end
    T0 = param(evalin('caller','length(A(1,:)) + length(b(1,:))+length(v(1,:))+1 : length(A(1,:)) + length(b(1,:))+length(v(1,:))+length(T0(1,:))'));
    if freeOrFixFlag(4) == SHARED_WITHIN_TARGET_FLAG
        T0 = [repmat(T0(1), 1, nCondition/2), repmat(T0(2), 1, nCondition/2)];
    elseif freeOrFixFlag(4) == ALL_SHARED_FLAG
        T0 = repmat(T0, 1, nCondition);
    end
    
    if sIsParameter
        s = param(evalin('caller','length(A(1,:)) + length(b(1,:))+length(v(1,:))+length(T0(1,:))+1 : length(param(1,:))'));
        if freeOrFixFlag(5) == SHARED_WITHIN_TARGET_FLAG
            s = [repmat(s(1), 1, nCondition/2), repmat(s(2), 1, nCondition/2)];
        elseif freeOrFixFlag(5) == ALL_SHARED_FLAG
            s = repmat(s, 1, nCondition);
        end
    else
        s = repmat(s, 1, nCondition);
    end
    
end





nTrialSim = 10000;

modelTargRT = cell(nGoStop, nCondition);
modelDistRT = cell(nGoStop, nCondition);

Parameter_Mat.A = nan(nTrialSim, 1);
Parameter_Mat.b = nan(nTrialSim, 1);
Parameter_Mat.v = nan(nTrialSim, 2);
Parameter_Mat.T0 = nan(nTrialSim, 1);
Parameter_Mat.s = nan(nTrialSim, 1);


%Parameterized based on what is free/fixed
quantileGo1Data = cell(nGoStop, nCondition);
quantileGo2Data = cell(nGoStop, nCondition);
quantileGo1Model = cell(nGoStop, nCondition);
quantileGo2Model = cell(nGoStop, nCondition);
chi2 = nan(nGoStop, nCondition);
cdfData = cell(nGoStop, nCondition);
cdfModel = cell(nGoStop, nCondition);
% correct = cell(nGoStop, nCondition);
% incorrect = cell(nGoStop, nCondition);
% t = cell(nGoStop, nCondition);

% if FIT_STOPS
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
        Parameter_Mat.s = ones(nTrialSim, 1) * s(iCond);
        
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
        %         vDist(vDist <= vLB) = vLB;
        %         vDist(vDist > vUP) = vUP;
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
        
        
        
        
        % TARGET ACCUMULATOR
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
        
        nExpTarg(nExpTarg == 0) = .0001;
        chi2Accum1 = sum((nObsTarg - nExpTarg).^2 ./ nExpTarg);
        disp('TARG')
        disp([nObsTarg;nExpTarg])
        
        
        
        
        % DISTRACTOR ACCUMULATOR
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
        
        nExpDist(nExpDist == 0) = .0001;
        chi2Accum2 = sum((nObsDist - nExpDist).^2 ./ nExpDist);
        
        disp('DIST')
        disp([nObsDist;nExpDist])
        
        chi2(jGoStop, iCond) = chi2Accum1 + chi2Accum2;
        
              
        
        % Get defective CDFs for the data
        cdfModel2{jGoStop, iCond} = getDefectiveCDF(ijModelTargTrial, ijModelDistTrial, modelRT, 0, pQuantile);
        cdfData{jGoStop, iCond} = getDefectiveCDF(find(ijTrialTarg), find(ijTrialDist), trialDataToFit.RT, 0, pQuantile);
        
%         cdfData{jGoStop, iCond}.correct(:,2)
%         cdfData{jGoStop, iCond}.err(:,2)
%         cdfModel{jGoStop, iCond}.correct
%         cdfModel{jGoStop, iCond}.err
    end
end

% Parameter_Mat.s(1:nTrialSim) = s(1);


discrepencyFn = sum(chi2(:));
% trialDataToFit(1:10, :)


%
% [AIC BIC] = aicbic(-LL,n_free,nTrial); %be sure to negate LL because we were minimizing the negative to maximize the positive









%this is kludgy, but will work for plotting.  Tile parameters as if they were free, but just replicate
%those that are fixed
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
    screenOrSave = 'screen';
    figureHandle = 923;
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = model_figure(nRow, nColumn, figureHandle, screenOrSave);
    for jGoStop = 1 : nGoStop
        for iCond = 1 : nCondition
            ax(jGoStop, iCond) = axes('units', 'centimeters', 'position', [xAxesPosition(jGoStop, iCond) yAxesPosition(jGoStop, iCond) axisWidth axisHeight]);
            hold(ax(jGoStop, iCond), 'on');
            
            ylim(ax(jGoStop, iCond), [0 1])
            xlim(ax(jGoStop, iCond), [0 1000])
            cla(ax(jGoStop, iCond))
                        plot(ax(jGoStop, iCond), cdfModel2{jGoStop, iCond}.correct(:,1), cdfModel2{jGoStop, iCond}.correct(:,2), '-k', cdfModel2{jGoStop, iCond}.err(:,1), cdfModel2{jGoStop, iCond}.err(:,2), '-r')
            plot(ax(jGoStop, iCond), cdfData{jGoStop, iCond}.correct(:,1), cdfModel{jGoStop, iCond}.correct(1:end-1), '*k', cdfData{jGoStop, iCond}.err(:,1), cdfModel{jGoStop, iCond}.err(1:end-1), '*r')
            plot(ax(jGoStop, iCond), cdfData{jGoStop, iCond}.correct(:,1), cdfData{jGoStop, iCond}.correct(:,2), 'ok', cdfData{jGoStop, iCond}.err(:,1), cdfData{jGoStop, iCond}.err(:,2), 'or')
            %             plot(ax(jGoStop, iCond), quantileGo1Data{jGoStop, iCond}, pQuantile, 'd', 'markerfacecolor', 'k', 'markeredgecolor', 'k');
            %             plot(ax(jGoStop, iCond), quantileGo2Data{jGoStop, iCond}, pQuantile, 'd', 'markerfacecolor', 'r', 'markeredgecolor', 'r');
            %             plot(ax(jGoStop, iCond), quantileGo1Model{jGoStop, iCond}, pQuantile, 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'k')
            %             plot(ax(jGoStop, iCond),quantileGo2Model{jGoStop, iCond}, pQuantile, 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'r')
            text(100,.95,['A = ' mat2str(round(A(jGoStop, iCond)*100)/100)])
            text(100,.90,['b = ' mat2str(round(b(jGoStop, iCond)*100)/100)])
            text(100,.85,['v = ' mat2str(round(v(jGoStop, iCond)*100)/100)])
            text(100,.75,['T0 = ' mat2str(round(T0(jGoStop, iCond)*100)/100)])
            text(100,.80,['s = ' mat2str(round(s(jGoStop, iCond)*100)/100)])
            ttl = sprintf('Signat Pct: %d', round(signalStrengthArray(iCond)*100));
            title(ttl)
            box off
        end
    end
    %     for jGoStop = 1 : nGoStop
    %         for iCond = 1 : nCondition
    %         end
    %     end
    
    
end
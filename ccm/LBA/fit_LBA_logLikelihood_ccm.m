function [LL] = fit_LBA_logLikelihood_ccm(param, trialDataToFit)

plotFlag = 0;

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

FIT_STOPS = evalin('caller','FIT_STOPS');

if FIT_STOPS
    nGoStop = 2;
else
    nGoStop = 1;
end

signalStrengthArray = evalin('caller', 'signalStrengthArray');
collapseDifficulty = evalin('caller', 'collapseDifficulty');

nSignalStrength = length(signalStrengthArray);
nCondition = nSignalStrength;

freeOrFix = evalin('caller','freeOrFix');
freeOrFixFlag = evalin('caller','freeOrFixFlag');


s = evalin('caller','s');
s = repmat(s,nGoStop, nCondition);

if FIT_STOPS
    A(1,:) = param(evalin('caller','1 : size(A, 2)'));
    b(1,:) = param(evalin('caller','size(A, 2)*2 + 1 : size(A, 2)*2 + size(b, 2)'));
    v(1,:) = param(evalin('caller','size(A, 2)*2 + size(b, 2)*2 + 1 : size(A, 2)*2 + size(b, 2)*2  + size(v, 2)'));
    T0(1,:) = param(evalin('caller','size(A, 2)*2 + size(b, 2)*2 + size(v, 2)*2 + 1 : size(A, 2)*2 + size(b, 2)*2 + size(v, 2)*2 + size(T0, 2)'));
    
    A(2,:) = param(evalin('caller','size(A, 2) + 1 : size(A, 2)*2'));
    b(2,:) = param(evalin('caller','size(A, 2)*2 + size(b, 2) + 1 : size(A, 2)*2 + size(b, 2)*2'));
    v(2,:) = param(evalin('caller','size(A, 2)*2 + size(b, 2)*2 + size(v, 2) + 1 : size(A, 2)*2 + size(b, 2)*2  + size(v, 2)*2'));
    T0(2,:) = param(evalin('caller','size(A, 2)*2 + size(b, 2)*2 + size(v, 2)*2 + size(T0, 2) + 1 : size(A, 2)*2 + size(b, 2)*2 + size(v, 2)*2 + size(T0, 2)*2'));
else
% Covert each parameter vector back to full nCondition length to make
% coding easier below in loop calls to the parameter values. It is assumed
% a parameter is free over all conditions (already has length =
% nCondition), and makes adjustments if that's not the case.
    A = param(evalin('caller','1:length(A(1,:))'));
    if freeOrFixFlag(1) == FREE_BTWN_TARGET_FLAG
        A = [repmat(A(1), 1, nCondition/2), repmat(A(2), 1, nCondition/2)];
    elseif freeOrFixFlag(1) == ALL_FIXED_FLAG
        A = repmat(A, 1, nCondition);
    end
    b = param(evalin('caller','length(A(1,:)) + 1:length(A(1,:))+length(b(1,:))'));
    if freeOrFixFlag(2) == FREE_BTWN_TARGET_FLAG
        b = [repmat(b(1), 1, nCondition/2), repmat(b(2), 1, nCondition/2)];
    elseif freeOrFixFlag(2) == ALL_FIXED_FLAG
        b = repmat(b, 1, nCondition);
    end
    v = param(evalin('caller','length(A(1,:)) + length(b(1,:))+1:length(A(1,:))+length(b(1,:))+length(v(1,:))'));
    if freeOrFixFlag(3) == FREE_BTWN_TARGET_FLAG
        v = [repmat(v(1), 1, nCondition/2), repmat(v(2), 1, nCondition/2)];
    elseif freeOrFixFlag(3) == ALL_FIXED_FLAG
        v = repmat(v, 1, nCondition);
    end
    T0 = param(evalin('caller','length(A(1,:)) + length(b(1,:))+length(v(1,:))+1:length(param(1,:))'));
    if freeOrFixFlag(4) == FREE_BTWN_TARGET_FLAG
        T0 = [repmat(T0(1), 1, nCondition/2), repmat(T0(2), 1, nCondition/2)];
    elseif freeOrFixFlag(4) == ALL_FIXED_FLAG
        T0 = repmat(T0, 1, nCondition);
    end
end





n_free = numel(param);
nTrial = size(trialData.trialOutcomeToFit,1);


Parameter_Mat.A = nan(nTrial, 1);
Parameter_Mat.b = nan(nTrial, 1);
Parameter_Mat.v = nan(nTrial, 2);
Parameter_Mat.T0 = nan(nTrial, 1);
Parameter_Mat.s = nan(nTrial, 1);


%Parameterized based on what is free/fixed


% if FIT_STOPS
for jGoStop = 1 : nGoStop
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
        
                Parameter_Mat.A(find(trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),1) = A(jGoStop, iCond);
        
                Parameter_Mat.b(find(trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),1) = b(jGoStop, iCond);
        
                Parameter_Mat.v(find(trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),1) = v(jGoStop, iCond);
                Parameter_Mat.v(find(trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),2) = 1-v(jGoStop, iCond);
                Parameter_Mat.v(find(trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),1) = 1-v(jGoStop, iCond);
                Parameter_Mat.v(find(trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),2) = v(jGoStop, iCond);
        
                Parameter_Mat.T0(find(trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),1) = T0(jGoStop, iCond);
        
        
        
        
%         Parameter_Mat.A(find(trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),1) = A(jGoStop, targetNumFlag);
%         
%         Parameter_Mat.b(find(trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop),1) = b(jGoStop, targetNumFlag);
%         
%         Parameter_Mat.v(find(trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),1) = v(jGoStop, iCond);
%         Parameter_Mat.v(find(trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),2) = 1-v(jGoStop, iCond);
%         Parameter_Mat.v(find(trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),1) = 1-v(jGoStop, iCond);
%         Parameter_Mat.v(find(trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),2) = v(jGoStop, iCond);
%         
%         Parameter_Mat.T0(find(trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop & trialDataToFit.targetNumber == targetNumFlag),1) = T0(jGoStop, targetNumFlag);
%         
        
        
        
    end
end
% else
%     for iCond = 1 : nCondition
%         iPct = signalStrengthArray(iCond);
%         %         x = trialDataToFit.signalStrength == iPct & trialDataToFit.goStop == jGoStop;
%         %         x(1:10)
%         Parameter_Mat.A(find(trialDataToFit.signalStrength == iPct),1) = A(iCond);
%         %         Parameter_Mat.A(1:10)
%
%         Parameter_Mat.b(find(trialDataToFit.signalStrength == iPct),1) = b(iCond);
%
%         Parameter_Mat.v(find(trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.signalStrength == iPct),1) = v(iCond);
%         Parameter_Mat.v(find(trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.signalStrength == iPct),2) = 1-v(iCond);
%         Parameter_Mat.v(find(trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.signalStrength == iPct),1) = v(iCond);
%         Parameter_Mat.v(find(trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.signalStrength == iPct),2) = 1-v(iCond);
%
%         Parameter_Mat.T0(find(trialDataToFit.signalStrength == iPct),1) = T0(iCond);
%     end
% end
Parameter_Mat.s(1:nTrial) = s(1);

trialDataToFit.RT = trialDataToFit.RT - Parameter_Mat.T0;

% trialDataToFit(1:10, :)




fF(:,1) = linearballisticPDF(trialDataToFit.RT,Parameter_Mat.A,Parameter_Mat.b,Parameter_Mat.v(:,1),Parameter_Mat.s);
fF(:,2) = linearballisticCDF(trialDataToFit.RT,Parameter_Mat.A,Parameter_Mat.b,Parameter_Mat.v(:,2),Parameter_Mat.s);

Likelihood = fF(:,1) .* (1-fF(:,2));

%set lower bounds
Likelihood(find(Likelihood < 1e-5)) = 1e-5;

%Note: log(1) = 0, and probabilities are bounded at 1, so a higher
%likelihood will push you closer and closer to 0, which is a larger and
%larger log(likelihood), where larger = closer to or greater than 0.  I.e., a REALLY well fitting model would have
%a likelihood of p = .99, and log(.99) = -.01. However, because log values between 0 and 1 are
%negative, to find the greatest log(likelihood) you have to negate it.
LL = -sum(log(Likelihood));

[AIC BIC] = aicbic(-LL,n_free,nTrial); %be sure to negate LL because we were minimizing the negative to maximize the positive




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
    if FIT_STOPS
        tttGo = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.goStop == GO_FLAG & trialDataToFit.targetNumber == targetNumFlag);
        dddGo = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.goStop == GO_FLAG & trialDataToFit.targetNumber == targetNumFlag);
        %return defective CDFs of current dataset
        CDF{1, iCond} = getDefectiveCDF(tttGo, dddGo, trialDataToFit.RT+Parameter_Mat.T0);
        tttStop = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.goStop == STOP_FLAG & trialDataToFit.targetNumber == targetNumFlag);
        dddStop = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.goStop == STOP_FLAG & trialDataToFit.targetNumber == targetNumFlag);
        %return defective CDFs of current dataset
        CDF{2, iCond} = getDefectiveCDF(tttStop, dddStop, trialDataToFit.RT+Parameter_Mat.T0);
    else
        tttGo = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == TARG_FLAG & trialDataToFit.targetNumber == targetNumFlag);
        dddGo = find(trialDataToFit.signalStrength == iPct & trialDataToFit.targetDistractor == DIST_FLAG & trialDataToFit.targetNumber == targetNumFlag);
        %return defective CDFs of current dataset
        CDF{1, iCond} = getDefectiveCDF(tttGo, dddGo, trialDataToFit.RT+Parameter_Mat.T0);
   end
end





%this is kludgy, but will work for plotting.  Tile parameters as if they were free, but just replicate
%those that are fixed
if all(freeOrFix == 1) %for ALL FIXED condition
    %     if ~include_med
    A = repmat(A,nGoStop, nCondition);
    b = repmat(b,nGoStop, nCondition);
    v = repmat(v,nGoStop, nCondition);
    T0 = repmat(T0,nGoStop, nCondition);
    s = repmat(s,nGoStop, nCondition);
    %     elseif include_med
    %         A = repmat(A,1,3);
    %         b = repmat(b,1,3);
    %         v = repmat(v,1,3);
    %         T0 = repmat(T0,1,3);
    %         s = repmat(s,1,3);
    %     end
elseif ~all(freeOrFix == 1)
    if freeOrFix(1) == 1; A = repmat(A,nGoStop,max(freeOrFix)); end
    if freeOrFix(2) == 1; b = repmat(b,nGoStop,max(freeOrFix)); end
    if freeOrFix(3) == 1; v = repmat(v,nGoStop,max(freeOrFix)); end
    if freeOrFix(4) == 1; T0 = repmat(T0,nGoStop,max(freeOrFix)); end
end




correct = cell(nGoStop, nCondition);
incorrect = cell(nGoStop, nCondition);
t = cell(nGoStop, nCondition);
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



if plotFlag
    %     figure(figureHandle)
    nRow = 2;
    nColumn = nCondition;
    screenOrSave = 'screen';
    figureHandle = 923;
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = model_figure(nRow, nColumn, figureHandle, screenOrSave);
    %     hold all
    %     fon
    for jGoStop = 1 : nGoStop
        fiftyPctIndex = 1;
        for iCond = 1 : nCondition
            ax(jGoStop, iCond) = axes('units', 'centimeters', 'position', [xAxesPosition(jGoStop, iCond) yAxesPosition(jGoStop, iCond) axisWidth axisHeight]);
            hold(ax(jGoStop, iCond), 'on');
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
%                         fiftyPctIndex = fiftyPctIndex + 1;
%                     else
%                         targetNumFlag = R_TARG_FLAG;
%                     end
%                 end
%             end
            ylim(ax(jGoStop, iCond), [0 1])
            xlim(ax(jGoStop, iCond), [0 1000])
            cla(ax(jGoStop, iCond))
            plot(ax(jGoStop, iCond), t{jGoStop, iCond}, correct{jGoStop, iCond}, 'k', t{jGoStop, iCond}, incorrect{jGoStop, iCond}, 'r');
            plot(ax(jGoStop, iCond), CDF{jGoStop, iCond}.correct(:,1), CDF{jGoStop, iCond}.correct(:,2), 'ok', CDF{jGoStop, iCond}.err(:,1), CDF{jGoStop, iCond}.err(:,2), 'or')
            text(100,.95,['A = ' mat2str(round(A(jGoStop, iCond)*100)/100)])
            text(100,.90,['b = ' mat2str(round(b(jGoStop, iCond)*100)/100)])
%             text(100,.95,['A = ' mat2str(round(A(jGoStop, targetNumFlag)*100)/100)])
%             text(100,.90,['b = ' mat2str(round(b(jGoStop, targetNumFlag)*100)/100)])
            text(100,.85,['v = ' mat2str(round(v(jGoStop, iCond)*100)/100)])
            text(100,.80,['s = ' mat2str(round(s(jGoStop, iCond)*100)/100)])
            text(100,.75,['T0 = ' mat2str(round(T0(jGoStop, iCond)*100)/100)])
%             text(100,.75,['T0 = ' mat2str(round(T0(jGoStop, targetNumFlag)*100)/100)])
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
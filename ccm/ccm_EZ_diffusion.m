function data = ccm_EZ_diffusion(subjectID, plotFlag, figureHandle)


%%
PLOT_CHRON = false;
sessionID = '_concat';   % Load the concatenated file for that subject

fprintf('\n\n\n\n')
disp('*******************************************************************************')
disp('EZ diffusion')

% subjectID = 'Human'
% subjectID = 'Xena'
subjectID = 'Broca'
if nargin < 2
    plotFlag = 1;
end
if nargin < 3
    figureHandle = 4895;
end

task = 'ccm';

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
pSignalArray = ExtraVar.pSignalArray;
ssdArray = ExtraVar.ssdArray;



% Flag to determine whether we want to include stop trial analyses for the
% session
DO_STOPS = 1;
MIN_RT = 120;
MAX_RT = 1200;
nSTD    = 3;

if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end



% Truncate RTs
[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTD);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);
allRT(outlierTrial) = [];





signalArray = unique(trialData.targ1CheckerProp);



if plotFlag
    screenOrSave = 'screen';
    figureHandle = 2348;
    nRow = 3;
    goAx = 1;
    stopAx = 2;
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = model_figure(nRow, length(signalArray), figureHandle, screenOrSave);
    for iAx = 1 : length(signalArray)
        ax(goAx, iAx) = axes('units', 'centimeters', 'position', [xAxesPosition(goAx, iAx) yAxesPosition(goAx, iAx) axisWidth axisHeight]);
        cla
        hold(ax(goAx, iAx), 'on')
        ax(stopAx, iAx) = axes('units', 'centimeters', 'position', [xAxesPosition(stopAx, iAx) yAxesPosition(stopAx, iAx) axisWidth axisHeight]);
        cla
        hold(ax(stopAx, iAx), 'on')
    end
    for iAx = [1 3 5]
        ax(3, iAx) = axes('units', 'centimeters', 'position', [xAxesPosition(3, iAx) yAxesPosition(3, iAx) axisWidth*2 axisHeight]);
        cla
        set(ax(3, iAx), 'xtick', 1:7, 'xticklabel', signalArray*100-50)
        xlabel(ax(3, iAx), 'Signal Strength %')
        hold(ax(3, iAx), 'on')
    end
    stopColor = [.3 .3 .3];
    goColor = [0 0 0];
    stopColor = [1 0 0];
    goColor = [0 1 0];
end









% ___________________   Fit the data to EZ Diffusion model   ___________________

s = 0.1;
goTargRT = cell(length(signalArray), 1);
goDistRT = cell(length(signalArray), 1);
goTargMean = nan(length(signalArray), 1);
goTargVar = nan(length(signalArray), 1);
goPCorrect = nan(length(signalArray), 1);

stopTargRT = cell(length(signalArray), 1);
stopDistRT = cell(length(signalArray), 1);
stopTargMean = nan(length(signalArray), 1);
stopTargVar = nan(length(signalArray), 1);
stopPCorrect = nan(length(signalArray), 1);

    
    rt = [];
    goStop = [];
    signalStrength = [];
    accuracy = [];
% Get default trial selection options
selectOpt       = ccm_trial_selection;
for iSignalInd = 1 : length(signalArray)
    iPct = signalArray(iSignalInd) * 100;
   selectOpt.rightCheckerPct = iPct;
    
    % For EZ diffusion, only use correct choice RTs:
    
    % Go trials
         selectOpt.ssd       = 'none';
         
      selectOpt.outcome     = {'goCorrectTarget', 'targetHoldAbort'};
    iGoTargTrial            = ccm_trial_selection(trialData,  selectOpt);
    iGoTargRT               = allRT(iGoTargTrial);
    goTargRT{iSignalInd}    = iGoTargRT;
    
      selectOpt.outcome     = {'goCorrectDistractor', 'distractorHoldAbort'};
    iGoDistTrial            = ccm_trial_selection(trialData,  selectOpt);
    iGoDistRT               = allRT(iGoDistTrial);
    goDistRT{iSignalInd}    = iGoDistRT;
    
    goTargMean(iSignalInd)  = nanmean(iGoTargRT);
    goTargVar(iSignalInd)   = nanvar(iGoTargRT);
    goPCorrect(iSignalInd)  = length(iGoTargRT) / (length(iGoTargRT) + length(iGoDistRT));
    
    
    % Stop trials
         selectOpt.ssd       = 'collapse';
         
            selectOpt.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
    iStopTargTrial          = ccm_trial_selection(trialData,  selectOpt);
    iStopTargRT             = allRT(iStopTargTrial);
    stopTargRT{iSignalInd}  = iStopTargRT;

    selectOpt.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
    iStopDistTrial          = ccm_trial_selection(trialData,  selectOpt);
    iStopDistRT             = allRT(iStopDistTrial);
    stopDistRT{iSignalInd}  = iStopDistRT;
    
    stopTargMean(iSignalInd) = nanmean(iStopTargRT);
    stopTargVar(iSignalInd) = nanvar(iStopTargRT);
    stopPCorrect(iSignalInd) = length(iStopTargRT) / (length(iStopTargRT) + length(iStopDistRT));
    
    
    % Organize data for ANOVA to check assumptions of EZ diffusion
    rt = [rt; iGoTargRT; iGoDistRT; iStopTargRT; iStopDistRT];
    goStop = [goStop; ones(length([iGoTargRT; iGoDistRT]), 1);  zeros(length([iStopTargRT; iStopDistRT]), 1)];
    signalStrength = [signalStrength; iPct * ones(length([iGoTargRT; iGoDistRT; iStopTargRT; iStopDistRT]), 1)];
    accuracy = [accuracy; ones(length(iGoTargRT), 1);  zeros(length(iGoDistRT), 1);  ones(length(iStopTargRT), 1);  zeros(length(iStopDistRT), 1)];
    
    
    
end

% EZ diffusion parameters:
[goV, goA, goTer]       = EZdiffusionfit(goPCorrect, goTargVar, goTargMean, s);
goV = abs(goV);
[stopV, stopA, stopTer] = EZdiffusionfit(stopPCorrect, stopTargVar, stopTargMean, s);
stopV = abs(stopV);




targMean = nanmean(rt(accuracy == 1));
distMean = nanmean(rt(accuracy == 0));









% ___________________   Perform EZ Diffusion model checks   ___________________
%
% Tests assumptions of EZ diffusion model

% 1: Are the RT distributions right-skewed? Uses D'agostino test
alpha = .05;
nNotSkewed = 0;
for iSignalInd = 1 : length(signalArray)
    iPct = signalArray(iSignalInd) * 100;

    goRT = [goTargRT{iSignalInd}; goDistRT{iSignalInd}];
    [X2, P] = DagosPtest(goRT, alpha, false);
    if P >= alpha
        fprintf('Go trial %d%% signal strength is not right-skewed \n', iPct)
        nNotSkewed = nNotSkewed + 1;
    end

    stopRT = [stopTargRT{iSignalInd}; stopDistRT{iSignalInd}];
    [X2, P] = DagosPtest(stopRT, alpha, false);
    if P >= alpha
        fprintf('Stop trial %d%% signal strength is not right-skewed \n', iPct)
        nNotSkewed = nNotSkewed + 1;
    end

end

if nNotSkewed == 0
    disp('EZ test 1 PASSED:    All conditions were right-skewed')
else
    fprintf('EZ test 1 FAILED:    %d of %d conditions were not right-skewed \n', nNotSkewed, length(signalArray)*2)
end
    




% 2: Are error RTs the same as correct RTs?
% fprintf('\n\n *******************  RT ANOVA  *******************  \n')
[p,table,stats] = anovan(rt,{goStop, accuracy, signalStrength}, 'varnames', {'Stop/Go', 'Targ/Dist', 'Signal'}, 'model', 'full', 'display', 'off');
% fprintf('\nStop vs. Go: \t\tp = %.3f\nTarg vs Dist: \tp = %.3f\nSignal Strength: \tp = %.3f\n', p(1), p(2), p(3))
% disp(table)

eta2Sig = table{4,2} / (table{4,2} + table{end-1,2});
eta2Targ = table{3,2} / (table{3,2} + table{end-1,2});
eta2InhTarg = table{5,2} / (table{2,2} + table{end-1,2});

% fprintf('\n\n *********  RT Multicompare  *********  \n')
c = multcompare(stats, 'dimension', 3, 'display', 'off');
% disp(c)
% disp(stats)

alpha = .05;
targDistP = table{3,end};
if targDistP > alpha
    fprintf('EZ test 2 PASSED:    Correct (%.0f ms) and error (%.0f ms) RTs were the same \n', targMean, distMean)
else
    if targMean > distMean
    fprintf('EZ test 2 FAILED:    Correct RTs (%.0f ms) were slower than errors (%.0f ms) \n', targMean, distMean)
    else
    fprintf('EZ test 2 FAILED:    Correct RTs (%.0f ms) were faster than errors (%.0f ms) \n', targMean, distMean)
    end
end






% 2: Are error RTs ivariant as a function of condition (signal strength)?
[p,table,stats] = anovan(rt(accuracy == 0), {signalStrength(accuracy == 0)}, 'varnames', {'Signal'}, 'model', 'full', 'display', 'off');
c = multcompare(stats, 'display', 'off');
% disp(table)
alpha = .05;
distP = table{2,end};
if distP > alpha
    fprintf('EZ test 3 PASSED:    Error RTs invariant as a function of signal strength \n')
else
    fprintf('EZ test 3 FAIELD:    Error RTs varied as a function of signal strength \n')
    display(c)
end










% ___________________   Plot the fit paramters   ___________________
%
% 1st row is decision variable and mean observed RT for go trials, one plot for signal strength level
% 2nd row is decision variable and mean observed RT for stop trials, one plot for signal strength level
% 3rd row is fit parameters, one plot per parameter
if plotFlag
    
    % Plot decision variable through time (avg), one axes per condition:
    for iAx = 1 : length(signalArray)
        set(ax(goAx, iAx), 'Ylim', [0 1.1*max(goA)/2])
        set(ax(stopAx, iAx), 'Ylim', [0 1.1*max(goA)/2])
        set(ax(goAx, iAx), 'Xlim', [0 max(goTargMean)+50])
        set(ax(stopAx, iAx), 'Xlim', [0 max(goTargMean)+50])
        
        x = [0, goTer(iAx), goTer(iAx)+goA(iAx)/2/goV(iAx)];
        y = [0 0 goA(iAx)/2];
        plot(ax(goAx, iAx), x, y, 'color', goColor, 'linewidth', 2);
        plot(ax(goAx, iAx), [0 max(x)], [goA(iAx)/2, goA(iAx)/2], '--', 'color', goColor, 'linewidth', 2)
        plot(ax(goAx, iAx), goTargMean(iAx), goA(iAx)/2, 'o', 'markersize', 8, 'markerfacecolor', goColor, 'markeredgecolor', 'k')
        goTitle = sprintf('%.0f %%', signalArray(iAx)*100-50);
    title(ax(goAx,iAx), goTitle)
        
        x = [0, stopTer(iAx), stopTer(iAx)+stopA(iAx)/2/stopV(iAx)];
        y = [0 0 stopA(iAx)/2];
        plot(ax(stopAx, iAx), x, y, 'color', stopColor, 'linewidth', 2);
        plot(ax(stopAx, iAx), [0 max(x)], [stopA(iAx)/2, stopA(iAx)/2], '--', 'color', stopColor, 'linewidth', 2)
        plot(ax(stopAx, iAx), stopTargMean(iAx), stopA(iAx)/2, 'o', 'markersize', 8, 'markerfacecolor', stopColor, 'markeredgecolor', 'k')
        
    end
    
    
    % Plot  parameter values:
    plot(ax(3,1), 1:length(signalArray), goTer,'o-', 'color', goColor, 'markerfacecolor', goColor, 'markeredgecolor', 'k')
    plot(ax(3,1), 1:length(signalArray), stopTer,'o-', 'color', stopColor, 'markerfacecolor', stopColor, 'markeredgecolor', 'k')
    title(ax(3,1), 'Non-decision Time')
    plot(ax(3,3), 1:length(signalArray), goV,'o-', 'color', goColor, 'markerfacecolor', goColor, 'markeredgecolor', 'k')
    plot(ax(3,3), 1:length(signalArray), stopV,'o-', 'color', stopColor, 'markerfacecolor', stopColor, 'markeredgecolor', 'k')
    title(ax(3,3), 'Drift rate')
    plot(ax(3,5), 1:length(signalArray), goA,'o-', 'color', goColor, 'markerfacecolor', goColor, 'markeredgecolor', 'k')
    plot(ax(3,5), 1:length(signalArray), stopA,'o-', 'color', stopColor, 'markerfacecolor', stopColor, 'markeredgecolor', 'k')
    title(ax(3,5), 'Boundary separation')
    
end  % plotFlag







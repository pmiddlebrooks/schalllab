function data = ccm_trial_history_population(subjectID, plotFlag, figureHandle)

%%
subjectID = 'human'
subjectID = 'xena'
subjectID = 'broca'
sessionSet = 'behavior1';
% if nargin < 2
    plotFlag = 1;
% end
% if nargin < 3
    figureHandle = 4445;
% end

task = 'ccm';

switch subjectID
    case 'human'
        signalStrength = [.35 .42 .46 .5 .54 .58 .65];
    case 'broca'
        signalStrength = [.41 .45 .48 .5 .52 .55 .59];
    case 'xena'
        signalStrength = [.35 .42 .47 .5 .53 .58 .65];
end
[sessionArray, subjectIDArray] = task_session_array(subjectID, task, sessionSet);
% signalStrength = unique([populationData.signalStrengthLeft{1}, populationData.signalStrengthRight{1}])';



nSession = length(sessionArray);


if plotFlag
    % axes names
    axXGoTarg = 1;
    axXGoTarg2 = 2;
    axXGo = 3;
    axXGo2 = 4;
    axGoDistX = 5;
    axGoDistX2 = 6;
    axStopP = 7;
    
    nRow = 3;
    nColumn = 2;
    screenOrSave = 'screen';
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = big_figure(nRow, nColumn, figureHandle, screenOrSave);
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
    
    goTargColor = [0 0 0];
    goDistColor = [.5 .5 .5];
    stopIncorrectColor = [1 0 0];
    stopStopColor = [.7 .4 .5];
end


goTargGoTargRT = {};
goTargRTGoTarg = {};
goTargGoTargGoTargRT = {};
goTargRTGoTargGoTarg = {};
goTargGoDistGoTargRT = {};
goTargRTGoDistGoTarg = {};
goTargStopIncorrectGoTargRT = {};
goTargRTStopIncorrectGoTarg = {};
goTargStopCorrectGoTargRT = {};
goTargRTStopCorrectGoTarg = {};

pGoDistStopCorrectCollapse = [];
pGoTargStopCorrectCollapse = [];
pStopCorrectStopCorrectCollapse = [];
pStopIncorrectStopCorrectCollapse = [];


for iSession = 1 : nSession
    
    
            iData = ccm_trial_history(subjectIDArray{iSession}, sessionArray{iSession}, 0);
    
    goTargGoTargRT = [goTargGoTargRT; iData.goTargGoTargRT];
    goTargRTGoTarg = [goTargRTGoTarg; iData.goTargRTGoTarg];    
    goTargGoTargGoTargRT = [goTargGoTargGoTargRT; iData.goTargGoTargGoTargRT];
    goTargRTGoTargGoTarg = [goTargRTGoTargGoTarg; iData.goTargRTGoTargGoTarg];    
    goTargGoDistGoTargRT = [goTargGoDistGoTargRT; iData.goTargGoDistGoTargRT];
    goTargRTGoDistGoTarg = [goTargRTGoDistGoTarg; iData.goTargRTGoDistGoTarg];
    goTargStopIncorrectGoTargRT = [goTargStopIncorrectGoTargRT; iData.goTargStopIncorrectGoTargRT];
    goTargRTStopIncorrectGoTarg = [goTargRTStopIncorrectGoTarg; iData.goTargRTStopIncorrectGoTarg];
    goTargStopCorrectGoTargRT = [goTargStopCorrectGoTargRT; iData.goTargStopCorrectGoTargRT];
    goTargRTStopCorrectGoTarg = [goTargRTStopCorrectGoTarg; iData.goTargRTStopCorrectGoTarg];
    
    pGoTargStopCorrectCollapse = [pGoTargStopCorrectCollapse; iData.pGoTargStopCorrectCollapse];
    pGoDistStopCorrectCollapse = [pGoDistStopCorrectCollapse; iData.pGoDistStopCorrectCollapse];
    pStopCorrectStopCorrectCollapse = [pStopCorrectStopCorrectCollapse; iData.pStopCorrectStopCorrectCollapse];
    pStopIncorrectStopCorrectCollapse = [pStopIncorrectStopCorrectCollapse; iData.pStopIncorrectStopCorrectCollapse];
    
    
end

goTargGoTargRT = cellfun(@(x) x', goTargGoTargRT, 'uniformOutput', false);
goTargRTGoTarg = cellfun(@(x) x', goTargRTGoTarg, 'uniformOutput', false);
% goTargGoTargGoTargRT = cellfun(@(x) x', goTargGoTargGoTargRT, 'uniformOutput', false);
% goTargRTGoTargGoTarg = cellfun(@(x) x', goTargRTGoTargGoTarg, 'uniformOutput', false);
% goTargGoDistGoTargRT = cellfun(@(x) x', goTargGoDistGoTargRT, 'uniformOutput', false);
% goTargRTGoDistGoTarg = cellfun(@(x) x', goTargRTGoDistGoTarg, 'uniformOutput', false);
% goTargStopIncorrectGoTargRT = cellfun(@(x) x', goTargStopIncorrectGoTargRT, 'uniformOutput', false);
% goTargRTStopIncorrectGoTarg = cellfun(@(x) x', goTargRTStopIncorrectGoTarg, 'uniformOutput', false);
% goTargStopCorrectGoTargRT = cellfun(@(x) x', goTargStopCorrectGoTargRT, 'uniformOutput', false);
% goTargRTStopCorrectGoTarg = cellfun(@(x) x', goTargRTStopCorrectGoTarg, 'uniformOutput', false);



% goTargGoTargRTMean = cellfun(@mean, goTargGoTargRT);
% goTargRTGoTargMean = cellfun(@mean, goTargRTGoTarg);
% goDistGoTargRTMean = cellfun(@mean, goDistGoTargRT);
% goDistRTGoTargMean = cellfun(@mean, goDistRTGoTarg);
%


goTargGoTargDiff1 = cellfun(@(x,y) x - y, goTargGoTargRT, goTargRTGoTarg, 'uniformOutput', false);
goTargGoTargDiff = cellfun(@(x,y) x - y, goTargGoTargGoTargRT, goTargRTGoTargGoTarg, 'uniformOutput', false);
goDistGoTargDiff = cellfun(@(x,y) x - y, goTargGoDistGoTargRT, goTargRTGoDistGoTarg, 'uniformOutput', false);
stopIncorrectGoTargDiff = cellfun(@(x,y) x - y, goTargStopIncorrectGoTargRT, goTargRTStopIncorrectGoTarg, 'uniformOutput', false);
stopStopGoTargDiff = cellfun(@(x,y) x - y, goTargStopCorrectGoTargRT, goTargRTStopCorrectGoTarg, 'uniformOutput', false);



goTargGoTargDiffMean1 = cellfun(@mean, goTargGoTargDiff1);
goTargGoTargDiffMean = cellfun(@mean, goTargGoTargDiff);
goDistGoTargDiffMean = cellfun(@mean, goDistGoTargDiff);
stopIncorrectGoTargDiffMean = cellfun(@mean, stopIncorrectGoTargDiff);
stopStopGoTargDiffMean = cellfun(@mean, stopStopGoTargDiff);



goTargGoTargDiffCollapse1 = cell(nSession, 1);
goTargGoTargDiffCollapse = cell(nSession, 1);
goDistGoTargDiffCollapse = cell(nSession, 1);
stopIncorrectGoTargDiffCollapse = cell(nSession, 1);
stopStopGoTargDiffCollapse = cell(nSession, 1);
for i = 1 : nSession
    goTargGoTargDiffCollapse1{i} = cell2mat(goTargGoTargDiff1(i,:));
    goTargGoTargDiffCollapse{i} = cell2mat(goTargGoTargDiff(i,:));
    goDistGoTargDiffCollapse{i} = cell2mat(goDistGoTargDiff(i,:));
    stopIncorrectGoTargDiffCollapse{i} = cell2mat(stopIncorrectGoTargDiff(i,:));
    stopStopGoTargDiffCollapse{i} = cell2mat(stopStopGoTargDiff(i,:));
end


goTargGoTargDiffMeanCollapse1 = cellfun(@mean, goTargGoTargDiffCollapse1);
goTargGoTargDiffMeanCollapse = cellfun(@mean, goTargGoTargDiffCollapse);
goDistGoTargDiffMeanCollapse = cellfun(@mean, goDistGoTargDiffCollapse);
stopIncorrectGoTargDiffMeanCollapse = cellfun(@mean, stopIncorrectGoTargDiffCollapse);
stopStopGoTargDiffMeanCollapse = cellfun(@mean, stopStopGoTargDiffCollapse);






% ****************************************************************
% ANOVA WITH GOTARG PAIRS

% ANOVA: RT previous trial effects
anovaData = [goTargGoTargDiffMeanCollapse1; goDistGoTargDiffMeanCollapse; stopStopGoTargDiffMeanCollapse; stopIncorrectGoTargDiffMeanCollapse];
groupTrial = [repmat({'goTarg'}, length(goTargGoTargDiffMeanCollapse1), 1); repmat({'goDist'}, length(goDistGoTargDiffMeanCollapse), 1); repmat({'stopC'}, length(stopStopGoTargDiffMeanCollapse), 1); repmat({'stopI'}, length(stopIncorrectGoTargDiffMeanCollapse), 1)];
% groupInh = [repmat({'go'}, length(goTargGoTargDiffMeanCollapse) + length(goDistGoTargDiffMeanCollapse), 1); repmat({'stop'}, length(stopStopGoTargDiffMeanCollapse) + length(stopIncorrectGoTargDiffMeanCollapse), 1)]

fprintf(' ************ GoTarg PAIRS RT ANOVA ************ \n')
[p,table,stats] = anovan(anovaData,{groupTrial}, 'varnames', {'Trial'}, 'display', 'off');
c = multcompare(stats, 'display', 'off');
fprintf('Trial Type: %.3f\n', p(1))
disp(table)
disp(c)


% ****************************************************************
% T-TESTS WITH GOTARG PAIRS

% t-tests: RT previous trial effects
fprintf('\n ************ RT T-tests ************ \n')
[h,p, ci, stats] = ttest2(goTargGoTargDiffMeanCollapse1, goDistGoTargDiffMeanCollapse);
% stats
fprintf('GoTarg PAIRS goDist-goTarg vs goTarg-goTarg t-test: t(%d) = %.2f \tp = %.5f \t CI: %.2f - %.2f\n', stats.df, stats.tstat, p, ci(1), ci(2));
[h,p, ci, stats] = ttest2(stopIncorrectGoTargDiffMeanCollapse, stopStopGoTargDiffMeanCollapse);
% stats
fprintf('GoTarg PAIRS stopCorr-goTarg vs stopIncorr-goTarg t-test: t(%d) = %.2f \tp = %.5f \t CI: %.2f - %.2f\n', stats.df, stats.tstat, p, ci(1), ci(2));



% ****************************************************************
% ANOVA WITH GOTARG TRIPLETS

% ANOVA: RT previous trial effects
anovaData = [goTargGoTargDiffMeanCollapse; goDistGoTargDiffMeanCollapse; stopStopGoTargDiffMeanCollapse; stopIncorrectGoTargDiffMeanCollapse];
groupTrial = [repmat({'goTarg'}, length(goTargGoTargDiffMeanCollapse), 1); repmat({'goDist'}, length(goDistGoTargDiffMeanCollapse), 1); repmat({'stopC'}, length(stopStopGoTargDiffMeanCollapse), 1); repmat({'stopI'}, length(stopIncorrectGoTargDiffMeanCollapse), 1)];
% groupInh = [repmat({'go'}, length(goTargGoTargDiffMeanCollapse) + length(goDistGoTargDiffMeanCollapse), 1); repmat({'stop'}, length(stopStopGoTargDiffMeanCollapse) + length(stopIncorrectGoTargDiffMeanCollapse), 1)]

fprintf('\n\n ************ ALL TRIPLET RT ANOVA ************ \n')
[p,table,stats] = anovan(anovaData,{groupTrial}, 'varnames', {'Trial'}, 'model', 'full', 'display', 'off');
c = multcompare(stats, 'display', 'off');
fprintf('Trial Type: %.3f\n', p(1))
disp(table)
disp(c)


% ****************************************************************
% T-TESTS WITH GOTARG TRIPLETS

% t-tests: RT previous trial effects
fprintf('\n\n\n************ ALL TRIPLET RT T-tests ************ \n')
[h,p, ci, stats] = ttest2(goTargGoTargDiffMeanCollapse, goDistGoTargDiffMeanCollapse);
% stats
fprintf('ALL TRIPLET goDist-goTarg vs goTarg-goTarg t-test: t(%d) = %.2f \tp = %.5f \t CI: %.2f - %.2f\n', stats.df, stats.tstat, p, ci(1), ci(2));
[h,p, ci, stats] = ttest2(stopIncorrectGoTargDiffMeanCollapse, stopStopGoTargDiffMeanCollapse);
% stats
fprintf('ALL TRIPLET stopCorr-goTarg vs stopIncorr-goTarg t-test: t(%d) = %.2f \tp = %.5f \t CI: %.2f - %.2f\n', stats.df, stats.tstat, p, ci(1), ci(2));



% ****************************************************************
% POST-ERROR SLOWING T-TEST
fprintf('\n\n************ POST-ERROR SLOWING T-TEST ************ \n')
[h, p, ci, stats] = ttest(goDistGoTargDiffMeanCollapse);
rTest = sqrt(stats.tstat^2 / (stats.tstat^2 + stats.df))
fprintf('Post-error slowing t-test: t(%d) = %.2f \tp = %.5f \t CI: %.2f - %.2f\n', stats.df, stats.tstat, p, ci(1), ci(2));
fprintf('Post-error RT diff: %.2f \tSD: %.2f\n', mean(goDistGoTargDiffMeanCollapse), std(goDistGoTargDiffMeanCollapse));

% ****************************************************************
% POST-ERROR SLOWING ANOVA ACROSS SIGNAL STRENGTH
anovaData = goTargGoTargDiffMean;
groupData = repmat(1:7, nSession, 1);
fprintf('\n\n ************ POST-ERROR SLOWING ANOVA ACROSS SIGNAL STRENGTH ************ \n')
[p,table,stats] = anovan(anovaData(:),groupData(:), 'varnames', {'Signal Strength'}, 'model', 'full', 'display', 'off');
c = multcompare(stats, 'display', 'off');
fprintf('Sig Strength: p = %.3f\n', p(1))
disp(table)
disp(c)



% ANOVA: probability of canceling a stop trial
anovaData = [pGoTargStopCorrectCollapse; pGoDistStopCorrectCollapse; pStopCorrectStopCorrectCollapse; pStopIncorrectStopCorrectCollapse];
groupTrial = [repmat({'goTarg'}, length(pGoTargStopCorrectCollapse), 1); repmat({'goDist'}, length(pGoDistStopCorrectCollapse), 1); repmat({'stopC'}, length(pStopCorrectStopCorrectCollapse), 1); repmat({'stopI'}, length(pStopIncorrectStopCorrectCollapse), 1)];
% groupInh = [repmat({'go'}, length(pGoTargStopCorrectCollapse) + length(pGoDistStopCorrectCollapse), 1); repmat({'stop'}, length(pStopCorrectStopCorrectCollapse) + length(pStopIncorrectStopCorrectCollapse), 1)];

fprintf('\n\n ************ p(Stop) ANOVA ************ \n')
[p,table,stats] = anovan(anovaData,{groupTrial}, 'varnames', {'Trial'}, 'display', 'off');
c = multcompare(stats, 'display', 'off');
fprintf('Trial Type: %.3f\n', p(1))
disp(table)
disp(c)


% t-tests: RT previous trial effects
fprintf('\n\n ************ p(Stop) T-tests ************ \n')
[h,p, ci, stats] = ttest2(stopStopGoTargDiffMeanCollapse, pGoDistStopCorrectCollapse);
% stats
fprintf('goTarg vs goDist t-test: t(%d) = %.2f \tp = %.5f \t CI: %.2f - %.2f\n', stats.df, stats.tstat, p, ci(1), ci(2));
[h,p, ci, stats] = ttest2(pStopCorrectStopCorrectCollapse, pStopIncorrectStopCorrectCollapse);
% stats
fprintf('stopCorr vs stopIncorr t-test: t(%d) = %.2f \tp = %.5f \t CI: %.2f - %.2f\n', stats.df, stats.tstat, p, ci(1), ci(2));












gt = plot(ax(axXGoTarg), signalStrength, nanmean(goTargGoTargDiffMean1, 1), '.--', 'color', goTargColor, 'linewidth', 2);
gt = plot(ax(axXGoTarg), signalStrength, nanmean(goTargGoTargDiffMean, 1), '.-', 'color', goTargColor, 'linewidth', 2);
gd = plot(ax(axXGoTarg), signalStrength, nanmean(goDistGoTargDiffMean, 1), '.-', 'color', goDistColor, 'linewidth', 2);
si = plot(ax(axXGoTarg), signalStrength, nanmean(stopIncorrectGoTargDiffMean, 1), '.-', 'color', stopIncorrectColor, 'linewidth', 2);
sc = plot(ax(axXGoTarg), signalStrength, nanmean(stopStopGoTargDiffMean, 1), '.-', 'color', stopStopColor, 'linewidth', 2);
set(ax(axXGoTarg), 'xlim', [signalStrength(1) - choicePlotXMargin signalStrength(end) + choicePlotXMargin])
plot(ax(axXGoTarg), signalStrength(1) - .02, nanmean(goTargGoTargDiffMeanCollapse1, 1), '*', 'markeredgecolor', goTargColor, 'markerfacecolor', goTargColor)
plot(ax(axXGoTarg), signalStrength(1) - .02, nanmean(goTargGoTargDiffMeanCollapse, 1), 'd', 'markeredgecolor', goTargColor, 'markerfacecolor', goTargColor)
plot(ax(axXGoTarg), signalStrength(1) - .02, nanmean(goDistGoTargDiffMeanCollapse, 1), 'd', 'markeredgecolor', goDistColor, 'markerfacecolor', goDistColor)
plot(ax(axXGoTarg), signalStrength(1) - .02, nanmean(stopIncorrectGoTargDiffMeanCollapse, 1), 'd', 'markeredgecolor', stopIncorrectColor, 'markerfacecolor', stopIncorrectColor)
plot(ax(axXGoTarg), signalStrength(1) - .02, nanmean(stopStopGoTargDiffMeanCollapse, 1), 'd', 'markeredgecolor', stopStopColor, 'markerfacecolor', stopStopColor)
plot(ax(axXGoTarg), xlim(ax(axXGoTarg)), [0 0], '--k')
legend([gt gd si sc], 'GoTarg-GoTarg', 'GoDist-GoTarg', 'StopCorr-GoTarg', 'StopIncorr-GoTarg')


plot(ax(axXGoTarg2), 0, nanmean(pGoTargStopCorrectCollapse, 1), 'o', 'markeredgecolor', goTargColor, 'markerfacecolor', goTargColor)
plot(ax(axXGoTarg2), 0, nanmean(pGoDistStopCorrectCollapse, 1), 'o', 'markeredgecolor', goDistColor, 'markerfacecolor', goDistColor)
plot(ax(axXGoTarg2), 0, nanmean(pStopCorrectStopCorrectCollapse, 1), 'o', 'markeredgecolor', stopIncorrectColor, 'markerfacecolor', stopIncorrectColor)
plot(ax(axXGoTarg2), 0, nanmean(pStopIncorrectStopCorrectCollapse, 1), 'o', 'markeredgecolor', stopStopColor, 'markerfacecolor', stopStopColor)

% goTargGoTargRTMean = cell2mat(cellfun(@mean, goTargGoTargRT));
% goTargRTGoTargMean = cell2mat(cellfun(@mean, goTargRTGoTarg));
% goDistGoTargRTMean = cell2mat(cellfun(@mean, goDistGoTargRT));
% goDistRTGoTargMean = cell2mat(cellfun(@mean, goDistRTGoTarg));
%
% goTargGoTargDiff = goTargGoTargRTMean - goTargRTGoTargMean;
% goDistGoTargDiff = goDistGoTargRTMean - goDistRTGoTargMean;






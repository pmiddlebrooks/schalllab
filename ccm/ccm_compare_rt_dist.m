
%%
nBin = 60;
minRT = 100;
localDataPath = local_data_path;
printFlag = 0;




%           COLLECT AND PROCESS THE RTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
h = load('~/schalllab/local_data/huAllSaccade.mat');
b = load('~/schalllab/local_data/Broca_concat.mat');
x = load('~/schalllab/local_data/Xena_concat.mat');

h = cell_to_mat(h.trialData);
b = cell_to_mat(b.trialData);
x = cell_to_mat(x.trialData);

% RTs
hRT = h.responseOnset - h.responseCueOn;
bRT = b.responseOnset - b.responseCueOn;
xRT = x.responseOnset - x.responseCueOn;

hRT(hRT < minRT) = nan;
bRT(bRT < minRT) = nan;
xRT(xRT < minRT) = nan;

% subsample for all correct rightward go trials
hRTR = hRT(h.targ1CheckerProp > .5 & strcmp(h.trialOutcome, 'goCorrectTarget'));
bRTR = bRT(b.targ1CheckerProp > .5 & strcmp(b.trialOutcome, 'goCorrectTarget'));
xRTR = xRT(x.targ1CheckerProp > .5 & strcmp(x.trialOutcome, 'goCorrectTarget'));

% subsample for all correct leftward go trials
hRTL = hRT(h.targ1CheckerProp < .5 & strcmp(h.trialOutcome, 'goCorrectTarget'));
bRTL = bRT(b.targ1CheckerProp < .5 & strcmp(b.trialOutcome, 'goCorrectTarget'));
xRTL = xRT(x.targ1CheckerProp < .5 & strcmp(x.trialOutcome, 'goCorrectTarget'));


downSample = false;
if downSample
   % Right trials
   nSubTrialR = length(hRTR);
   bRand = randperm(length(bRTR));
   bRTR = bRTR(bRand(1:nSubTrialR));
   xRand = randperm(length(xRTR));
   xRTR = xRTR(xRand(1:nSubTrialR));
   
   % Left trials
   nSubTrialL = length(hRTL);
   bRand = randperm(length(bRTL));
   bRTL = bRTL(bRand(1:nSubTrialL));
   xRand = randperm(length(xRTL));
   xRTL = xRTL(xRand(1:nSubTrialL));
end







%       DO EX-GAUSSIAN FITS TO OBTAIN PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From Myerson, Robertson, Hale 2007:
%
% "Unfortunately, tau is not an adequate measure of the skew of an
% ex-Gaussian distribution. Tau reflects the absolute length of the
% right-hand tail of a distribution, whereas in statistics, skew is a
% measure of the asymmetry of a distribution (i.e., the length of the
% right-hand tail relative to the length of the left-hand tail). One
% relevant, but potentially troublesome, consequence is that simple
% slowing, in which all RTs are multiplied by a constant, will produce an
% increase in tau but leave skew unchanged. If one wants to measure skew
% using the parameters of an ex-Gaussian distribution and avoid this
% problem, a more appropriate measure of skew is the ratio of tau to sigma
% (Heathcote, Brown, & Mewhort, 2002). The larger this ratio, the more
% skewed the distribution. To the best of our knowledge, however, previous
% studies have not compared the tau/sigma ratio in younger and older
% adults. The only study to actually measure skew as defined in
% mathematical statistics (i.e., as the normalized third central moment of
% a distribution or skewness, which, like the tau/sigma ratio, remains
% constant over changes in scale; Hays, 1994) did not examine asymptotic
% performance (Salthouse, 1993)."

localFigurePath = local_figure_path;

[hPramR hMinvalR hAICR hBICR] = fitModel(hRTR, nanmean(hRTR), nanstd(hRTR), nanstd(hRTR), 'exGauss');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'exGauss-HumanRight'],'-dpdf', '-r300'); end;

[bPramR bMinvalR bAICR bBICR] = fitModel(bRTR, nanmean(bRTR), nanstd(bRTR), nanstd(bRTR), 'exGauss');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'exGauss-Monkey1Right'],'-dpdf', '-r300'); end;

[xPramR xMinvalR xAICR xBICR] = fitModel(xRTR, nanmean(xRTR), nanstd(xRTR), nanstd(xRTR), 'exGauss');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'exGauss-Monkey2Right'],'-dpdf', '-r300'); end;

[hPramL hMinvalL hAICL hBICL] = fitModel(hRTL, nanmean(hRTL), nanstd(hRTL), nanstd(hRTL), 'exGauss');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'exGauss-HumanLeft'],'-dpdf', '-r300'); end;

[bPramL bMinvalL bAICL bBICL] = fitModel(bRTL, nanmean(bRTL), nanstd(bRTL), nanstd(bRTL), 'exGauss');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'exGauss-Monkey1Left'],'-dpdf', '-r300'); end;

[xPramL xMinvalL xAICL xBICL] = fitModel(xRTL, nanmean(xRTL), nanstd(xRTL), nanstd(xRTL), 'exGauss');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'exGauss-Monkey2Left'],'-dpdf', '-r300'); end;


disp('  **************   EX-GAUSSIAN    **************')
fprintf('Human R:\tmu: %.2f\tsigma: %.2f\tTau: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\tTau/Sigma: %.2f\n', hPramR(1), hPramR(2), hPramR(3), hAICR, hBICR, hMinvalR, hPramR(3)/hPramR(2));
fprintf('Broca R:\tmu: %.2f\tsigma: %.2f\tTau: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\tTau/Sigma: %.2f\n', bPramR(1), bPramR(2), bPramR(3), bAICR, bBICR, bMinvalR, bPramR(3)/bPramR(2));
fprintf('Xena  R:\tmu: %.2f\tsigma: %.2f\tTau: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\tTau/Sigma: %.2f\n\n', xPramR(1), xPramR(2), xPramR(3), xAICR, xBICR, xMinvalR, xPramR(3)/xPramR(2));

fprintf('Human L:\tmu: %.2f\tsigma: %.2f\tTau: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\tTau/Sigma: %.2f\n', hPramL(1), hPramL(2), hPramL(3), hAICL, hBICL, hMinvalL, hPramL(3)/hPramL(2));
fprintf('Broca L:\tmu: %.2f\tsigma: %.2f\tTau: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\tTau/Sigma: %.2f\n', bPramL(1), bPramL(2), bPramL(3), bAICL, bBICL, bMinvalL, bPramL(3)/bPramL(2));
fprintf('Xena  L:\tmu: %.2f\tsigma: %.2f\tTau: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\tTau/Sigma: %.2f\n\n', xPramL(1), xPramL(2), xPramL(3), xAICL, xBICL, xMinvalL, xPramL(3)/xPramL(2));


human.exGauss.right.RT = hRTR;
human.exGauss.right.Tau = hPramR(3);
human.exGauss.right.TauSigmaRatio = hPramR(3)/hPramR(2);

monkey1.exGauss.right.RT = bRTR;
monkey1.exGauss.right.Tau = bPramR(3);
monkey1.exGauss.right.TauSigmaRatio = bPramR(3)/bPramR(2);

monkey2.exGauss.right.RT = xRTR;
monkey2.exGauss.right.Tau = xPramR(3);
monkey2.exGauss.right.TauSigmaRatio = xPramR(3)/xPramR(2);

human.exGauss.left.RT = hRTR;
human.exGauss.left.Tau = hPramR(3);
human.exGauss.left.TauSigmaRatio = hPramR(3)/hPramR(2);

monkey1.exGauss.left.RT = bRTR;
monkey1.exGauss.left.Tau = bPramR(3);
monkey1.exGauss.left.TauSigmaRatio = bPramR(3)/bPramR(2);

monkey2.exGauss.left.RT = xRTR;
monkey2.exGauss.left.Tau = xPramR(3);
monkey2.exGauss.left.TauSigmaRatio = xPramR(3)/xPramR(2);








%       DO GAUSSIAN FITS TO OBTAIN PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hPramR hMinvalR hAICR hBICR] = fitModel(hRTR, nanmean(hRTR), nanstd(hRTR), 'norm');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'Gauss-HumanRight'],'-dpdf', '-r300'); end;

[bPramR bMinvalR bAICR bBICR] = fitModel(bRTR, nanmean(bRTR), nanstd(bRTR), 'norm');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'Gauss-Monkey1Right'],'-dpdf', '-r300'); end;

[xPramR xMinvalR xAICR xBICR] = fitModel(xRTR, nanmean(xRTR), nanstd(xRTR), 'norm');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'Gauss-Monkey2Right'],'-dpdf', '-r300'); end;

[hPramL hMinvalL hAICL hBICL] = fitModel(hRTL, nanmean(hRTL), nanstd(hRTL), 'norm');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'Gauss-HumanLeft'],'-dpdf', '-r300'); end;

[bPramL bMinvalL bAICL bBICL] = fitModel(bRTL, nanmean(bRTL), nanstd(bRTL), 'norm');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'Gauss-Monkey1Left'],'-dpdf', '-r300'); end;

[xPramL xMinvalL xAICL xBICL] = fitModel(xRTL, nanmean(xRTL), nanstd(xRTL), 'norm');
set(gcf, 'PaperOrientation','landscape');
set(gcf, 'paperUnits', 'centimeters', 'paperposition', [2 4 26 10])
if printFlag; print([localFigurePath, 'Gauss-Monkey2Left'],'-dpdf', '-r300'); end;


disp('  **************   GAUSSIAN    **************')
fprintf('Human R:\tmu: %.2f\tsigma: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\n', hPramR(1), hPramR(2), hAICR, hBICR, hMinvalR);
fprintf('Broca R:\tmu: %.2f\tsigma: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\n', bPramR(1), bPramR(2), bAICR, bBICR, bMinvalR);
fprintf('Xena  R:\tmu: %.2f\tsigma: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\n\n', xPramR(1), xPramR(2), xAICR, xBICR, xMinvalR);

fprintf('Human L:\tmu: %.2f\tsigma: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\n', hPramL(1), hPramL(2), hAICL, hBICL, hMinvalL);
fprintf('Broca L:\tmu: %.2f\tsigma: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\n', bPramL(1), bPramL(2), bAICL, bBICL, bMinvalL);
fprintf('Xena  L:\tmu: %.2f\tsigma: %.2f\tAIC: %.1f\tBIC: %.1f\tLL: %.1f\n\n', xPramL(1), xPramL(2), xAICL, xBICL, xMinvalL);








% GET THE PHYSICAL DISTRIBUTIONS FOR PLOTTING PURPOSES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Go RT Distribution Right
timeStep = (max(hRTR) - min(hRTR)) / nBin;
hRTBinValuesR = hist(hRTR, nBin);
distributionArea = sum(hRTBinValuesR * timeStep);
hPDFR = hRTBinValuesR / distributionArea;
hBinCentersR = min(hRTR)+timeStep/2 : timeStep : max(hRTR)-timeStep/2;

% Go RT Distribution Right
timeStep = (max(bRTR) - min(bRTR)) / nBin;
bRTBinValuesR = hist(bRTR, nBin);
distributionArea = sum(bRTBinValuesR * timeStep);
bPDFR = bRTBinValuesR / distributionArea;
bBinCentersR = min(bRTR)+timeStep/2 : timeStep : max(bRTR)-timeStep/2;

% Go RT Distribution Right
timeStep = (max(xRTR) - min(xRTR)) / nBin;
xRTBinValuesR = hist(xRTR, nBin);
distributionArea = sum(xRTBinValuesR * timeStep);
xPDFR = xRTBinValuesR / distributionArea;
xBinCentersR = min(xRTR)+timeStep/2 : timeStep : max(xRTR)-timeStep/2;


% Go RT Distribution Left
timeStep = (max(hRTL) - min(hRTL)) / nBin;
hRTBinValuesL = hist(hRTL, nBin);
distributionArea = sum(hRTBinValuesL * timeStep);
hPDFL = hRTBinValuesL / distributionArea;
hBinCentersL = min(hRTL)+timeStep/2 : timeStep : max(hRTL)-timeStep/2;

% Go RT Distribution Left
timeStep = (max(bRTL) - min(bRTL)) / nBin;
bRTBinValuesL = hist(bRTL, nBin);
distributionArea = sum(bRTBinValuesL * timeStep);
bPDFL = bRTBinValuesL / distributionArea;
bBinCentersL = min(bRTL)+timeStep/2 : timeStep : max(bRTL)-timeStep/2;

% Go RT Distribution Left
timeStep = (max(xRTL) - min(xRTL)) / nBin;
xRTBinValuesL = hist(xRTL, nBin);
distributionArea = sum(xRTBinValuesL * timeStep);
xPDFL = xRTBinValuesL / distributionArea;
xBinCentersL = min(xRTL)+timeStep/2 : timeStep : max(xRTL)-timeStep/2;






plotFlag = 0;
if plotFlag
   %                       PLOTTING
   % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   figure(12);
   hold on;
   title('Right')
   
   plot(hBinCentersR, hPDFR, 'k')
   plot(bBinCentersR, bPDFR, 'b')
   plot(xBinCentersR, xPDFR, 'r')
   
   figure(13);
   hold on;
   title('Left')
   
   plot(hBinCentersL, hPDFL, 'k')
   plot(bBinCentersL, bPDFL, 'b')
   plot(xBinCentersL, xPDFL, 'r')
end % if plotFlag


save([localDataPath,'choiceRTs.mat'], 'human', 'monkey1', 'monkey2');



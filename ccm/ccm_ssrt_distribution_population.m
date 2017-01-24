function data = ccm_ssrt_distribution_population(subjectID, sessionSet)

if nargin < 2
   sessionSet = 'neural1';
end
%%

% subjectID = 'broca';
% sessionSet = 'behavior';
% sessionSet = 'neural1';
task = 'ccm';

if iscell(sessionSet)
    sessionArray = sessionSet;
else
[sessionArray, subjectIDArray] = task_session_array(subjectID, task, sessionSet);
end
nSession = length(sessionArray);

switch lower(subjectID)
   case 'human'
      pSignalArray = [.35 .42 .46 .5 .54 .58 .65];
   case 'broca'
      switch sessionSet
         case 'behavior'
            pSignalArray = [.41 .45 .48 .5 .52 .55 .59];
         case 'neural1'
            pSignalArray = [.41 .44 .47 .53 .56 .59];
         case 'neural2'
            pSignalArray = [.42 .44 .46 .54 .56 .58];
      end
   case 'xena'
      pSignalArray = [.35 .42 .47 .5 .53 .58 .65];
end
nSignal = length(pSignalArray);



ssrtIntW = nan(length(sessionArray), nSignal);

optInh              = ccm_inhibition;
optInh.plotFlag     = 0;
optInh.printPlot     = 0;
optInh.collapseTarg     = 1;

% Determine whether there is a 50% condition
if mod(nSignal, 2)
   optInh.include50    = 1;
else
   optInh.include50    = 0;
end



for i = 1 : nSession
   
   
      sessionArray{i}
   iInh               = ccm_inhibition(subjectIDArray{i}, sessionArray{i}, optInh);
   iInh.ssrtIntegration = cellfun(@nanmean, iInh.ssrtIntegration);
   
   ssrtIntW(i,:) = iInh.ssrtIntegrationWeighted';
   
end


% ssrtIntWMean = mean(ssrtIntW, 1);
ssrtIntWMean = median(ssrtIntW, 1);
%%
figureHandle = 61;
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nSignal, 1, 'portrait', figureHandle);
cMap = ccm_colormap(pSignalArray);

for j = 1 : nSignal
   
   ax(j) = axes('units', 'centimeters', 'position', [xAxesPosition(j) yAxesPosition(j) axisWidth axisHeight]);
   hold on;
   set(ax(j), 'xlim', [min(ssrtIntW(:)) max(ssrtIntW(:))])
   set(ax(j), 'ylim', [0 .5])
   
   
   [jY, jX] = hist(ssrtIntW(:,j), nSession);
   bar(jX, jY/nSession, 1, 'facecolor', cMap(j,:))
   
   plot([ssrtIntWMean(j) ssrtIntWMean(j)], ylim)
end


data.ssrtIntW = ssrtIntW;


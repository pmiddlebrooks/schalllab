% function ccm_sam_simulation
%
% Makes plots, etc, comparing observed data (from an observations file) to
% predicted data (in a prd structure created using sam_spec_job.m).
%%
includeStop = false;
[a, b] = ismember('pStopFailureCorr', prd.Properties.VarNames);
if a; includeStop = true; end
iSTOP = 3;

    

[trialData, SessionData, ExtraVar] = load_data(iSubj, [iSubj,'_concat']);
% pSignalArray = ExtraVar.pSignalArray;
% pSignalArray(pSignalArray == .5) = [];
sessionSet = 'behavior';
switch lower(iSubj)
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
pSignalArray(pSignalArray == .5) = [];
nSignal = length(pSignalArray);
ssdArray = ExtraVar.ssdArray;
% trialData.ssd = cell2mat(trialData.stopSignalOn) - cell2mat(trialData.responseCueOn);
% ssdArray = unique(trialData.ssd(~isnan(trialData.ssd)));
minSSDTrial = 20;
% Get rid of SSDs that have less than minSSDTrial trials
for i = 1 : length(ssdArray)
   iSSD = ssdArray(i);
   if sum(trialData.ssd == iSSD) < minSSDTrial
      trialData(trialData.ssd == iSSD, :) = [];
   end
end
ssdArray = unique(trialData.ssd(~isnan(trialData.ssd)));
% trialData.ssd = cell2mat(trialData.stopSignalOn) - cell2mat(trialData.responseCueOn);


% SAM.io.jobDir                     = '/Users/paulmiddlebrooks/matlab/local_data/';
% SAM.io.jobName                    = 'test';
% SAM.io.outDir                     = ['/Users/paulmiddlebrooks/matlab/local_data/sam/',iSubj,'/output/'];

if subSampleSSDFlag
   SAM.io.obsFile                    = ['/Users/paulmiddlebrooks/matlab/local_data/sam/',iSubj,'/obs_',iSubj,'_concat_sam_sub.mat'];
else
   SAM.io.obsFile                    = ['/Users/paulmiddlebrooks/matlab/local_data/sam/',iSubj,'/obs_',iSubj,'_concat_sam.mat'];
end
load(SAM.io.obsFile);


% Number of simulated trials
nSimObs            = length(prd.rtGoCorr{1}) + length(prd.rtGoComm{1});
nSimPrd          = SAM.sim.nSim(2 : end);

pQuantile = .1 : .1 : .9;       % Quantiles for chi2 statistic


% Pre-allocate cell arrays
% ==========================
nCond           = size(obs, 1);

obsQuantiles    = cell(nCond, 1);
obsTargCDF      = cell(nCond, 1);
obsDistCDF      = cell(nCond, 1);
obsTargT        = cell(nCond, 1);
obsDistT        = cell(nCond, 1);
prdQuantiles    = cell(nCond, 1);
prdTargCDF      = cell(nCond, 1);
prdDistCDF      = cell(nCond, 1);
prdTargCDFDef      = cell(nCond, 1);
prdDistCDFDef      = cell(nCond, 1);
prdTargT        = cell(nCond, 1);
prdDistT        = cell(nCond, 1);

rtObsStopCorr       = cell(nCond, 1);
rtObsStopComm       = cell(nCond, 1);
obsRTStop       = cell(nCond, 1);
rtPrdStopCorr       = cell(nCond, 1);
rtPrdStopComm       = cell(nCond, 1);
obsQuantilesStop    = cell(nCond, 1);
obsTargCDFStop      = cell(nCond, 1);
obsDistCDFStop      = cell(nCond, 1);
obsTargTStop        = cell(nCond, 1);
obsDistTStop        = cell(nCond, 1);
prdQuantilesStop    = cell(nCond, 1);
prdTargCDFStop      = cell(nCond, 1);
prdDistCDFStop      = cell(nCond, 1);
prdTargCDFDefStop      = cell(nCond, 1);
prdDistCDFDefStop      = cell(nCond, 1);
prdTargTStop        = cell(nCond, 1);
prdDistTStop        = cell(nCond, 1);


obsQuantileGoCorr   = cell(nCond, 1);
obsQuantileGoComm   = cell(nCond, 1);

chi2GoCorr          = nan(nCond, 1);
chi2GoComm          = nan(nCond, 1);

obsSsdArray        = obs.ssd(1,:);
ssdTimePoints   = obsSsdArray(1) : obsSsdArray(end);
prdSsdArray      = obsSsdArray;
obsInh          = nan(nCond, length(ssdTimePoints));
prdInh          = nan(nCond, length(ssdTimePoints));

rtMax = nan(nCond, 2);
for iCond = 1 : nCond
   
   
   % GO TRIALS
   % ====================================================================
   
   % Get quantiles and full CDF for observed data
   % ====================================================================
   
   % Establish the trials and rt vector for input to getDefectiveCDF
   obsTrialTarg        = 1 : length(obs.rtGoCorr{iCond});
   obsTrialDist        = length(obsTrialTarg) + 1 : length(obsTrialTarg) + length(obs.rtGoComm{iCond});
   obsRT               = [obs.rtGoCorr{iCond}(:); obs.rtGoComm{iCond}(:)];
   
   % Ignorre omission errors for now when plotting
   pGoCorrObs          = length(obsTrialTarg) / length(obsRT);
   pGoCommObs          = 1 - pGoCorrObs;
   
   % The Quantiles
   obsQuantiles{iCond}                     = getDefectiveCDF(obsTrialTarg, obsTrialDist, obsRT);
   
   
   % The Full CDFs
   if ~isempty(obs.rtGoCorr{iCond})
      [obsTargCDF{iCond}, obsTargT{iCond}]    = ecdf(obs.rtGoCorr{iCond}(:));
      % Make them defective like the quantile data
      obsTargCDF{iCond}                       = obsTargCDF{iCond} .* pGoCorrObs;
   end
   
   if ~isempty(obs.rtGoComm{iCond})
      [obsDistCDF{iCond}, obsDistT{iCond}]    = ecdf(obs.rtGoComm{iCond}(:));
      % Make them defective like the quantile data
      obsDistCDF{iCond}                       = obsDistCDF{iCond} .* pGoCommObs;
   end
   
   
   
   
   
   % Get quantiles and full CDF for predicted data
   % ====================================================================
   
   % Establish the trials and rt vector for input to getDefectiveCDF
   prdTrialTarg        = 1 : length(prd.rtGoCorr{iCond});
   prdTrialDist        = length(prdTrialTarg) + 1 : length(prdTrialTarg) + length(prd.rtGoComm{iCond});
   prdRT               = [prd.rtGoCorr{iCond}(:); prd.rtGoComm{iCond}(:)];
   
   % Ignorre omission errors for now when plotting
   pGoCorrPrd          = length(prdTrialTarg) / length(prdRT);
   pGoCommPrd          = 1 - pGoCorrPrd;
   
   % The Quantiles
   prdQuantiles{iCond}                     = getDefectiveCDF(prdTrialTarg, prdTrialDist, prdRT);
   
   % The Full CDFs
   if ~isempty(prd.rtGoCorr{iCond})
      [prdTargCDF{iCond}, prdTargT{iCond}]    = ecdf(prd.rtGoCorr{iCond}(:));
      % Make them defective like the quantile data
      prdTargCDFDef{iCond}   = prdTargCDF{iCond} .* pGoCorrPrd;
   end
   
   if ~isempty(prd.rtGoComm{iCond})
      [prdDistCDF{iCond}, prdDistT{iCond}]    = ecdf(prd.rtGoComm{iCond}(:));
      % Make them defective like the quantile data
      prdDistCDFDef{iCond}   = prdDistCDF{iCond} .* pGoCommPrd;
   end
   
   
   rtMax(iCond, :) = [max(obsRT) max(prdRT)];
   
   
   
   
   
   % CHI 2 CALCULATION
   % ====================================================================
   
   % TARGET ACCUMULATOR
   nObsTargPrev = 0;
   nPrdTargPrev = 0;
   nObsTarg = nan(length(pQuantile)+1, 1);
   nPrdTarg = nan(length(pQuantile)+1, 1);
   obsQuantileGoCorr{iCond} = quantile(obs.rtGoCorr{iCond}(:), pQuantile);
   prdQuantileGoCorr{iCond} = quantile(prd.rtGoCorr{iCond}(:), pQuantile);
   for k = 1 : length(pQuantile)
      nObsTarg(k) = sum(obs.rtGoCorr{iCond}(:) <= obsQuantileGoCorr{iCond}(k)) - nObsTargPrev;
      pPrd        = sum(prd.rtGoCorr{iCond}(:) <= obsQuantileGoCorr{iCond}(k)) / length(prd.rtGoCorr{iCond}(:));
      nPrdTarg(k) = pPrd * length(obs.rtGoCorr{iCond}(:)) - nPrdTargPrev;
      
      nObsTargPrev = nObsTarg(k) + nObsTargPrev;
      nPrdTargPrev = nPrdTarg(k) + nPrdTargPrev;
   end
   nObsTarg(k+1) = sum(obs.rtGoCorr{iCond}(:) > obsQuantileGoCorr{iCond}(k));
   pPrd        = sum(prd.rtGoCorr{iCond}(:) <= obsQuantileGoCorr{iCond}(k)) / length(prd.rtGoCorr{iCond}(:));
   nPrdTarg(k+1) = pPrd * length(obs.rtGoCorr{iCond}(:)) ;
   
   
   nPrdTarg(nPrdTarg == 0) = .0001;
   chi2GoCorr(iCond) = sum((nObsTarg - nPrdTarg).^2 ./ nPrdTarg);
   
   
   % DISTRACTOR ACCUMULATOR
   nObsDistPrev = 0;
   nPrdDistPrev = 0;
   nObsDist = nan(length(pQuantile)+1, 1);
   nPrdDist = nan(length(pQuantile)+1, 1);
   obsQuantileGoComm{iCond} = quantile(obs.rtGoComm{iCond}(:), pQuantile);
   prdQuantileGoComm{iCond} = quantile(prd.rtGoComm{iCond}(:), pQuantile);
   for k = 1 : length(pQuantile)
      nObsDist(k) = sum(obs.rtGoComm{iCond}(:) <= obsQuantileGoComm{iCond}(k)) - nObsDistPrev;
      pPrd        = sum(prd.rtGoComm{iCond}(:) <= obsQuantileGoComm{iCond}(k)) / length(prd.rtGoComm{iCond}(:));
      nPrdDist(k) = pPrd * length(obs.rtGoComm{iCond}(:)) - nPrdDistPrev;
      
      nObsDistPrev = nObsDist(k) + nObsDistPrev;
      nPrdDistPrev = nPrdDist(k) + nPrdDistPrev;
   end
   nObsDist(k+1) = sum(obs.rtGoComm{iCond}(:) > obsQuantileGoComm{iCond}(k));
   pPrd        = sum(prd.rtGoComm{iCond}(:) <= obsQuantileGoComm{iCond}(k)) / length(prd.rtGoComm{iCond}(:));
   nPrdDist(k+1) = pPrd * length(obs.rtGoComm{iCond}(:)) ;
   
   
   nPrdDist(nPrdDist == 0) = .0001;
   chi2GoComm(iCond) = sum((nObsDist - nPrdDist).^2 ./ nPrdDist);
   
   
   
   
   
   
   
   if includeStop
      % Get Inhibiton functions for observed  and predicted data
      % ====================================================================
      
      [fitParameters, lowestSSE] = Weibull(obsSsdArray, obs.pStopFailureCorr(iCond,:), repmat(nSimObs, length(obsSsdArray), 1));
      obsInh(iCond, :) = weibull_curve(fitParameters, ssdTimePoints);
      
      [fitParameters, lowestSSE] = Weibull(prdSsdArray, prd.pStopFailureCorr(iCond,:), nSimPrd);
      prdInh(iCond, :) = weibull_curve(fitParameters, ssdTimePoints);
      
      %         obs.pStopFailureCorr(iCond,:)
      %         [fitParameters, lowestSSE] = Weibull(ssdArray, obs.pStopFailureCorr(iCond,:), repmat(nSim, length(ssdArray), 1));
      %         obsInh(iCond, :) = weibull_curve(fitParameters, ssdTimePoints);
      %
      %         prd.pStopFailureCorr(iCond,:)
      %         [fitParameters, lowestSSE] = Weibull(ssdArray, prd.pStopFailureCorr(iCond,:), repmat(nSim, length(ssdArray), 1));
      %         prdInh(iCond, :) = weibull_curve(fitParameters, ssdTimePoints);
      
      
      
      % STOP TRIALS
      % ====================================================================
      
      % Get quantiles and full CDF for observed data
      % ====================================================================
[m, ssdInd] = max(sum(cellfun(@length, obs.rtStopFailureCorr), 1))
      
      % Establish the trials and rt vector for input to getDefectiveCDF
      rtObsStopCorr{iCond}              = cell2mat(obs.rtStopFailureCorr(iCond,ssdInd));
      rtObsStopComm{iCond}              = cell2mat(obs.rtStopFailureComm(iCond,ssdInd));
%       rtObsStopCorr{iCond}              = cell2mat(obs.rtStopFailureCorr(iCond,:));
%       rtObsStopComm{iCond}              = cell2mat(obs.rtStopFailureComm(iCond,:));
      obsTrialTargStop        = 1 : length(rtObsStopCorr{iCond});
      obsTrialDistStop        = length(obsTrialTargStop) + 1 : length(obsTrialTargStop) + length(rtObsStopComm{iCond});
      obsRTStop{iCond}               = [rtObsStopCorr{iCond}(:); rtObsStopComm{iCond}(:)];
      
      % Ignorre omission errors for now when plotting
      pStopCorrObs          = length(obsTrialTargStop) / length(obsRTStop{iCond});
      pStopCommObs          = 1 - pStopCorrObs;
      
      % The Quantiles
      obsQuantilesStop{iCond}                     = getDefectiveCDF(obsTrialTargStop, obsTrialDistStop, obsRTStop{iCond});
      
      
      % The Full CDFs
      if ~isempty(rtObsStopCorr{iCond})
         [obsTargCDFStop{iCond}, obsTargTStop{iCond}]    = ecdf(rtObsStopCorr{iCond}(:));
         % Make them defective like the quantile data
         obsTargCDFDefStop{iCond}                       = obsTargCDFStop{iCond} .* pStopCorrObs;
      end
      
      if ~isempty(rtObsStopComm{iCond})
         [obsDistCDFStop{iCond}, obsDistTStop{iCond}]    = ecdf(rtObsStopComm{iCond}(:));
         % Make them defective like the quantile data
         obsDistCDFDefStop{iCond}                       = obsDistCDFStop{iCond} .* pStopCommObs;
      end
      
      % Get quantiles and full CDF for predicted data
      % ====================================================================
      
      % Establish the trials and rt vector for input to getDefectiveCDF
      rtPrdStopCorr{iCond}              = cell2mat(prd.rtStopFailureCorr(iCond,:));
      rtPrdStopComm{iCond}              = cell2mat(prd.rtStopFailureComm(iCond,:));
      prdTrialTargStop        = 1 : length(rtPrdStopCorr{iCond});
      prdTrialDistStop        = length(prdTrialTargStop) + 1 : length(prdTrialTargStop) + length(rtPrdStopComm{iCond});
      prdRTStop               = [rtPrdStopCorr{iCond}(:); rtPrdStopComm{iCond}(:)];
      
      % Ignorre omission errors for now when plotting
      pStopCorrPrd          = length(prdTrialTargStop) / length(prdRTStop);
      pStopCommPrd          = 1 - pStopCorrPrd;
      
      % The Quantiles
      prdQuantilesStop{iCond}                     = getDefectiveCDF(prdTrialTargStop, prdTrialDistStop, prdRTStop);
      
      
      % The Full CDFs
      if ~isempty(rtPrdStopCorr{iCond})
         [prdTargCDFStop{iCond}, prdTargTStop{iCond}]    = ecdf(rtPrdStopCorr{iCond}(:));
         % Make them defective like the quantile data
         prdTargCDFDefStop{iCond}                       = prdTargCDFStop{iCond} .* pStopCorrPrd;
      end
      
      if ~isempty(rtPrdStopComm{iCond})
         [prdDistCDFStop{iCond}, prdDistTStop{iCond}]    = ecdf(rtPrdStopComm{iCond}(:));
         % Make them defective like the quantile data
         prdDistCDFDefStop{iCond}                       = prdDistCDFStop{iCond} .* pStopCommPrd;
      end
      
      
      obsQuantileStopFailureCorr{iCond} = quantile(rtObsStopCorr{iCond}, pQuantile);
      prdQuantileStopFailureCorr{iCond} = quantile(rtPrdStopCorr{iCond}, pQuantile);
      obsQuantileStopFailureComm{iCond} = quantile(rtObsStopComm{iCond}, pQuantile);
      prdQuantileStopFailureComm{iCond} = quantile(rtPrdStopComm{iCond}, pQuantile);
      
      
   end % if includeStop
end % for iCond = 1 : nCond
rtMax = max(rtMax(:));






chi2Tot = sum(chi2GoComm + chi2GoComm)













%       PLOTTING: Defective CDFs
% ====================================================================
cMap = ccm_colormap(pSignalArray);
obsColor    = 'k';
prdColor    = 'r';
goColor    = 'k';
stopColor    = 'r';
inhSize     = 8;
rtSize     = 8;
lwCorr      = 2;

switch iSubj
   case 'human'
rtLim       = 900;      
xMax        = rtLim; % rtMax;
xMin        = 350; % rtMax;
rtLimY      = [300 900];
inhMax      = 1000;
ssrtLimY    = [0 400];
   case {'broca'}      
rtLim       = 400;      
xMax        = rtLim; % rtMax;
xMin        = 100; % rtMax;
rtLimY      = [100 400];
inhMax      = 700;
ssrtLimY    = [0 200];
   case {'xena'}      
rtLim       = 450;      
xMax        = rtLim; % rtMax;
xMin        = 200; % rtMax;
rtLimY      = [100 400];
inhMax      = 700;
ssrtLimY    = [0 200];
end

switch SAM.des.choiceMech.type;
   case 'race'
      titleMech = 'Race';
   case 'ffi'
      titleMech = 'Feed Forward Inhibition';
   case 'li'
      titleMech = 'Lateral Inhibition';
end
figureHandle = 55;
nColumn     = nCond;
nRow        = 4;   % for now, plot 2 rows-- need to add more with stop trials
rowTarg     = 1;
rowStop     = 2;
rowInh      = 3;
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
clf
for iCond = 1 : nCond
   
   
   % Set up plot axes
   % ----------------------------------------------------------------
   
   % Go trials: Correct and Error
   ax(rowTarg, iCond) = axes('units', 'centimeters', 'position', [xAxesPosition(rowTarg, iCond) yAxesPosition(rowTarg, iCond) axisWidth axisHeight]);
   cla
   set(ax(rowTarg, iCond), 'ylim', [0 1], 'xlim', [xMin xMax])
   hold(ax(rowTarg, iCond), 'on')
   %             ttl = sprintf('SSD:  %d', iSSD);
   %             title(ttl)
   
   if includeStop
      
      % Stop Trials: Correct and Error
      ax(rowStop, iCond) = axes('units', 'centimeters', 'position', [xAxesPosition(rowStop, iCond) yAxesPosition(rowStop, iCond) axisWidth axisHeight]);
      cla
      set(ax(rowStop, iCond), 'ylim', [0 1], 'xlim', [xMin xMax])
      hold(ax(rowStop, iCond), 'on')
      
      % Inhibition functions
      ax(rowInh, iCond) = axes('units', 'centimeters', 'position', [xAxesPosition(rowInh, iCond) yAxesPosition(rowInh, iCond) axisWidth axisHeight]);
      cla
      set(ax(rowInh, iCond), 'ylim', [0 1], 'xlim', [0 inhMax])
      hold(ax(rowInh, iCond), 'on')
   end
   if iCond == 1
      title(ax(rowTarg, iCond), titleMech)
   end
   if iCond > 1
      set(ax(rowTarg, iCond), 'yticklabel', [])
      if includeStop
         set(ax(rowStop, iCond), 'yticklabel', [])
         set(ax(rowInh, iCond), 'yticklabel', [])
      end
   end
   
   
   
   % Plot the data
   % ----------------------------------------------------------------
   
   % Plot the full CDFs for Target (correct)
   %       plot(ax(rowTarg, iCond), obsTargT{iCond}, obsTargCDFDef{iCond}, 'b')
%    plot(ax(rowTarg, iCond), prdTargT{iCond}, prdTargCDFDef{iCond}, goColor, 'lineWidth', lwCorr)
   
   plot(ax(rowTarg, iCond), obsQuantiles{iCond}.correct(:,1), obsQuantiles{iCond}.correct(:,2), 'o', 'color', goColor, 'markerSize', rtSize, 'lineWidth', lwCorr)
       plot(ax(rowTarg, iCond), prdQuantiles{iCond}.correct(:,1), prdQuantiles{iCond}.correct(:,2), 'color', goColor, 'lineWidth', lwCorr)
   
   
   
   % Plot the full CDFs for Distacotr (errors)
   %     plot(ax(rowTarg, iCond), obsDistT{iCond}, obsDistCDF{iCond}, obsColor)
%    plot(ax(rowTarg, iCond), prdDistT{iCond}, prdDistCDFDef{iCond}, '--', 'color', goColor, 'lineWidth', lwCorr)
   
   plot(ax(rowTarg, iCond), obsQuantiles{iCond}.err(:,1), obsQuantiles{iCond}.err(:,2), 'd', 'color', goColor, 'markerSize', rtSize, 'lineWidth', lwCorr)
       plot(ax(rowTarg, iCond), prdQuantiles{iCond}.err(:,1), prdQuantiles{iCond}.err(:,2), '--', 'color', goColor, 'lineWidth', lwCorr)
   
   
   if includeStop
      
      
      
      % Plot the full CDFs for Target (correct)
      %     plot(ax(rowStop, iCond), obsTargTStop{iCond}, obsTargCDFDefStop{iCond}, obsColor)
%       plot(ax(rowStop, iCond), prdTargTStop{iCond}, prdTargCDFDefStop{iCond}, stopColor, 'lineWidth', lwCorr)
      
      plot(ax(rowStop, iCond), obsQuantilesStop{iCond}.correct(:,1), obsQuantilesStop{iCond}.correct(:,2), 'o', 'color', stopColor, 'markerSize', rtSize, 'lineWidth', lwCorr)
          plot(ax(rowStop, iCond), prdQuantilesStop{iCond}.correct(:,1), prdQuantilesStop{iCond}.correct(:,2), 'color', stopColor, 'lineWidth', lwCorr)
      
      
      
      % Plot the full CDFs for Distacotr (errors)
      %     plot(ax(rowStop, iCond), obsDistTStop{iCond}, obsDistCDFStop{iCond}, obsColor)
%       plot(ax(rowStop, iCond), prdDistTStop{iCond}, prdDistCDFDefStop{iCond}, '--', 'color', stopColor, 'lineWidth', lwCorr)
      
      plot(ax(rowStop, iCond), obsQuantilesStop{iCond}.err(:,1), obsQuantilesStop{iCond}.err(:,2), 'd', 'color', stopColor, 'markerSize', rtSize, 'lineWidth', lwCorr)
          plot(ax(rowStop, iCond), prdQuantilesStop{iCond}.err(:,1), prdQuantilesStop{iCond}.err(:,2), '--', 'color', stopColor, 'lineWidth', lwCorr)
      
      
      
      % Also plot inhibition functions
          plot(ax(rowInh, iCond), obs.ssd(iCond,:), obs.pStopFailureCorr(iCond,:), 'o-', 'color', obsColor, 'markerSize', inhSize, 'lineWidth', lwCorr)
          plot(ax(rowInh, iCond), obs.ssd(iCond,:), prd.pStopFailureCorr(iCond,:), 'o-', 'color', prdColor, 'markerSize', inhSize, 'lineWidth', lwCorr)
      %
      %     % Plot the Inhibition functions
      plot(ax(rowInh, iCond), ssdTimePoints, obsInh(iCond, :), '--b')
      plot(ax(rowInh, iCond), ssdTimePoints, prdInh(iCond, :), '-b')
      %
   end
end % for iCond = 1 : nCond












%%



overlapCDF = true



%       PLOTTING: Full CDFs
% ====================================================================
inhSize     = 10;
rtSize     = 10;
lwCorr      = 2;

switch SAM.des.choiceMech.type;
   case 'race'
      titleMech = 'Race';
   case 'ffi'
      titleMech = 'Feed Forward Inhibition';
   case 'li'
      titleMech = 'Lateral Inhibition';
end
figureHandle = 57;


if overlapCDF
   
   nColumn     = 2;
   nRow        = 3;   % for now, plot 2 rows-- need to add more with stop trials
   rowTargF     = 1;
   rowStopF     = 2;
   rowStopI     = 3;
   [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nRow, nColumn, 'portrait', figureHandle);
   clf
   
   % Set up plot axes
   % ----------------------------------------------------------------
   % Left Go trials: Correct and Error
   ax(rowTargF, 1) = axes('units', 'centimeters', 'position', [xAxesPosition(rowTargF, 1) yAxesPosition(rowTargF, 1) axisWidth axisHeight]);
   cla
   set(ax(rowTargF, 1), 'ylim', [0 1], 'yticklabel', [], 'xlim', [xMin xMax])
   hold(ax(rowTargF, 1), 'on')
   
   % Right Go trials: Correct and Error
   ax(rowTargF, 2) = axes('units', 'centimeters', 'position', [xAxesPosition(rowTargF, 2) yAxesPosition(rowTargF, 2) axisWidth axisHeight]);
   cla
   set(ax(rowTargF, 2), 'ylim', [0 1], 'yticklabel', [], 'xlim', [xMin xMax])
   hold(ax(rowTargF, 2), 'on')
   
   if includeStop
      
      % Left Stop Trials: Correct and Error
      ax(rowStopF, 1) = axes('units', 'centimeters', 'position', [xAxesPosition(rowStopF, 1) yAxesPosition(rowStopF, 1) axisWidth axisHeight]);
      cla
      set(ax(rowStopF, 1), 'ylim', [0 1], 'yticklabel', [], 'xlim', [xMin xMax])
      hold(ax(rowStopF, 1), 'on')
      
      % Left Inhibition Function
      ax(rowStopI, 1) = axes('units', 'centimeters', 'position', [xAxesPosition(rowStopI, 1) yAxesPosition(rowStopI, 1) axisWidth axisHeight]);
      cla
      set(ax(rowStopI, 1), 'ylim', [0 1], 'yticklabel', [], 'xlim', [0 inhMax])
      hold(ax(rowStopI, 1), 'on')
      
      % Right Stop Trials: Correct and Error
      ax(rowStopF, 2) = axes('units', 'centimeters', 'position', [xAxesPosition(rowStopF, 2) yAxesPosition(rowStopF, 2) axisWidth axisHeight]);
      cla
      set(ax(rowStopF, 2), 'ylim', [0 1], 'yticklabel', [], 'xlim', [xMin xMax])
      hold(ax(rowStopF, 2), 'on')
      
      % Right Inhibition Function
      ax(rowStopI, 2) = axes('units', 'centimeters', 'position', [xAxesPosition(rowStopI, 2) yAxesPosition(rowStopI, 2) axisWidth axisHeight]);
      cla
      set(ax(rowStopI, 2), 'ylim', [0 1], 'yticklabel', [], 'xlim', [0 inhMax])
      hold(ax(rowStopI, 2), 'on')
      
   end
   
   
   
   for i = 1 : nSignal/2
      
      % PLOT RIGHT TARGET TRIALS
      iPropIndexR = i + nSignal/2;  % Reverse order of plotting to keep color overlays similar between left and right target
      
      % Target (errors)
%       plot(ax(rowTargF, 2), prdTargT{iPropIndexR}, prdTargCDF{iPropIndexR}, 'color', cMap(iPropIndexR,:), 'lineWidth', lwCorr)
      plot(ax(rowTargF, 2), prdQuantileGoCorr{iPropIndexR}, pQuantile, 'color', cMap(iPropIndexR,:), 'markerSize', rtSize, 'lineWidth', lwCorr)
      plot(ax(rowTargF, 2), obsQuantileGoCorr{iPropIndexR}, pQuantile, 'o', 'color', cMap(iPropIndexR,:), 'markerSize', rtSize, 'lineWidth', lwCorr)

      %       % Distractor (errors)
%       plot(ax(rowTargF, 2), prdDistT{iPropIndexR}, prdDistCDF{iPropIndexR}, '--', 'color', cMap(iPropIndexR,:))
%       plot(ax(rowTargF, 2), obsQuantileGoComm{iPropIndexR}, pQuantile, 'd', 'color', cMap(iPropIndexR,:))
%       
      
      if includeStop
         %Stop to Target
%          plot(ax(rowStopF, 2), prdTargTStop{iPropIndexR}, prdTargCDFStop{iPropIndexR}, 'color', cMap(iPropIndexR,:), 'lineWidth', lwCorr)
         plot(ax(rowStopF, 2), prdQuantileStopFailureCorr{iPropIndexR}, pQuantile, 'color', cMap(iPropIndexR,:), 'markerSize', rtSize, 'lineWidth', lwCorr)
         plot(ax(rowStopF, 2), obsQuantileStopFailureCorr{iPropIndexR}, pQuantile, 'o', 'color', cMap(iPropIndexR,:), 'markerSize', rtSize, 'lineWidth', lwCorr)
         
%          % Stop to Distractor
%          plot(ax(rowStopF, 2), prdDistTStop{iPropIndexR}, prdDistCDFStop{iPropIndexR}, '--', 'color', cMap(iPropIndexR,:))
%          plot(ax(rowStopF, 2), obsQuantileStopFailureComm{iPropIndexR}, pQuantile, 'd', 'color', cMap(iPropIndexR,:))
         
         % Inhibition Functions
%          plot(ax(rowStopI, 2), ssdTimePoints, obsInh(iPropIndexR, :), '--', 'color', cMap(iPropIndexR,:), 'lineWidth', lwCorr)
%          plot(ax(rowStopI, 2), ssdTimePoints, prdInh(iPropIndexR, :), 'color', cMap(iPropIndexR,:), 'lineWidth', lwCorr)
         
          plot(ax(rowStopI, 2), obs.ssd(iPropIndexR,:), obs.pStopFailure(iPropIndexR,:), 'o', 'color', cMap(iPropIndexR,:), 'markerSize', inhSize, 'lineWidth', lwCorr)
          plot(ax(rowStopI, 2), obs.ssd(iPropIndexR,:), prd.pStopFailure(iPropIndexR,:), '-', 'color', cMap(iPropIndexR,:), 'markerSize', inhSize, 'lineWidth', lwCorr)
      end % if includeStop
      
      
      % PLOT LEFT TARGET TRIALS
      iPropIndexL = nSignal/2 + 1 - i;
      
      
      % Target (errors)
%       plot(ax(rowTargF, 1), prdTargT{iPropIndexL}, prdTargCDF{iPropIndexL}, 'color', cMap(iPropIndexL,:), 'lineWidth', lwCorr)
      plot(ax(rowTargF, 1), prdQuantileGoCorr{iPropIndexL}, pQuantile, 'color', cMap(iPropIndexL,:), 'markerSize', rtSize, 'lineWidth', lwCorr)
      plot(ax(rowTargF, 1), obsQuantileGoCorr{iPropIndexL}, pQuantile, 'o', 'color', cMap(iPropIndexL,:), 'markerSize', rtSize, 'lineWidth', lwCorr)
      
%       % Distractor (errors)
%       plot(ax(rowTargF, 1), prdDistT{iPropIndexL}, prdDistCDF{iPropIndexL}, '--', 'color', cMap(iPropIndexL,:), 'lineWidth', lwCorr)
%       plot(ax(rowTargF, 1), obsQuantileGoComm{iPropIndexL}, pQuantile, 'd', 'color', cMap(iPropIndexL,:), 'lineWidth', lwCorr)
      
      
      if includeStop
         %Stop to Target
%          plot(ax(rowStopF, 1), prdTargTStop{iPropIndexL}, prdTargCDFStop{iPropIndexL}, 'color', cMap(iPropIndexL,:), 'lineWidth', lwCorr)
         plot(ax(rowStopF, 1), prdQuantileStopFailureCorr{iPropIndexL}, pQuantile, 'color', cMap(iPropIndexL,:), 'markerSize', rtSize, 'lineWidth', lwCorr)
         plot(ax(rowStopF, 1), obsQuantileStopFailureCorr{iPropIndexL}, pQuantile, 'o', 'color', cMap(iPropIndexL,:), 'markerSize', rtSize, 'lineWidth', lwCorr)
         
%          % Stop to Distractor
%          plot(ax(rowStopF, 1), prdDistTStop{iPropIndexL}, prdDistCDFStop{iPropIndexL}, '--', 'color', cMap(iPropIndexL,:), 'lineWidth', lwCorr)
%          plot(ax(rowStopF, 1), obsQuantileStopFailureComm{iPropIndexL}, pQuantile, 'd', 'color', cMap(iPropIndexL,:), 'lineWidth', lwCorr)
         
         % Inhibition Functions
%          plot(ax(rowStopI, 1), ssdTimePoints, obsInh(iPropIndexL, :), '--', 'color', cMap(iPropIndexL,:), 'lineWidth', lwCorr)
%          plot(ax(rowStopI, 1), ssdTimePoints, prdInh(iPropIndexL, :), 'color', cMap(iPropIndexL,:), 'lineWidth', lwCorr)
         
          plot(ax(rowStopI, 1), obs.ssd(iPropIndexL,:), obs.pStopFailure(iPropIndexL,:), 'o', 'color', cMap(iPropIndexL,:), 'markerSize', inhSize, 'lineWidth', lwCorr)
          plot(ax(rowStopI, 1), obs.ssd(iPropIndexL,:), prd.pStopFailure(iPropIndexL,:), '-', 'color', cMap(iPropIndexL,:), 'markerSize', inhSize, 'lineWidth', lwCorr)
      %
      end % if includeStop
      
      
      
   end % for i = 1 : nSignal/2
   
   
   
   
else
   
   obsColor    = 'k';
   prdColor    = 'r';
   goColor    = 'k';
   stopColor    = 'r';
   
   
   nColumn     = nCond;
   nRow        = 4;   % for now, plot 2 rows-- need to add more with stop trials
   rowTargF     = 1;
   rowStopF     = 2;
   [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
   clf
   for iCond = 1 : nCond
      
      
      
      % Set up plot axes
      % ----------------------------------------------------------------
      
      % Go trials: Correct and Error
      ax(rowTargF, iCond) = axes('units', 'centimeters', 'position', [xAxesPosition(rowTargF, iCond) yAxesPosition(rowTargF, iCond) axisWidth axisHeight]);
      cla
      set(ax(rowTargF, iCond), 'ylim', [0 1], 'xlim', [xMin xMax])
      hold(ax(rowTargF, iCond), 'on')
      %             ttl = sprintf('SSD:  %d', iSSD);
      %             title(ttl)
      
      if includeStop
         
         % Stop Trials: Correct and Error
         ax(rowStopF, iCond) = axes('units', 'centimeters', 'position', [xAxesPosition(rowStopF, iCond) yAxesPosition(rowStopF, iCond) axisWidth axisHeight]);
         cla
         set(ax(rowStopF, iCond), 'ylim', [0 1], 'xlim', [xMin xMax])
         hold(ax(rowStopF, iCond), 'on')
         
      end
      if iCond == 1
         title(ax(rowTargF, iCond), titleMech)
      end
      if iCond > 1
         set(ax(rowTargF, iCond), 'yticklabel', [])
         if includeStop
            set(ax(rowStopF, iCond), 'yticklabel', [])
            set(ax(rowInh, iCond), 'yticklabel', [])
         end
      end
      
      
      
      % Plot the data
      % ----------------------------------------------------------------
      
      % Plot the full CDFs for Target (correct)
      %       plot(ax(rowTargF, iCond), obsTargT{iCond}, obsTargCDF{iCond}, 'b')
      plot(ax(rowTargF, iCond), prdTargT{iCond}, prdTargCDF{iCond}, goColor)
      
      plot(ax(rowTargF, iCond), obsQuantileGoCorr{iCond}, pQuantile, 'o', 'color', goColor)
      %     plot(ax(rowTargF, iCond), prdQuantileGoCorr{iCond}, pQuantile, 'o', 'color', prdColor)
      
      
      
      % Plot the full CDFs for Distacotr (errors)
      %     plot(ax(rowTargF, iCond), obsDistT{iCond}, obsDistCDF{iCond}, obsColor)
      plot(ax(rowTargF, iCond), prdDistT{iCond}, prdDistCDF{iCond}, '--', 'color', goColor)
      
      plot(ax(rowTargF, iCond), obsQuantileGoComm{iCond}, pQuantile, 'd', 'color', goColor)
      %     plot(ax(rowTargF, iCond), prdQuantileGoCorr{iCond}, pQuantile, 'o', 'color', prdColor)
      
      
      if includeStop
         
         
         
         % Plot the full CDFs for Target (correct)
         %     plot(ax(rowStopF, iCond), obsTargTStop{iCond}, obsTargCDFStop{iCond}, obsColor)
         plot(ax(rowStopF, iCond), prdTargTStop{iCond}, prdTargCDFStop{iCond}, stopColor)
         
         plot(ax(rowStopF, iCond), obsQuantileStopFailureCorr{iCond}, pQuantile, 'o', 'color', stopColor)
         %     plot(ax(rowStopF, iCond), prdQuantileStopFailureCorr{iCond}.correct(:,1), pQuantile, 'o', 'color', prdColor)
         
         
         
         % Plot the full CDFs for Distacotr (errors)
         %     plot(ax(rowStopF, iCond), obsDistTStop{iCond}, obsDistCDFStop{iCond}, obsColor)
         plot(ax(rowStopF, iCond), prdDistTStop{iCond}, prdDistCDFStop{iCond}, '--', 'color', stopColor)
         
         plot(ax(rowStopF, iCond), obsQuantileStopFailureComm{iCond}, pQuantile, 'd', 'color', stopColor)
         %     plot(ax(rowStopF, iCond), prdQuantileStopFailureComm{iCond}, pQuantile, 'o', 'color', prdColor)
         
      end
   end % for iCond = 1 : nCond
   
   
end













%% Plot mean RTs:
% ----------------------------------------------------------------
figureHandle = 56;
rtMeanAx = 99;
nRow = 3;
nColumn = 3;
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
clf
axM = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2) yAxesPosition(2, 2) axisWidth axisHeight]);
cla
choicePlotXMargin = .03;
set(axM, 'xlim', [pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
set(axM, 'ylim', rtLimY)
    set(axM, 'xtick', pSignalArray)
    set(axM, 'xtickLabel', pSignalArray*100)
hold(axM, 'on')



% Go trials: Correct
rtObsMeanGoCorr = cellfun(@nanmean, obs.rtGoCorr);
rtPrdMeanGoCorr = cellfun(@mean, prd.rtGoCorr);

rtObsMeanGoComm = cellfun(@nanmean, obs.rtGoComm);
rtPrdMeanGoComm = cellfun(@mean, prd.rtGoComm);

plot(axM, pSignalArray, rtObsMeanGoCorr, 'o', 'markeredgecolor', goColor, 'markerfacecolor', [1 1 1], 'markersize', 12, 'lineWidth', lwCorr)
plot(axM, pSignalArray, rtObsMeanGoComm, 'd', 'markeredgecolor', goColor, 'markerfacecolor', [1 1 1], 'markersize', 8, 'lineWidth', lwCorr)
plot(axM, pSignalArray(1 : nSignal/2), rtPrdMeanGoCorr(1 : nSignal/2), 'color', goColor, 'lineWidth', lwCorr);
plot(axM, pSignalArray(1+nSignal/2 : end), rtPrdMeanGoCorr(1+nSignal/2 : end), 'color', goColor, 'lineWidth', lwCorr);
plot(axM, pSignalArray(1 : nSignal/2), rtPrdMeanGoComm(1 : nSignal/2), '--', 'color', goColor, 'lineWidth', lwCorr);
plot(axM, pSignalArray(1+nSignal/2 : end), rtPrdMeanGoComm(1+nSignal/2 : end), '--', 'color', goColor, 'lineWidth', lwCorr);


if includeStop
   % Stop trials: Correct
   rtObsMeanStopCorr = cellfun(@nanmean, rtObsStopCorr);
   rtPrdMeanStopCorr = cellfun(@nanmean, rtPrdStopCorr);
   
   rtObsMeanStopComm = cellfun(@nanmean, rtObsStopComm);
   rtPrdMeanStopComm = cellfun(@nanmean, rtPrdStopComm);
   
   
   
   plot(axM, pSignalArray, rtObsMeanStopCorr, 'o', 'markeredgecolor', stopColor, 'markerfacecolor', [1 1 1], 'markersize', 12, 'lineWidth', lwCorr)
   plot(axM, pSignalArray(1 : nSignal/2), rtPrdMeanStopCorr(1 : nSignal/2), 'color', stopColor, 'lineWidth', lwCorr);
   plot(axM, pSignalArray(1+nSignal/2 : end), rtPrdMeanStopCorr(1+nSignal/2 : end), 'color', stopColor, 'lineWidth', lwCorr);
   plot(axM, pSignalArray, rtObsMeanStopComm, 'd', 'markeredgecolor', stopColor, 'markerfacecolor', [1 1 1], 'markersize', 8, 'lineWidth', lwCorr)
   plot(axM, pSignalArray(1 : nSignal/2), rtPrdMeanStopComm(1 : nSignal/2), '--', 'color', stopColor, 'lineWidth', lwCorr);
   plot(axM, pSignalArray(1+nSignal/2 : end), rtPrdMeanStopComm(1+nSignal/2 : end), '--', 'color', stopColor, 'lineWidth', lwCorr);
   set(axM, 'xtick', pSignalArray)
   set(axM, 'xtickLabel', pSignalArray*100)
end













%% Plot mean SSRTs:
% ----------------------------------------------------------------

if includeStop

axSSRT = axes('units', 'centimeters', 'position', [xAxesPosition(2, 3) yAxesPosition(2, 3) axisWidth axisHeight]);
cla
set(axSSRT, 'xlim', [pSignalArray(1) - choicePlotXMargin pSignalArray(end) + choicePlotXMargin])
set(axSSRT, 'ylim', ssrtLimY)
    set(axSSRT, 'xtick', pSignalArray)
    set(axSSRT, 'xtickLabel', pSignalArray*100)
hold(axSSRT, 'on')



% iSSD = SAM.des.expt.stimOns{iCnd, iStopTrType}(iSTOP);
% ssd = obs.ssd(1,:);

% SSRTs collapse over SSDs
prdRTStopSuccess = cell(nSignal, 1);
for iCond = 1 : nSignal
   prdRTStopSuccess{iCond} = cell2mat(prd.rtStopSuccess(iCond,:));
end
prdSsrtStopSuccelsMean = cellfun(@nanmean,prdRTStopSuccess);
prdSsrtStopSuccessStd = cellfun(@nanstd,prdRTStopSuccess);


plot(axSSRT, pSignalArray, prdSsrtStopSuccelsMean, 'o', 'markeredgecolor', goColor, 'markerfacecolor', [1 1 1], 'markersize', 12, 'lineWidth', lwCorr)
    errorbar(axSSRT, pSignalArray ,prdSsrtStopSuccelsMean, prdSsrtStopSuccessStd, 'linestyle' , 'none', 'color', goColor, 'linewidth' , 2)

% plot(axSSRT, pSignalArray, rtObsMeanGoComm, 'd', 'markeredgecolor', goColor, 'markerfacecolor', [1 1 1], 'markersize', 8, 'lineWidth', lwCorr)
% plot(axSSRT, pSignalArray(1 : nSignal/2), rtPrdMeanGoCorr(1 : nSignal/2), 'color', goColor, 'lineWidth', lwCorr);
% plot(axSSRT, pSignalArray(1+nSignal/2 : end), rtPrdMeanGoCorr(1+nSignal/2 : end), 'color', goColor, 'lineWidth', lwCorr);
% plot(axSSRT, pSignalArray(1 : nSignal/2), rtPrdMeanGoComm(1 : nSignal/2), '--', 'color', goColor, 'lineWidth', lwCorr);
% plot(axSSRT, pSignalArray(1+nSignal/2 : end), rtPrdMeanGoComm(1+nSignal/2 : end), '--', 'color', goColor, 'lineWidth', lwCorr);
% 
% 
% if includeStop
%    % Stop trials: Correct
%    rtObsMeanStopCorr = cellfun(@mean, rtObsStopCorr);
%    rtPrdMeanStopCorr = cellfun(@mean, rtPrdStopCorr);
%    
%    rtObsMeanStopComm = cellfun(@mean, rtObsStopComm);
%    rtPrdMeanStopComm = cellfun(@mean, rtPrdStopComm);
%    
%    
%    
%    plot(axSSRT, pSignalArray, rtObsMeanStopCorr, 'o', 'markeredgecolor', stopColor, 'markerfacecolor', [1 1 1], 'markersize', 12, 'lineWidth', lwCorr)
%    plot(axSSRT, pSignalArray(1 : nSignal/2), rtPrdMeanStopCorr(1 : nSignal/2), 'color', stopColor, 'lineWidth', lwCorr);
%    plot(axSSRT, pSignalArray(1+nSignal/2 : end), rtPrdMeanStopCorr(1+nSignal/2 : end), 'color', stopColor, 'lineWidth', lwCorr);
%    plot(axSSRT, pSignalArray, rtObsMeanStopComm, 'd', 'markeredgecolor', stopColor, 'markerfacecolor', [1 1 1], 'markersize', 8, 'lineWidth', lwCorr)
%    plot(axSSRT, pSignalArray(1 : nSignal/2), rtPrdMeanStopComm(1 : nSignal/2), '--', 'color', stopColor, 'lineWidth', lwCorr);
%    plot(axSSRT, pSignalArray(1+nSignal/2 : end), rtPrdMeanStopComm(1+nSignal/2 : end), '--', 'color', stopColor, 'lineWidth', lwCorr);
%    set(axSSRT, 'xtick', pSignalArray)
%    set(axSSRT, 'xtickLabel', pSignalArray*100)
% end
end
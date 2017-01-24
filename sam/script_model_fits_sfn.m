%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTINGS (THESE MAY BE ADJUSTED)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables;close all;clc;

% Fixed variables
% =========================================================================

% Column indices
% -------------------------------------------------------------------------
iGoCorr         = 1;
iGoComm         = 2;
iStopFail       = 4; % SSD3

% Quantiles
% -------------------------------------------------------------------------
qntlsPrd        = [0.01:0.01:0.99,0.9999];
qntlsObs        = 0.1:0.1:0.9;

% Plotting
% -------------------------------------------------------------------------
clrCnd        = {[0 1 1],[1 0 1],[0.5 0 1]};
lnCndCorr     = {'-','-','-'};
lnCndError    = {'--','--','--'};

mrkSize         = 5;
mrkCndCorr    = 'o';
mrkCndError   = '^';

% clrGoCorr       = [0 0.5 0];
% clrGoCorrTransp = [0 0.25 0];
% clrGoComm       = [0.5 0.5 0];
% clrGoCommTransp = [0.25 0.25 0];
% clrStop         = [1 0 0];
% clrStopTransp   = [0.5 0 0];
% lnGoCorr        = '-';
% lnGoCommS       = '-';
% lnGoCommQ       = '--';
% lnStop          = '-';


clrGoCorr       = {[0 1 1],[1 0 1],[0.5 0 1]};
clrGoComm       = {[0 1 1],[1 0 1],[0.5 0 1]};
clrStopFail     = {[0 1 1],[1 0 1],[0.5 0 1]};
clrStopSucc     = {[0 1 1],[1 0 1],[0.5 0 1]};

lnGoCorr        = {'-','-','-'};
lnGoComm        = {'--','--','--'};
lnStopFail      = {'--','--','--'};
lnStopSucc      = {'-','-','-'};
lnGoCommQ        = {'--','--','--'};
lnStopFail      = {'-','-','-'};
lnStopSucc      = {'-','-','-'};
% 
% mrkSize         = 5;
% mrkTypeGoCorr   = 'o';
% mrkTypeGoComm   = '^';
% mrkTypeStopFail = 's';
% mrkTypeStopSucc = 's';

subj            = 8:13;
nSubj           = 6;
nCnd            = 3;
nSsd            = 5;
nQPrd           = numel(qntlsPrd);
nQObs           = numel(qntlsObs);

% Figure panels
figure;
p               = panel;
p.margin        = [5 5 2 2];
p.pack(8,3);

% Figure directory (for dynamics only)
figDir          = '/Users/bramzandbelt/Dropbox/SAM/Neuroscience/Figures/Raw/ModelDynamics';

% Figure layout
dynType         = 'qaverage';    % 'indiv' or 'qaverage'
normalizeDyn    = false;
xDataLim        = [0 1200];
yDataLimRt      = [0 1];
if normalizeDyn
  yDataLimDyn     = [0 1.25];
else
  yDataLimDyn     = [0 100];
end

edgeAlpha       = 0.1;


% Dynamic variables
% =========================================================================

% Choice mechanism
% -------------------------------------------------------------------------
choiceMechType  = 'race';

% Inhibition mechanism
% -------------------------------------------------------------------------
inhibMechType   = 'race';

% Parameter varying across conditions
% -------------------------------------------------------------------------
condParam       = {'t0','v','zc'};

for iCondParam = 1:numel(condParam)

  fprintf('Working on condParam %d ... \n',iCondParam);
  
  % Tag
  % -------------------------------------------------------------------------
  figTag = [choiceMechType(1),inhibMechType(1),condParam{iCondParam}];
  
  % Example subject index for trial dynamics
  % -----------------------------------------------------------------------
  switch condParam{iCondParam}
    case 't0'
      iExSubj = 11;
    case 'v'
%       iExSubj = 8;
      iExSubj = 12;
    case 'zc'
      iExSubj = 10;
  end
  
  % Data file
  % -------------------------------------------------------------------------
  dataFile        = sprintf('/Users/bramzandbelt/Dropbox/SAM/data/qaveragedata/qaverage_data_allsubjects_alltrials_c%s_i%s_p%s.mat',choiceMechType,inhibMechType,condParam);
%   dataFile        = sprintf('/Users/bramzandbelt/Documents/PROJECTS/SAM/output/qaverage_wsmall_error_in_iCnd2/qaverage_data_allsubjects_alltrials_c%s_i%s_p%s.mat',choiceMechType,inhibMechType,condParam{iCondParam});

  % Plot type
  % -------------------------------------------------------------------------
  plotType        = 'qaverage';     % 'indiv' or 'qaverage'

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PLOTTING
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Load data file
  % =========================================================================
  load(dataFile);

  % Panel 1 - RT distributions: correct go and choice error
  % =========================================================================
  
  fprintf('Panel 1 -  RT distributions: correct go and choice error \n');
  
  % Pre-allocate matrices for logging
  % -------------------------------------------------------------------------
  rtGoCorrPrd = nan(nCnd,nQPrd,nSubj);
  rtGoCorrObs = nan(nCnd,nQObs,nSubj);

  rtGoCommPrd = nan(nCnd,nQPrd,nSubj);
  rtGoCommObs = nan(nCnd,nQObs,nSubj);

  % Extract observed and predicted quantiles
  % -------------------------------------------------------------------------
  for iCnd = 1:nCnd
    for iS = 1:nSubj

      % Subject index
      iSubj = subj(iS);

      % Observations
      rtGoCorrObs(iCnd,:,iS) = quantile(obsData{iSubj}.rt{iCnd,iGoCorr},qntlsObs);

      rtGoCommObs(iCnd,:,iS) = quantile(obsData{iSubj}.rt{iCnd,iGoComm},qntlsObs);

      % Predictions
      rtGoCorrPrd(iCnd,:,iS) = quantile(prdData{iSubj}.rt{iCnd,iGoCorr},qntlsPrd);

      rtGoCommPrd(iCnd,:,iS) = quantile(prdData{iSubj}.rt{iCnd,iGoComm},qntlsPrd);

    end
  end

  % Compute quantile averages (average across subjects)
  % -------------------------------------------------------------------------
  switch plotType
    case 'qaverage'

      qRtGoCorrObs = nanmean(rtGoCorrObs,3);
      qRtGoCorrPrd = nanmean(rtGoCorrPrd,3);

      qRtGoCommObs = nanmean(rtGoCommObs,3);
      qRtGoCommPrd = nanmean(rtGoCommPrd,3);

  end

  % Plot observations and predictions
  % -------------------------------------------------------------------------
  switch plotType
    case 'indiv'

      for iS = 1:nSubj

        figure; hold on;
        set(gcf,'Name',sprintf('c%s, i%s, p%s, subject %d',choiceMechType,inhibMechType,condParam{iCondParam},subj(iS)));

        for iCnd = 1:nCnd

          % Plot GoCorr data
          col   = clrCnd{iCnd};
          ln    = lnCndCorr{iCnd};
          mrk   = mrkCndCorr;
          plot(rtGoCorrPrd(iCnd,:,iS),qntlsPrd,'Color',col,'LineStyle',ln,'LineWidth',0.5);
          scatter(rtGoCorrObs(iCnd,:,iS),qntlsObs,mrkSize,col,mrk);

          % Plot GoComm data
          col   = clrCnd{iCnd};
          ln    = lnCndError{iCnd};
          mrk   = mrkCndError;
          plot(rtGoCommPrd(iCnd,:,iS),qntlsPrd,'Color',col,'LineStyle',ln,'LineWidth',0.5);
          scatter(rtGoCommObs(iCnd,:,iS),qntlsObs,mrkSize,col,mrk);

        end
      end


    case 'qaverage'

      % Select panel
      % -------------------------------------------------------------------
      p(1,iCondParam).select();
      p(1,iCondParam).hold('on');
      
      for iCnd = 1:3

        % Plot GoCorr data
        col   = clrCnd{iCnd};
        ln    = lnCndCorr{iCnd};
        mrk   = mrkCndCorr;
        plot(qRtGoCorrPrd(iCnd,:),qntlsPrd,'Color',col,'LineStyle',ln,'LineWidth',0.5);
        scatter(qRtGoCorrObs(iCnd,:),qntlsObs,mrkSize,col,mrk);

        % Plot GoComm data
        col   = clrCnd{iCnd};
        ln    = lnCndError{iCnd};
        mrk   = mrkCndError;
        plot(qRtGoCommPrd(iCnd,:),qntlsPrd,'Color',col,'LineStyle',ln,'LineWidth',0.5);
        scatter(qRtGoCommObs(iCnd,:),qntlsObs,mrkSize,col,mrk);

      end

  end

  % Tag the panel
  text(xDataLim(2),0.2,figTag,'FontSize',8,'HorizontalAlignment','right');
  
  % Figure layout
%   title('Correct and choice error RTs');
  set(gca,'PlotBoxAspectRatio',[1.61,1,1],'XLim',xDataLim,'XTick',[],'YLim',yDataLimRt,'YTick',[]);

  % Panel 2 - RT distributions: correct go and inhibition error
  % =========================================================================

  fprintf('Panel 2 - RT distributions: correct go and inhibition error \n');
  
  % Pre-allocate matrices for logging
  % -------------------------------------------------------------------------
  rtGoCorrPrd = nan(nCnd,nQPrd,nSubj);
  rtGoCorrObs = nan(nCnd,nQObs,nSubj);

  rtStopFailPrd = nan(nCnd,nQPrd,nSubj);
  rtStopFailObs = nan(nCnd,nQObs,nSubj);

  % Extract observed and predicted quantiles
  % -------------------------------------------------------------------------
  for iCnd = 1:3
    for iS = 1:6

      % Subject index
      iSubj = subj(iS);

      % Observations
      rtGoCorrObs(iCnd,:,iS) = quantile(obsData{iSubj}.rt{iCnd,iGoCorr},qntlsObs);

      rtStopFailObs(iCnd,:,iS) = quantile(obsData{iSubj}.rt{iCnd,iStopFail},qntlsObs);

      % Predictions
      rtGoCorrPrd(iCnd,:,iS) = quantile(prdData{iSubj}.rt{iCnd,iGoCorr},qntlsPrd);

      rtStopFailPrd(iCnd,:,iS) = quantile(prdData{iSubj}.rt{iCnd,iStopFail},qntlsPrd);

    end
  end

  % Compute quantile averages (average across subjects)
  % -------------------------------------------------------------------------
  switch plotType
    case 'qaverage'

      qRtGoCorrObs = nanmean(rtGoCorrObs,3);
      qRtGoCorrPrd = nanmean(rtGoCorrPrd,3);

      qRtStopFailObs = nanmean(rtStopFailObs,3);
      qRtStopFailPrd = nanmean(rtStopFailPrd,3);

  end

  % Plot observations and predictions
  % -------------------------------------------------------------------------
  switch plotType
    case 'indiv'

      for iS = 1:nSubj

        figure; hold on;
        set(gcf,'Name',sprintf('c%s, i%s, p%s, subject %d',choiceMechType,inhibMechType,condParam{iCondParam},subj(iS)));

        for iCnd = 1:nCnd
  
          % Plot GoCorr data
          col   = clrCnd{iCnd};
          ln    = lnCndCorr{iCnd};
          mrk   = mrkCndCorr;
          plot(rtGoCorrPrd(iCnd,:,iS),qntlsPrd,'Color',col,'LineStyle',ln,'LineWidth',0.5);
          scatter(rtGoCorrObs(iCnd,:,iS),qntlsObs,mrkSize,col,mrk);

          % Plot StopFail data
          col   = clrCnd{iCnd};
          ln    = lnCndError{iCnd};
          mrk   = mrkCndError;
          plot(rtStopFailPrd(iCnd,:,iS),qntlsPrd,'Color',col,'LineStyle',ln,'LineWidth',0.5);
          scatter(rtStopFailObs(iCnd,:,iS),qntlsObs,mrkSize,col,mrk);

        end
      end

    case 'qaverage'

      % Select panel
      % -------------------------------------------------------------------
      p(2,iCondParam).select();
      p(2,iCondParam).hold('on');

      for iCnd = 1:3

        % Plot GoCorr data
        col   = clrCnd{iCnd};
        ln    = lnCndCorr{iCnd};
        mrk   = mrkCndCorr;
        plot(qRtGoCorrPrd(iCnd,:),qntlsPrd,'Color',col,'LineStyle',ln,'LineWidth',0.5);
        scatter(qRtGoCorrObs(iCnd,:),qntlsObs,mrkSize,col,mrk);

        % Plot StopFail data
        col   = clrCnd{iCnd};
        ln    = lnCndError{iCnd};
        mrk   = mrkCndError;
        plot(qRtStopFailPrd(iCnd,:),qntlsPrd,'Color',col,'LineStyle',ln,'LineWidth',0.5);
        scatter(qRtStopFailObs(iCnd,:),qntlsObs,mrkSize,col,mrk);

      end

  end

  % Tag the panel
  text(xDataLim(2),0.2,figTag,'FontSize',8,'HorizontalAlignment','right');
  
  % Figure layout
%   title('Correct and inhibition error RTs');
  set(gca,'PlotBoxAspectRatio',[1.61,1,1],'XLim',xDataLim,'XTick',[],'YLim',yDataLimRt,'YTick',[]);

  % Panel 3 - Inhibition function
  % =========================================================================

  fprintf('Panel 3 - Inhibition function \n');
  
  % Pre-allocate matrices for logging
  % -------------------------------------------------------------------------
  ssd         = nan(nCnd,nSsd,nSubj);
  inhibFunPrd = nan(nCnd,nSsd,nSubj);
  inhibFunObs = nan(nCnd,nSsd,nSubj);

  % Extract observed and predicted quantiles
  % -------------------------------------------------------------------------

  for iS = 1:6

    % Subject index
    iSubj = subj(iS);

    % Load observations (for SSD)
    obsFile = sprintf('/Users/bramzandbelt/Dropbox/SAM/data/subj%.2d/obs.mat',iSubj);
%     obsFile = sprintf('/Users/bramzandbelt/Documents/Dropbox/SAM/data/subj%.2d/obs.mat',iSubj);
    load(obsFile);

    % SSDs
    ssd(:,:,iS)         = obs.ssd;

    % Observations
    inhibFunObs(:,:,iS) = obsData{iSubj}.P(:,3:end);

    % Predictions
    inhibFunPrd(:,:,iS) = prdData{iSubj}.P(:,3:end);

  end


  % Compute quantile averages (average across subjects)
  % -------------------------------------------------------------------------
  switch plotType
    case 'qaverage'

        meanSsd         = mean(ssd,3);
        meanInhibFunObs = mean(inhibFunObs,3);
        meanInhibFunPrd = mean(inhibFunPrd,3);

  end

  % Plot observations and predictions
  % -------------------------------------------------------------------------
  switch plotType

    case 'indiv'

      for iS = 1:nSubj

        figure; hold on;
        set(gcf,'Name',sprintf('c%s, i%s, p%s, subject %d',choiceMechType,inhibMechType,condParam{iCondParam},subj(iS)));

        for iCnd = 1:nCnd

          % Plot inhibition functions
          col   = clrCnd{iCnd};
          ln    = lnCndError{iCnd};
          mrk   = mrkCndError;
          plot(ssd(iCnd,:,iS),inhibFunPrd(iCnd,:,iS),'Color',col,'LineStyle',ln,'LineWidth',0.5);
          scatter(ssd(iCnd,:,iS),inhibFunObs(iCnd,:,iS),mrkSize,col,mrk);
          
        end
      end

    case 'qaverage'

      % Select panel
      % -------------------------------------------------------------------
      p(3,iCondParam).select();
      p(3,iCondParam).hold('on');

      for iCnd = 1:3

        % Plot inhibition functions
        col   = clrCnd{iCnd};
        ln    = lnCndError{iCnd};
        mrk   = mrkCndError;
        plot(meanSsd(iCnd,:),meanInhibFunPrd(iCnd,:),'Color',col,'LineStyle',ln,'LineWidth',0.5);
        scatter(meanSsd(iCnd,:),meanInhibFunObs(iCnd,:),mrkSize,col,mrk);
  
      end
  end

  % Tag the panel
  text(xDataLim(2),0.2,figTag,'FontSize',8,'HorizontalAlignment','right');
  
  % Figure layout
%   title('Inhibition function');
  set(gca,'PlotBoxAspectRatio',[1.61,1,1],'XLim',xDataLim,'XTick',[],'YLim',yDataLimRt,'YTick',[]);

  % Panel 4 - SSRT distributions
  % =========================================================================

  fprintf('Panel 4 - SSRT distributions \n');
  
  % Pre-allocate matrices for logging
  % -------------------------------------------------------------------------
  ssrtPrd     = nan(nCnd,nQPrd,nSubj);

  % Extract observed and predicted quantiles
  % -------------------------------------------------------------------------
  for iCnd = 1:3
    for iS = 1:6

      % Subject index
      iSubj = subj(iS);

      % Predictions (for one SSD)
      ssrtPrd(iCnd,:,iS) = quantile(prd{iSubj}.rtStopSuccess{iCnd,iStopFail-2},qntlsPrd);

    end
  end

  % Compute quantile averages (average across subjects)
  % -------------------------------------------------------------------------
  switch plotType
    case 'qaverage'

      qSsrtPrd     = nanmean(ssrtPrd,3);

  end

  % Plot observations and predictions
  % -------------------------------------------------------------------------
  switch plotType

    case 'indiv'

      for iS = 1:nSubj

        figure; hold on;
        set(gcf,'Name',sprintf('c%s, i%s, p%s, subject %d',choiceMechType,inhibMechType,condParam{iCondParam},subj(iS)));

        for iCnd = 1:nCnd

          % Plot predicted SSRT distributions
          col   = clrCnd{iCnd};
          ln    = lnCndCorr{iCnd};
          mrk   = mrkCndCorr;
          plot(ssrtPrd(iCnd,:,iS),qntlsPrd,'Color',col,'LineStyle',ln,'LineWidth',0.5);

        end
      end

    case 'qaverage'

      % Select panel
      % -------------------------------------------------------------------
      p(4,iCondParam).select();
      p(4,iCondParam).hold('on');

      for iCnd = 1:3

        % Plot inhibition functions
        col   = clrCnd{iCnd};
        ln    = lnCndCorr{iCnd};
        mrk   = mrkCndCorr;
        plot(qSsrtPrd(iCnd,:),qntlsPrd,'Color',col,'LineStyle',ln,'LineWidth',0.5);

      end

  end

  % Tag the panel
  text(xDataLim(2),0.2,figTag,'FontSize',8,'HorizontalAlignment','right');
  
  % Figure layout
%   title('SSRT distributions');
  set(gca,'PlotBoxAspectRatio',[1.61,1,1],'XLim',xDataLim,'XTick',[],'YLim',yDataLimRt,'YTick',[]);

  % for iCnd = 1:3
  %   for iS = 1:6
  % 
  %     iSubj = subj(iS);
  %     
  %     % Load observations (for SSD)
  %     obsFile = sprintf('/Users/bramzandbelt/Dropbox/SAM/data/subj%.2d/obs.mat',iSubj);
  %     load(obsFile);
  %    
  %     % Compute SSRT according to integration method
  %     pStopFail   = obsData{iSubj}.P(iCnd,iStopFail);
  %     ssd         = obs.ssd(iCnd,iStopFail-2);
  %     finishTime  = quantile(obsData{iSubj}.rt{iCnd,iGoCorr},pStopFail);
  %     ssrt        = finishTime - ssd;
  %        
  %   end
  % end

  % Panel 5 - Model dynamics: correct no-signal
  % =========================================================================

  fprintf('Panel 5 - Model dynamics: correct no-signal \n');
  
  iGoTrType = 1;

  % Pre-allocate matrices for logging
  % -------------------------------------------------------------------------
  switch dynType
    case 'qaverage'
      qXGOT   = nan(nCnd,nQPrd,nSubj);
      qYGOT   = nan(nCnd,nQPrd,nSubj);
      qXGONT  = nan(nCnd,nQPrd,nSubj);
      qYGONT  = nan(nCnd,nQPrd,nSubj);
  end
  

  % Extract, quantile average, and normalize predicted dynamics
  % -------------------------------------------------------------------------
  for iCnd = 1:3
    
    switch dynType
      case 'qaverage'
        for iS = 1:6

          % Subject index
          iSubj = subj(iS);

          % Get quantile averaged dynamics
          xI = prd{iSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GOT.qX;
          yI = prd{iSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GOT.qY;
          xO = quantile(prd{iSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GOT.qX,qntlsPrd);
          yO = interp1(xI,yI,xO);
          
          qXGOT(iCnd,:,iS) = xO;
          qYGOT(iCnd,:,iS) = yO;
          
          xI = prd{iSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GONT.qX;
          yI = prd{iSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GONT.qY;
          xO = quantile(prd{iSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GONT.qX,qntlsPrd);
          yO = interp1(xI,yI,xO);
          
          qXGONT(iCnd,:,iS) = xO;
          qYGONT(iCnd,:,iS) = yO;
          
          % Normalize dynamics
          if normalizeDyn
            zcGo = max(cell2mat(modelMat{iSubj}.ZC(:))); % Maximum threshold (for normalization)
            qYGOT(iCnd,:,iS) = qYGOT(iCnd,:,iS)./zcGo;
            qYGONT(iCnd,:,iS) = qYGONT(iCnd,:,iS)./zcGo;
          end

        end
        
        % Compute group average
        meanQXGOT = nanmean(qXGOT,3);
        meanQYGOT = nanmean(qYGOT,3);
        meanQXGONT = nanmean(qXGONT,3);
        meanQYGONT = nanmean(qYGONT,3);
        
      case 'indiv' 
        
        if iCnd == 2
        
          % Get individual trial dynamics
          sXGOT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GOT.sX;
          sYGOT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GOT.sY;
          sXGONT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GONT.sX;
          sYGONT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GONT.sY;

          % Get quantile averaged dynamics
          qXGOT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GOT.qX;
          qYGOT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GOT.qY;
          qXGONT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GONT.qX;
          qYGONT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoCorr.goStim.GONT.qY;

          % Normalize dynamics
          if normalizeDyn
            zcGo = modelMat{iExSubj}.ZC{iCnd}(3);

            sYGOT = cellfun(@(a) a./zcGo,sYGOT,'Uni',0);
            sYGONT = cellfun(@(a) a./zcGo,sYGONT,'Uni',0);

            qYGOT = qYGOT./zcGo;
            qYGONT = qYGONT./zcGo;
          end
        
        end
        
    end
  end
  
  % Plot predictions
  % -------------------------------------------------------------------------
  switch plotType
    case 'indiv'

      for iS = 1:nSubj

        figure; hold on;
        set(gcf,'Name',sprintf('c%s, i%s, p%s, subject %d',choiceMechType,inhibMechType,condParam{iCondParam},subj(iS)));

        for iCnd = 1:nCnd

          % Plot dynamics target GO
          col = clrGoCorr;
          ln = lnGoCorr;
          plot(qXGOT(iCnd,:,iS),qYGOT(iCnd,:,iS),'Color',col,'LineStyle',ln,'LineWidth',0.5);

          % Plot dynamics non-target GO
          col = clrGoComm;
          ln = lnGoCommQ;
          plot(qXGONT(iCnd,:,iS),qYGONT(iCnd,:,iS),'Color',col,'LineStyle',ln,'LineWidth',0.5);

        end
      end

    case 'qaverage'

      switch dynType
        
        case 'qaverage'
          
          % Select panel
          % ---------------------------------------------------------------
          p(5,iCondParam).select();
          p(5,iCondParam).hold('on');
      
          for iCnd = 1:3

              % Plot dynamics target GO
              col = clrGoCorr{iCnd};
              ln = lnGoCorr{iCnd};
              plot(meanQXGOT(iCnd,:),meanQYGOT(iCnd,:),'Color',col,'LineStyle',ln,'LineWidth',0.5);


              % Plot dynamics non-target GO
              col = clrGoComm{iCnd};
              ln = lnGoCommQ{iCnd};
              plot(meanQXGONT(iCnd,:),meanQYGONT(iCnd,:),'Color',col,'LineStyle',ln,'LineWidth',0.5);

          end
          
          % Plot threshold
          zCGoAllSubj = nan(nSubj,nCnd);
          for iS = 1:nSubj

            % Subject index
            iSubj = subj(iS);
            allZc = cell2mat(modelMat{iSubj}.ZC(:)');

            if normalizeDyn
              zcGo = max(max(allZc(1:6,:)));
              zCGoAllSubj(iS,:) = allZc(1,:)./zcGo;
            else
              zCGoAllSubj(iS,:) = allZc(1,:);
            end
          end  

          meanZcGo = nanmean(zCGoAllSubj,1);
          
          
          if iCondParam == 3

            for iCnd = 1:nCnd
              col = clrGoCorr{iCnd};
              line([xDataLim(1),xDataLim(2)],[meanZcGo(iCnd),meanZcGo(iCnd)], 'Color',col,'LineWidth',0.5);
            end

          else

            if normalizeDyn
              line([xDataLim(1),xDataLim(2)],[1 1], 'Color','k','LineWidth',0.5);
            else
              line([xDataLim(1),xDataLim(2)],[meanZcGo(1) meanZcGo(1)], 'Color','k','LineWidth',0.5);
            end

          end
          
          % Tag the panel
          text(xDataLim(2),0.2,figTag,'FontSize',8,'HorizontalAlignment','right');

        case 'indiv'
            
            if iCondParam == 2
              hFig = set_figure({20,20/1.61,'centimeters'},{'USLetter','landscape'});
              set(gcf,'PaperPositionMode','auto')
              pnl = panel;
              pnl.margin = [0 0 0 0];
              hold on;

              iCnd = 2;

              % Plot data
              col = clrGoCorrTransp;
              ln = lnGoCorr;
  %             cellfun(@(a,b) patchline(a,b,'EdgeColor',col,'LineWidth',1,'LineStyle',ln,'EdgeAlpha',edgeAlpha),sXGOT,sYGOT,'Uni',0);
              cellfun(@(a,b) plot(a,b,'Color',col,'LineWidth',1,'LineStyle',ln),sXGOT,sYGOT,'Uni',0);

              col = clrGoCommTransp;
              ln = lnGoCommS;
  %             cellfun(@(a,b) patchline(a,b,'EdgeColor',col,'LineWidth',1,'LineStyle',ln,'EdgeAlpha',edgeAlpha),sXGONT,sYGONT,'Uni',0);
              cellfun(@(a,b) plot(a,b,'Color',col,'LineWidth',1,'LineStyle',ln),sXGONT,sYGONT,'Uni',0);

              col = clrGoCorr;
              ln = lnGoCorr;
              plot(qXGOT,qYGOT,'Color',col,'LineWidth',3,'LineStyle',ln);

              col = clrGoComm;
              ln = lnGoCommQ;
              plot(qXGONT,qYGONT,'Color',col,'LineWidth',3,'LineStyle',ln);

              % Plot the threshold
              if normalizeDyn
                line([xDataLim(1),xDataLim(2)],[1 1],'Color','k','LineWidth',3);
              else
                thisThres = modelMat{iExSubj}.ZC{iCnd}(1); %Go threshold
                line([xDataLim(1),xDataLim(2)],[thisThres thisThres],'Color','k','LineWidth',3);
              end

              set(gca,'LineWidth',2);
              
              % File name
              if normalizeDyn
                normStr = 'normThres';
              else
                normStr = 'absThres';
              end
              fName = sprintf('DynamicsCorrectChoice_c%s_i%s_condParam%s_subj%.2d_cond%.2d_%s.eps',choiceMechType,inhibMechType,condParam{iCnd},iExSubj,iCnd,normStr);
  
            end
      end 
  end
  
  % Figure layout
  set(gca,'PlotBoxAspectRatio',[1.61,1,1],'XLim',xDataLim,'XTick',[],'YLim',yDataLimDyn,'YTick',[]);
%   ylabel('Activation (a.u.)');
%   title('Dynamics on correct no-signal trials')
  
  if iCondParam == 2
    switch dynType
      case 'indiv'
        print('-depsc',fullfile(figDir,fName));
    end
  end
  
  % Panel 6 - Model dynamics: incorrect no-signal (choice error)
  % =========================================================================

  fprintf('Panel 6 - Model dynamics: incorrect no-signal (choice error) \n');
  
  iGoTrType = 1;

  % Pre-allocate matrices for logging
  % -------------------------------------------------------------------------
  qXGOT = nan(nCnd,nQPrd,nSubj);
  qYGOT = nan(nCnd,nQPrd,nSubj);
  qXGONT = nan(nCnd,nQPrd,nSubj);
  qYGONT = nan(nCnd,nQPrd,nSubj);

  % Extract, quantile average, and normalize predicted dynamics
  % -------------------------------------------------------------------------
  for iCnd = 1:3
    
    switch dynType
      case 'qaverage'
        for iS = 1:6

          % Subject index
          iSubj = subj(iS);

          try % Try-catch, because there may be no errors

            % Get quantile averaged dynamics
            xI = prd{iSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GOT.qX;
            yI = prd{iSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GOT.qY;
            xO = quantile(prd{iSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GOT.qX,qntlsPrd);
            yO = interp1(xI,yI,xO);

            qXGOT(iCnd,:,iS) = xO;
            qYGOT(iCnd,:,iS) = yO;

            xI = prd{iSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GONTE.qX;
            yI = prd{iSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GONTE.qY;
            xO = quantile(prd{iSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GONTE.qX,qntlsPrd);
            yO = interp1(xI,yI,xO);

            qXGONTE(iCnd,:,iS) = xO;
            qYGONTE(iCnd,:,iS) = yO;
            
            % Normalize dynamics
            if normalizeDyn
              zcGo = max(cell2mat(modelMat{iSubj}.ZC(:))); % Maximum threshold (for normalization)
              qYGOT(iCnd,:,iS) = qYGOT(iCnd,:,iS)./zcGo;
              qYGONTE(iCnd,:,iS) = qYGONTE(iCnd,:,iS)./zcGo;
            end

          catch
          end

        end
        
        % Compute group average
        meanQXGOT = nanmean(qXGOT,3);
        meanQYGOT = nanmean(qYGOT,3);
        meanQXGONTE = nanmean(qXGONTE,3);
        meanQYGONTE = nanmean(qYGONTE,3);
        
      case 'indiv' 
        
        if iCnd == 2
        
          % Get individual trial dynamics
          sXGOT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GOT.sX;
          sYGOT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GOT.sY;
          sXGONTE = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GONTE.sX;
          sYGONTE = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GONTE.sY;

          % Get quantile averaged dynamics
          qXGOT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GOT.qX;
          qYGOT = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GOT.qY;
          qXGONTE = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GONTE.qX;
          qYGONTE = prd{iExSubj}.dyn{iCnd,iGoTrType}.GoComm.goStim.GONTE.qY;

          % Normalize dynamics
          if normalizeDyn
            zcGo = modelMat{iExSubj}.ZC{iCnd}(3);

            sYGOT = cellfun(@(a) a./zcGo,sYGOT,'Uni',0);
            sYGONTE = cellfun(@(a) a./zcGo,sYGONTE,'Uni',0);

            qYGOT = qYGOT./zcGo;
            qYGONTE = qYGONTE./zcGo;
          end
        
        end
    end
  end
    
  % Plot predictions
  % -------------------------------------------------------------------------
  switch plotType
    case 'indiv'

      for iS = 1:nSubj

        figure; hold on;
        set(gcf,'Name',sprintf('c%s, i%s, p%s, subject %d',choiceMechType,inhibMechType,condParam{iCondParam},subj(iS)));

        for iCnd = 1:nCnd

          % Plot dynamics target GO
          col = clrGoCorr{iCnd};
          ln = lnGoCorr{iCnd};
          plot(qXGOT(iCnd,:,iS),qYGOT(iCnd,:,iS),'Color',col,'LineStyle',ln,'LineWidth',0.5);


          % Plot dynamics non-target GO
          col = clrGoComm{iCnd};
          ln = lnGoCommQ{iCnd};
          plot(qXGONTE(iCnd,:,iS),qYGONTE(iCnd,:,iS),'Color',col,'LineStyle',ln,'LineWidth',0.5);

        end
      end

    case 'qaverage'
      
      switch dynType
        
        case 'qaverage'
      

          % Select panel
          % -------------------------------------------------------------------
          p(6,iCondParam).select();
          p(6,iCondParam).hold('on');

          for iCnd = 1:3

              % Plot dynamics target GO
              col = clrGoCorr{iCnd};
              ln = lnGoCorr{iCnd};
              plot(meanQXGOT(iCnd,:),meanQYGOT(iCnd,:),'Color',col,'LineStyle',ln,'LineWidth',0.5);


              % Plot dynamics non-target GO
              col = clrGoComm{iCnd};
              ln = lnGoCommQ{iCnd};
              plot(meanQXGONTE(iCnd,:),meanQYGONTE(iCnd,:),'Color',col,'LineStyle',ln,'LineWidth',0.5);

          end

          % Plot threshold
          zCGoAllSubj = nan(nSubj,nCnd);
          for iS = 1:nSubj

            % Subject index
            iSubj = subj(iS);
            allZc = cell2mat(modelMat{iSubj}.ZC(:)');

            if normalizeDyn
              zcGo = max(max(allZc(1:6,:)));
              zCGoAllSubj(iS,:) = allZc(1,:)./zcGo;
            else
              zCGoAllSubj(iS,:) = allZc(1,:);
            end
          end  

          meanZcGo = nanmean(zCGoAllSubj,1);
          
          
          if iCondParam == 3

            for iCnd = 1:nCnd
              col = clrGoCorr{iCnd};
              line([xDataLim(1),xDataLim(2)],[meanZcGo(iCnd),meanZcGo(iCnd)], 'Color',col,'LineWidth',0.5);
            end

          else

            if normalizeDyn
              line([xDataLim(1),xDataLim(2)],[1 1], 'Color','k','LineWidth',0.5);
            else
              line([xDataLim(1),xDataLim(2)],[meanZcGo(1) meanZcGo(1)], 'Color','k','LineWidth',0.5);
            end

          end
          
          % Tag the panel
          text(xDataLim(2),0.2,figTag,'FontSize',8,'HorizontalAlignment','right');

        case 'indiv'
          
          if iCondParam == 2
          
            hFig = set_figure({20,20/1.61,'centimeters'},{'USLetter','landscape'});
            set(gcf,'PaperPositionMode','auto')
            pnl = panel;
            pnl.margin = [0 0 0 0];
            hold on;

            iCnd = 2;

            % Plot data
            col = clrGoCorrTransp;
            ln = lnGoCorr;
  %           cellfun(@(a,b) patchline(a,b,'EdgeColor',col,'LineWidth',1,'LineStyle',ln,'EdgeAlpha',edgeAlpha),sXGOT,sYGOT,'Uni',0);
            cellfun(@(a,b) plot(a,b,'Color',col,'LineWidth',1,'LineStyle',ln),sXGOT,sYGOT,'Uni',0);

            col = clrGoCommTransp;
            ln = lnGoCommS;
  %           cellfun(@(a,b) patchline(a,b,'EdgeColor',col,'LineWidth',1,'LineStyle',ln,'EdgeAlpha',edgeAlpha),sXGONTE,sYGONTE,'Uni',0);
            cellfun(@(a,b) plot(a,b,'Color',col,'LineWidth',1,'LineStyle',ln),sXGONTE,sYGONTE,'Uni',0);

            col = clrGoCorr;
            ln = lnGoCorr;
            plot(qXGOT,qYGOT,'Color',col,'LineWidth',3,'LineStyle',ln);

            col = clrGoComm;
            ln = lnGoCommQ;
            plot(qXGONTE,qYGONTE,'Color',col,'LineWidth',3,'LineStyle',ln);

           % Plot the threshold
            if normalizeDyn
              line([xDataLim(1),xDataLim(2)],[1 1],'Color','k','LineWidth',3);
            else
              thisThres = modelMat{iExSubj}.ZC{iCnd}(1); %Go threshold
              line([xDataLim(1),xDataLim(2)],[thisThres thisThres],'Color','k','LineWidth',3);
            end

            set(gca,'LineWidth',2);
            
            % File name
            if normalizeDyn
              normStr = 'normThres';
            else
              normStr = 'absThres';
            end
            fName = sprintf('DynamicsErrorChoice_c%s_i%s_condParam%s_subj%.2d_cond%.2d_%s.eps',choiceMechType,inhibMechType,condParam{iCnd},iExSubj,iCnd,normStr);
            
          end
          
      end
  end
      
  % Figure layout
  set(gca,'PlotBoxAspectRatio',[1.61,1,1],'XLim',xDataLim,'XTick',[],'YLim',yDataLimDyn,'YTick',[]);
%   ylabel('Activation (a.u.)');
%   title('Dynamics on error no-signal trials')
  
  if iCondParam == 2
    switch dynType
      case 'indiv'
        print('-depsc',fullfile(figDir,fName));
    end
  end

  % Panel 7 - Model dynamics: signal-inhibit trials
  % =========================================================================

  fprintf('Panel 7 - Model dynamics: signal-inhibit trials \n');
  
  % Pre-allocate matrices for logging
  % -------------------------------------------------------------------------
  qXSTOP = nan(nCnd,nQPrd,nSubj);
  qYSTOP = nan(nCnd,nQPrd,nSubj);
  qXGOT = nan(nCnd,nQPrd,nSubj);
  qYGOT = nan(nCnd,nQPrd,nSubj);

  % Extract, quantile average, and normalize predicted dynamics
  % -------------------------------------------------------------------------
  for iCnd = 1:nCnd
    
    switch dynType
      case 'qaverage'
        for iS = 1:nSubj

          % Subject index
          iSubj = subj(iS);

          % Get quantile averaged dynamics
          xI = prd{iSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.STOP.qX;
          yI = prd{iSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.STOP.qY;
          xO = quantile(prd{iSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.STOP.qX,qntlsPrd);
          yO = interp1(xI,yI,xO);
          
          qXSTOP(iCnd,:,iS) = xO;
          qYSTOP(iCnd,:,iS) = yO;
          
          xI = prd{iSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.GOT.qX;
          yI = prd{iSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.GOT.qY;
          xO = quantile(prd{iSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.GOT.qX,qntlsPrd);
          yO = interp1(xI,yI,xO);
          
          qXGOT(iCnd,:,iS) = xO;
          qYGOT(iCnd,:,iS) = yO;
          
          % Normalize dynamics
          if normalizeDyn
            allZc = cell2mat(modelMat{iSubj}.ZC(:)');
            zcGo = max(max(allZc(1:6,:)));
            zcStop = max(max(allZc(7,:)));
            qYSTOP(iCnd,:,iS) = qYSTOP(iCnd,:,iS)./zcStop;
            qYGOT(iCnd,:,iS) = qYGOT(iCnd,:,iS)./zcGo;
          end

        end

        meanQXSTOP = nanmean(qXSTOP,3);
        meanQYSTOP = nanmean(qYSTOP,3);
        meanQXGOT = nanmean(qXGOT,3);
        meanQYGOT = nanmean(qYGOT,3);
        
      case 'indiv' 
        
        if iCnd == 2
        
          % Get individual trial dynamics
          sXSTOP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.STOP.sX;
          sYSTOP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.STOP.sY;
          sXGOT = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.GOT.sX;
          sYGOT = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.GOT.sY;
          
          % Get quantile averaged dynamics
          qXSTOP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.STOP.qX;
          qYSTOP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.STOP.qY;
          qXGOT = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.GOT.qX;
          qYGOT = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopSuccess.goStim.GOT.qY;
          
          % Normalize dynamics
          if normalizeDyn
            zcStop = modelMat{iExSubj}.ZC{iCnd}(7);
            zcGo = modelMat{iExSubj}.ZC{iCnd}(3);

            sYSTOP = cellfun(@(a) a./zcStop,sYSTOP,'Uni',0);
            sYGOT = cellfun(@(a) a./zcGo,sYGOT,'Uni',0);

            qYSTOP = qYSTOP./zcStop;
            qYGOT = qYGOT./zcGo;
          end
          
        end
        
    end
  end
    
  % Plot predictions
  % -------------------------------------------------------------------------
  switch plotType
    case 'indiv'

      for iS = 1:nSubj

        figure; hold on;
        set(gcf,'Name',sprintf('c%s, i%s, p%s, subject %d',choiceMechType,inhibMechType,condParam{iCondParam},subj(iS)));

        for iCnd = 1:nCnd

          % Plot dynamics STOP
          col = clrGoCorr{iCnd};
          ln = lnGoCorr{iCnd};
          plot(qXSTOP(iCnd,:,iS),qYSTOP(iCnd,:,iS),'Color',col,'LineStyle',ln,'LineWidth',0.5);

          % Plot dynamics target GO
          col = clrGoComm{iCnd};
          ln = lnGoCommQ{iCnd};
          plot(qXGOT(iCnd,:,iS),qYGOT(iCnd,:,iS),'Color',col,'LineStyle',ln,'LineWidth',0.5);

        end
      end

    case 'qaverage'
      
      switch dynType
        
        case 'qaverage'
          
          % Select panel
          % -------------------------------------------------------------------
          p(7,iCondParam).select();
          p(7,iCondParam).hold('on');

          for iCnd = 1:3

              % Plot dynamics STOP
              col = clrGoCorr{iCnd};
              ln = lnGoCorr{iCnd};
              plot(meanQXSTOP(iCnd,:),meanQYSTOP(iCnd,:),'Color',col,'LineStyle',ln,'LineWidth',0.5);

              % Plot dynamics target GO
              col = clrGoComm{iCnd};
              ln = lnGoCommQ{iCnd};
              plot(meanQXGOT(iCnd,:),meanQYGOT(iCnd,:),'Color',col,'LineStyle',ln,'LineWidth',0.5);
              
              % Plot mean SSD
              theseSsds = reshape(ssd(iCnd,iStopFail-2,:),1,nSubj);
              
              line([mean(theseSsds),mean(theseSsds)],[0 100],'Color','k','LineStyle',':');

              % Plot mean SSRT
              thesePStopFails = num2cell(reshape(inhibFunObs(iCnd,iStopFail-2,:),1,nSubj));
              theseFinishTimes = cell2mat(cellfun(@(a,b) quantile(a.rt{iCnd,iGoCorr},b),obsData(8:13),thesePStopFails,'Uni',0));
                            
              line([mean(theseFinishTimes),mean(theseFinishTimes)],[0 100],'Color','k','LineStyle','--');
              
          end
          
          % Plot threshold
          zCGoAllSubj = nan(nSubj,nCnd);
          zCStopAllSubj = nan(nSubj,nCnd);
          for iS = 1:nSubj

            % Subject index
            iSubj = subj(iS);
            allZc = cell2mat(modelMat{iSubj}.ZC(:)');

            if normalizeDyn
              zcGo = max(max(allZc(1:6,:)));
              zCGoAllSubj(iS,:) = allZc(1,:)./zcGo;
            else
              zCGoAllSubj(iS,:) = allZc(1,:);
              zCStopAllSubj(iS,:) = allZc(7,:);
            end
          end  

          meanZcGo = nanmean(zCGoAllSubj,1);
          meanZcStop = nanmean(zCStopAllSubj,1);
          
          if iCondParam == 3

            for iCnd = 1:nCnd
              col = clrGoCorr{iCnd};
              line([xDataLim(1),xDataLim(2)],[meanZcGo(iCnd),meanZcGo(iCnd)], 'Color',col,'LineWidth',0.5);
              line([xDataLim(1),xDataLim(2)],[meanZcStop(iCnd),meanZcStop(iCnd)], 'Color',col,'LineWidth',0.5);
            end

          else
            
            if normalizeDyn
              line([xDataLim(1),xDataLim(2)],[1 1], 'Color','k','LineWidth',0.5);
            else
              line([xDataLim(1),xDataLim(2)],[meanZcGo(1) meanZcGo(1)], 'Color','k','LineWidth',0.5);
              line([xDataLim(1),xDataLim(2)],[meanZcStop(1) meanZcStop(1)], 'Color','k','LineWidth',0.5);
            end

          end

          % Tag the panel
          text(xDataLim(2),0.2,figTag,'FontSize',8,'HorizontalAlignment','right');
          
        case 'indiv'
            
            if iCondParam == 2
              hFig = set_figure({20,20/1.61,'centimeters'},{'USLetter','landscape'});
              set(gcf,'PaperPositionMode','auto')
              pnl = panel;
              pnl.margin = [0 0 0 0];
              hold on;

              iCnd = 2;

              % Plot data
              col = clrStopTransp;
              ln = lnStop;
  %             cellfun(@(a,b) patchline(a,b,'EdgeColor',col,'LineWidth',1,'LineStyle',ln,'EdgeAlpha',edgeAlpha),sXSTOP,sYSTOP,'Uni',0);
              cellfun(@(a,b) plot(a,b,'Color',col,'LineWidth',1,'LineStyle',ln),sXSTOP,sYSTOP,'Uni',0);

              col = clrGoCorrTransp;
              ln = lnGoCorr;
  %             cellfun(@(a,b) patchline(a,b,'EdgeColor',col,'LineWidth',1,'LineStyle',ln,'EdgeAlpha',edgeAlpha),sXGOT,sYGOT,'Uni',0);
              cellfun(@(a,b) plot(a,b,'Color',col,'LineWidth',1,'LineStyle',ln),sXGOT,sYGOT,'Uni',0);


              col = clrStop;
              ln = lnStop;
              plot(qXSTOP,qYSTOP,'Color',col,'LineWidth',3,'LineStyle',ln);

              col = clrGoCorr;
              ln = lnGoCorr;
              plot(qXGOT,qYGOT,'Color',col,'LineWidth',3,'LineStyle',ln);

              % Plot the threshold
              if normalizeDyn
                line([xDataLim(1),xDataLim(2)],[1 1],'Color','k','LineWidth',3);
              else
                thisGoThres = modelMat{iExSubj}.ZC{iCnd}(1); %Go threshold
                thisStopThres = modelMat{iExSubj}.ZC{iCnd}(7); %Go threshold
                line([xDataLim(1),xDataLim(2)],[thisGoThres thisGoThres],'Color','k','LineWidth',3);
                line([xDataLim(1),xDataLim(2)],[thisStopThres thisStopThres],'Color','k','LineWidth',3);
              end

              set(gca,'LineWidth',2);
              
              % Plot SSD & SSRT
              thisSsd = ssd(iCnd,iStopFail-2,find(subj == iExSubj));
              thisPStopFail = inhibFunObs(iCnd,iStopFail-2,find(subj == iExSubj));
              thisRt = obsData{iExSubj}.rt{iCnd,iGoCorr};
              thisFinishTime = quantile(thisRt,thisPStopFail);

              if normalizeDyn
                line([thisSsd,thisSsd],[0 1],'Color','k','LineWidth',1.5,'LineStyle','--');
                line([thisFinishTime,thisFinishTime],[0 1],'Color','k','LineWidth',1.5,'LineStyle','--');
              else
                line([thisSsd,thisSsd],[0 thisGoThres],'Color','k','LineWidth',1.5,'LineStyle','--');
                line([thisFinishTime,thisFinishTime],[0 thisGoThres],'Color','k','LineWidth',1.5,'LineStyle','--');
              end
              
              % File name
              if normalizeDyn
                normStr = 'normThres';
              else
                normStr = 'absThres';
              end
              fName = sprintf('DynamicsStopSuccess_c%s_i%s_condParam%s_subj%.2d_cond%.2d_%s.eps',choiceMechType,inhibMechType,condParam{iCnd},iExSubj,iCnd,normStr);
              
            end
      end
  end

  % Figure layout
  set(gca,'PlotBoxAspectRatio',[1.61,1,1],'XLim',xDataLim,'XTick',[],'YLim',yDataLimDyn,'YTick',[]);
%   ylabel('Activation (a.u.)');
%   title('Dynamics on signal-respond trials');

  if iCondParam == 2
    switch dynType
      case 'indiv'
        print('-depsc',fullfile(figDir,fName));
    end
  end

  % Panel 8 - Model dynamics: signal-respond trials
  % =========================================================================

  fprintf('Panel 8 - Model dynamics: signal-respond trials \n');
  
  % Pre-allocate matrices for logging
  % -------------------------------------------------------------------------
  qXSTOP = nan(nCnd,nQPrd,nSubj);
  qYSTOP = nan(nCnd,nQPrd,nSubj);
  qXGORESP = nan(nCnd,nQPrd,nSubj);
  qYGORESP = nan(nCnd,nQPrd,nSubj);

  % Extract, quantile average, and normalize predicted dynamics
  % -------------------------------------------------------------------------
  for iCnd = 1:nCnd
    
    switch dynType
      case 'qaverage'
    
        for iS = 1:nSubj

          % Subject index
          iSubj = subj(iS);

          % Get quantile averaged dynamics
          xI = prd{iSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.STOP.qX;
          yI = prd{iSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.STOP.qY;
          xO = quantile(prd{iSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.STOP.qX,qntlsPrd);
          yO = interp1(xI,yI,xO);
          
          qXSTOP(iCnd,:,iS) = xO;
          qYSTOP(iCnd,:,iS) = yO;
          
          xI = prd{iSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.GORESP.qX;
          yI = prd{iSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.GORESP.qY;
          xO = quantile(prd{iSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.GORESP.qX,qntlsPrd);
          yO = interp1(xI,yI,xO);
          
          qXGORESP(iCnd,:,iS) = xO;
          qYGORESP(iCnd,:,iS) = yO;
          
          % Normalize dynamics
          if normalizeDyn
            allZc = cell2mat(modelMat{iSubj}.ZC(:)');
            zcGo = max(max(allZc(1:6,:)));
            zcStop = max(max(allZc(7,:)));
            qYSTOP(iCnd,:,iS) = qYSTOP(iCnd,:,iS)./zcStop;
            qYGORESP(iCnd,:,iS) = qYGORESP(iCnd,:,iS)./zcGo;
          end

        end
 
        meanQXSTOP = nanmean(qXSTOP,3);
        meanQYSTOP = nanmean(qYSTOP,3);
        meanQXGORESP = nanmean(qXGORESP,3);
        meanQYGORESP = nanmean(qYGORESP,3);

      case 'indiv'
        
        if iCnd == 2
        
          % Get individual trial dynamics
          sXSTOP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.STOP.sX;
          sYSTOP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.STOP.sY;
          sXGORESP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.GORESP.sX;
          sYGORESP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.GORESP.sY;
          
          % Get quantile averaged dynamics
          qXSTOP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.STOP.qX;
          qYSTOP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.STOP.qY;
          qXGORESP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.GORESP.qX;
          qYGORESP = prd{iExSubj}.dyn{iCnd,iStopFail-1}.StopFailure.goStim.GORESP.qY;
          
          % Normalize dynamics
          if normalizeDyn
            zcStop = modelMat{iExSubj}.ZC{iCnd}(7);
            zcGo = modelMat{iExSubj}.ZC{iCnd}(3);

            sYSTOP = cellfun(@(a) a./zcStop,sYSTOP,'Uni',0);
            sYGORESP = cellfun(@(a) a./zcGo,sYGORESP,'Uni',0);

            qYSTOP = qYSTOP./zcStop;
            qYGORESP = qYGORESP./zcGo;
          end
          
        end
        
    end
  end
        
        
  % Plot predictions
  % -------------------------------------------------------------------------
  switch plotType
    case 'indiv'

      for iS = 1:nSubj

        figure; hold on;
        set(gcf,'Name',sprintf('c%s, i%s, p%s, subject %d',choiceMechType,inhibMechType,condParam{iCondParam},subj(iS)));

        for iCnd = 1:nCnd

          % Plot dynamics STOP
          col = clrGoCorr{iCnd};
          ln = lnGoCorr{iCnd};
          plot(qXSTOP(iCnd,:,iS),qYSTOP(iCnd,:,iS),'Color',col,'LineStyle',ln,'LineWidth',0.5);

          % Plot dynamics target GORESP
          col = clrGoComm{iCnd};
          ln = lnGoCommQ{iCnd};
          plot(qXGORESP(iCnd,:,iS),qYGORESP(iCnd,:,iS),'Color',col,'LineStyle',ln,'LineWidth',0.5);

        end
      end

    case 'qaverage'
      
      switch dynType
        
        case 'qaverage'
      
          % Select panel
          % -------------------------------------------------------------------
          p(8,iCondParam).select();
          p(8,iCondParam).hold('on');

          for iCnd = 1:3

              % Plot dynamics STOP
              col = clrGoCorr{iCnd};
              ln = lnGoCorr{iCnd};
              plot(meanQXSTOP(iCnd,:),meanQYSTOP(iCnd,:),'Color',col,'LineStyle',ln,'LineWidth',0.5);

              % Plot dynamics target GO
              col = clrGoComm{iCnd};
              ln = lnGoCommQ{iCnd};
              plot(meanQXGORESP(iCnd,:),meanQYGORESP(iCnd,:),'Color',col,'LineStyle',ln,'LineWidth',0.5);

          end

          % Plot threshold
          zCGoAllSubj = nan(nSubj,nCnd);
          zCStopAllSubj = nan(nSubj,nCnd);
          for iS = 1:nSubj

            % Subject index
            iSubj = subj(iS);
            allZc = cell2mat(modelMat{iSubj}.ZC(:)');

            if normalizeDyn
              zcGo = max(max(allZc(1:6,:)));
              zCGoAllSubj(iS,:) = allZc(1,:)./zcGo;
            else
              zCGoAllSubj(iS,:) = allZc(1,:);
              zCStopAllSubj(iS,:) = allZc(7,:);
            end
          end  

          meanZcGo = nanmean(zCGoAllSubj,1);
          meanZcStop = nanmean(zCStopAllSubj,1);
          
          if iCondParam == 3

            for iCnd = 1:nCnd
              col = clrGoCorr{iCnd};
              line([xDataLim(1),xDataLim(2)],[meanZcGo(iCnd),meanZcGo(iCnd)], 'Color',col,'LineWidth',0.5);
              line([xDataLim(1),xDataLim(2)],[meanZcStop(iCnd),meanZcStop(iCnd)], 'Color',col,'LineWidth',0.5);
            end

          else
            
            if normalizeDyn
              line([xDataLim(1),xDataLim(2)],[1 1], 'Color','k','LineWidth',0.5);
            else
              line([xDataLim(1),xDataLim(2)],[meanZcGo(1) meanZcGo(1)], 'Color','k','LineWidth',0.5);
              line([xDataLim(1),xDataLim(2)],[meanZcStop(1) meanZcStop(1)], 'Color','k','LineWidth',0.5);
            end

          end
          
          % Tag the panel
          text(xDataLim(2),0.2,figTag,'FontSize',8,'HorizontalAlignment','right');
 
        case 'indiv'
          
          if iCondParam == 2
          
            hFig = set_figure({20,20/1.61,'centimeters'},{'USLetter','landscape'});
            set(gcf,'PaperPositionMode','auto')
            pnl = panel;
            pnl.margin = [0 0 0 0];
            hold on;

            iCnd = 2;

            % Plot data
            col = clrStopTransp;
            ln = lnStop;
  %           cellfun(@(a,b) patchline(a,b,'EdgeColor',col,'LineWidth',1,'LineStyle',ln,'EdgeAlpha',edgeAlpha),sXSTOP,sYSTOP,'Uni',0);
            cellfun(@(a,b) plot(a,b,'Color',col,'LineWidth',1,'LineStyle',ln),sXSTOP,sYSTOP,'Uni',0);

            col = clrGoCorrTransp;
            ln = lnGoCorr;
  %           cellfun(@(a,b) patchline(a,b,'EdgeColor',col,'LineWidth',1,'LineStyle',ln,'EdgeAlpha',edgeAlpha),sXGORESP,sYGORESP,'Uni',0);
            cellfun(@(a,b) plot(a,b,'Color',col,'LineWidth',1,'LineStyle',ln),sXGORESP,sYGORESP,'Uni',0);

            col = clrStop;
            ln = lnStop;
            plot(qXSTOP,qYSTOP,'Color',col,'LineWidth',3,'LineStyle',ln);

            col = clrGoCorr;
            ln = lnGoCorr;
            plot(qXGORESP,qYGORESP,'Color',col,'LineWidth',3,'LineStyle',ln);

            % Plot the threshold
            if normalizeDyn
              line([xDataLim(1),xDataLim(2)],[1 1],'Color','k','LineWidth',3);
            else
              thisGoThres = modelMat{iExSubj}.ZC{iCnd}(1); %Go threshold
              thisStopThres = modelMat{iExSubj}.ZC{iCnd}(7); %Go threshold
              line([xDataLim(1),xDataLim(2)],[thisGoThres thisGoThres],'Color','k','LineWidth',3);
              line([xDataLim(1),xDataLim(2)],[thisStopThres thisStopThres],'Color','k','LineWidth',3);
            end
            
            set(gca,'LineWidth',2);
            
            % Plot SSD & SSRT
            thisSsd = ssd(iCnd,iStopFail-2,find(subj == iExSubj));
            thisPStopFail = inhibFunObs(iCnd,iStopFail-2,find(subj == iExSubj));
            thisRt = obsData{iExSubj}.rt{iCnd,iGoCorr};
            thisFinishTime = quantile(thisRt,thisPStopFail);
            
            if normalizeDyn
              line([thisSsd,thisSsd],[0 1],'Color','k','LineWidth',1.5,'LineStyle','--');
              line([thisFinishTime,thisFinishTime],[0 1],'Color','k','LineWidth',1.5,'LineStyle','--');
            else
              line([thisSsd,thisSsd],[0 thisGoThres],'Color','k','LineWidth',1.5,'LineStyle','--');
              line([thisFinishTime,thisFinishTime],[0 thisGoThres],'Color','k','LineWidth',1.5,'LineStyle','--');
            end
               
            % File name
            if normalizeDyn
              normStr = 'normThres';
            else
              normStr = 'absThres';
            end
            fName = sprintf('DynamicsStopFailure_c%s_i%s_condParam%s_subj%.2d_cond%.2d_%s.eps',choiceMechType,inhibMechType,condParam{iCnd},iExSubj,iCnd,normStr);
                      
          end
      end    
  end

   
  % Figure layout
  set(gca,'PlotBoxAspectRatio',[1.61,1,1],'XLim',xDataLim,'XTick',[],'YLim',yDataLimDyn,'YTick',[]);
%   ylabel('Activation (a.u.)');
%   title('Dynamics on signal-respond trials');

  if iCondParam == 2
    switch dynType
      case 'indiv'
        print('-depsc',fullfile(figDir,fName));
    end
  end

end
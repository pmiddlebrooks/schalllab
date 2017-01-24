%%
iCnd = 6;
iGOT = 1;
iSTOP = 3;



%%

figureHandle = 650;
nColumn     = 3;
nRow        = 2;   % for now, plot 2 rows-- need to add more with stop trials
rowGo     = 1;
rowStop     = 2;
rowInh      = 3;
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
clf


% Go trials: Correct and Error
  ax(rowGo, 1) = axes('units', 'centimeters', 'position', [xAxesPosition(rowGo, 1) yAxesPosition(rowGo, 1) axisWidth axisHeight]);
  cla
  set(ax(rowGo, 1), 'ylim', [0 1], 'xlim', [0 xMax])
  hold(ax(rowGo, 1), 'on')
  %             ttl = sprintf('SSD:  %d', iSSD);
  %             title(ttl)
  
  if includeStop
    
    % Stop Trials: Correct and Error
    ax(rowStop, 1) = axes('units', 'centimeters', 'position', [xAxesPosition(rowStop, 1) yAxesPosition(rowStop, 1) axisWidth axisHeight]);
    cla
    set(ax(rowStop, 1), 'ylim', [0 1], 'xlim', [0 xMax])
    hold(ax(rowStop, 1), 'on')
  end


% EdgeAlpha level:
edgeAlpha = 0.05;

% Colors
colGOT = [0 0.5 0];
colGONT = [0.5 0.5 0];
colGONTE = [0.5 0.5 0];
colSTOP = [1 0 0];

% Settings of axes
xStimLim = [250 2250];
xDataLim = [0 2000];


simScope = SAM.sim.scope;
simScope = 'go';
switch lower(simScope)
  case 'go'
    N = 2;
  case 'all'
    N = [2 1];
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GO TRIALS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iGoTrType = 1;


% Get RT data
rtGoCorr = sort(prd.rtGoCorr{iCnd});
rtGoComm = sort(prd.rtGoComm{iCnd});

% Cumulative probabilities
FGoCorr = mtb_edf(rtGoCorr(:),rtGoCorr(:));
FGoComm = mtb_edf(rtGoComm(:),rtGoComm(:));

% Plot predicted data
plot(ax(rowGo, 1), rtGoCorr,FGoCorr,'Color',colGOT,'LineWidth',3);
plot(ax(rowGo, 1), FGoComm,'Color',colGONTE,'LineWidth',3);

% Plot observed quantiles
% scatter(quantile(obs.rtGoCorr{iCnd},[.1 .3 .5 .7 .9]),[.1 .3 .5 .7 .9],50,'ko','MarkerEdgeColor',colGOT,'LineWidth',2);
% scatter(quantile(obs.rtGoComm{iCnd},[.1 .3 .5 .7 .9]),[.1 .3 .5 .7 .9],50,'ko','MarkerEdgeColor',colGONTE,'LineWidth',2);

% Adjust axes
set(gca,'XLim',xDataLim, ...
         'YTick',[0 1]);
ylabel('Cumulative probability');
title('RTs and probabilities on no-signal trials');
axes('Position',[0.4 0.6 0.1 0.1])


clear rtGoCorr rtGoComm FGoCorr FGoComm



% Dynamics Go correct trial
% =========================================================================

if prd.pGoCorr(iCnd) > 0


  % Get individual trial dynamics
  sXGOT = prd.dyn{iCnd,iGoTrType}.GoCorr.goStim.GOT.sX;
  sYGOT = prd.dyn{iCnd,iGoTrType}.GoCorr.goStim.GOT.sY;
  sXGONT = prd.dyn{iCnd,iGoTrType}.GoCorr.goStim.GONT.sX;
  sYGONT = prd.dyn{iCnd,iGoTrType}.GoCorr.goStim.GONT.sY;

  % Get quantile averaged dynamics
  qXGOT = prd.dyn{iCnd,iGoTrType}.GoCorr.goStim.GOT.qX;
  qYGOT = prd.dyn{iCnd,iGoTrType}.GoCorr.goStim.GOT.qY;
  qXGONT = prd.dyn{iCnd,iGoTrType}.GoCorr.goStim.GONT.qX;
  qYGONT = prd.dyn{iCnd,iGoTrType}.GoCorr.goStim.GONT.qY;

  % Plot data
  cellfun(@(a,b) patchline(a,b,'EdgeColor',colGOT,'LineWidth',1,'EdgeAlpha',edgeAlpha),sXGOT,sYGOT,'Uni',0);
  cellfun(@(a,b) patchline(a,b,'EdgeColor',colGONT,'LineWidth',1,'EdgeAlpha',edgeAlpha),sXGONT,sYGONT,'Uni',0);

  plot(ax(rowGo, 1), qXGOT,qYGOT,'Color',colGOT,'LineWidth',3);
  plot(ax(rowGo, 1), qXGONT,qYGONT,'Color',colGONT,'LineWidth',3);

  line([xDataLim(1),xDataLim(2)], ...
       [modelMat.ZC{iCnd}(iGOT),modelMat.ZC{iCnd}(iGOT)], 'Color','k','LineWidth',3);

  % Adjust axes
  set(gca,'XLim',xDataLim);
  ylabel('Activation (a.u.)');
  title('Dynamics on correct no-signal trials')

  clear sXGOT sYGOT sXGONT sYGONT qXGOT qYGOT qXGONT qYGONT
  
end

% Dynamics Go commission error trial
% =========================================================================

if prd.pGoComm(iCnd) > 0


  % Get individual trial dynamics
  sXGOT = prd.dyn{iCnd,iGoTrType}.GoComm.goStim.GOT.sX;
  sYGOT = prd.dyn{iCnd,iGoTrType}.GoComm.goStim.GOT.sY;
  sXGONTE = prd.dyn{iCnd,iGoTrType}.GoComm.goStim.GONTE.sX;
  sYGONTE = prd.dyn{iCnd,iGoTrType}.GoComm.goStim.GONTE.sY;

  % Get quantile averaged dynamics
  qXGOT = prd.dyn{iCnd,iGoTrType}.GoComm.goStim.GOT.qX;
  qYGOT = prd.dyn{iCnd,iGoTrType}.GoComm.goStim.GOT.qY;
  qXGONTE = prd.dyn{iCnd,iGoTrType}.GoComm.goStim.GONTE.qX;
  qYGONTE = prd.dyn{iCnd,iGoTrType}.GoComm.goStim.GONTE.qY;

  % Plot data
  cellfun(@(a,b) patchline(a,b,'EdgeColor',colGOT,'LineWidth',1,'EdgeAlpha',edgeAlpha),sXGOT,sYGOT,'Uni',0);
  cellfun(@(a,b) patchline(a,b,'EdgeColor',colGONTE,'LineWidth',1,'EdgeAlpha',edgeAlpha),sXGONTE,sYGONTE,'Uni',0);

  plot(ax(rowGo, 1), qXGOT,qYGOT,'Color',colGOT,'LineWidth',3);
  plot(ax(rowGo, 1), qXGONTE,qYGONTE,'Color',colGONTE,'LineWidth',3);

  line([xDataLim(1),xDataLim(2)], ...
       [modelMat.ZC{iCnd}(iGOT),modelMat.ZC{iCnd}(iGOT)], 'Color','k','LineWidth',3);

  % Adjust axes
  set(gca,'XLim',xDataLim);
  xlabel('Time from go-signal (ms');
  ylabel('Activation (a.u.)');
  title('Dynamics on choice error no-signal trials');

  clear sXGOT sYGOT sXGONTE sYGONTE qXGOT qYGOT qXGONTE qYGONTE
  
end


%%

iGoTrType = 1;


% Get RT data
rtGoCorr = sort(prd.rtGoCorr{iCnd});
rtGoComm = sort(prd.rtGoComm{iCnd});

% Cumulative probabilities
FGoCorr = mtb_edf(rtGoCorr(:),rtGoCorr(:));
FGoComm = mtb_edf(rtGoComm(:),rtGoComm(:));


figure;hold on;
plot(rtGoCorr,FGoCorr,'Color',colGOT,'LineWidth',3);
plot(FGoComm,'Color',colGONTE,'LineWidth',3);








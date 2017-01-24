function sam_plot_obs_prd(SAM,obsOptimData,prdOptimData,obs,prd)
% SAM_PLOT_OBS_PRD <Synopsis of what this function does> 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% SAM_PLOT_OBS_PRD; 
%  
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Mon 30 Sep 2013 09:22:57 CDT by bram 
% $Modified: Mon 30 Sep 2013 09:22:57 CDT by bram 

 
% CONTENTS 

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% #.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 
% #.#.
% ========================================================================= 

 
% #.#.#.
% ------------------------------------------------------------------------- 

nSim = SAM.sim.nSim;


switch size(obsOptimData.rt,2)
  case 2

    qntls = cell(3,2);

    clrGo = cell(3,2);
    clrGo{1,1} = [0 0.5 1];
    clrGo{1,2} = [0 0.5 1];
    clrGo{2,1} = [0 0.5 0];
    clrGo{2,2} = [0 0.5 0];
    clrGo{3,1} = [1 0 0];
    clrGo{3,2} = [1 0 0];

    lnGo = cell(3,2);
    lnGo{1,1} = '-';
    lnGo{1,2} = '--';
    lnGo{2,1} = '-';
    lnGo{2,2} = '--';
    lnGo{3,1} = '-';
    lnGo{3,2} = '--';

    mrkGo = cell(3,2);
    mrkGo{1,1} = 'o';
    mrkGo{1,2} = '^';
    mrkGo{2,1} = 'o';
    mrkGo{2,2} = '^';
    mrkGo{3,1} = 'o';
    mrkGo{3,2} = '^';
    
  case 7
    
    qntls = cell(3,7);

    clrStop = cell(3,7);
    clrStop{1,1} = [0 0.5 1];
    clrStop{1,2} = [0 0.5 1];
    clrStop{1,3} = [0 0.5 1];
    clrStop{1,4} = [0 0.5 1];
    clrStop{1,5} = [0 0.5 1];
    clrStop{1,6} = [0 0.5 1];
    clrStop{1,7} = [0 0.5 1];
        
    clrStop{2,1} = [0 0.5 0];
    clrStop{2,2} = [0 0.5 0];
    clrStop{2,3} = [0 0.5 0];
    clrStop{2,4} = [0 0.5 0];
    clrStop{2,5} = [0 0.5 0];
    clrStop{2,6} = [0 0.5 0];
    clrStop{2,7} = [0 0.5 0];
    
    clrStop{3,1} = [1 0 0];
    clrStop{3,2} = [1 0 0];
    clrStop{3,3} = [1 0 0];
    clrStop{3,4} = [1 0 0];
    clrStop{3,5} = [1 0 0];
    clrStop{3,6} = [1 0 0];
    clrStop{3,7} = [1 0 0];
    
    lnStop = cell(3,7);
    lnStop{1,1} = '-';
    lnStop{1,2} = '--';
    lnStop{1,3} = '-.';
    lnStop{1,4} = '-.';
    lnStop{1,5} = '-.';
    lnStop{1,6} = '-.';
    lnStop{1,7} = '-.';
    
    lnStop{2,1} = '-';
    lnStop{2,2} = '--';
    lnStop{2,3} = '-.';
    lnStop{2,4} = '-.';
    lnStop{2,5} = '-.';
    lnStop{2,6} = '-.';
    lnStop{2,7} = '-.';
    
    lnStop{3,1} = '-';
    lnStop{3,2} = '--';
    lnStop{3,3} = '-.';
    lnStop{3,4} = '-.';
    lnStop{3,5} = '-.';
    lnStop{3,6} = '-.';
    lnStop{3,7} = '-.';

    mrkStop = cell(3,7);
    mrkStop{1,1} = 'o';
    mrkStop{1,2} = '^';
    mrkStop{1,3} = 's';
    mrkStop{1,4} = 's';
    mrkStop{1,5} = 's';
    mrkStop{1,6} = 's';
    mrkStop{1,7} = 's';
    
    mrkStop{2,1} = 'o';
    mrkStop{2,2} = '^';
    mrkStop{2,3} = 's';
    mrkStop{2,4} = 's';
    mrkStop{2,5} = 's';
    mrkStop{2,6} = 's';
    mrkStop{2,7} = 's';
    
    mrkStop{3,1} = 'o';
    mrkStop{3,2} = '^';
    mrkStop{3,3} = 's';
    mrkStop{3,4} = 's';
    mrkStop{3,5} = 's';
    mrkStop{3,6} = 's';
    mrkStop{3,7} = 's';
    
    prdOptimDataForStop = prdOptimData;
    
end


% lnGo{2,1} = '-';
% lnGo{2,2} = '--';
% lnGo{3,1} = '-';
% lnGo{3,2} = '--';
% 
% lnGo{1,2} = '--';
% 
% mrkGo = cell(3,2);
% mrkGo{1,1} = 'o';
% mrkGo{1,2} = '^';
% mrkGo{2,1} = 'o';
% mrkGo{2,2} = '^';
% mrkGo{3,1} = 'o';
% mrkGo{3,2} = '^';

% Plot correct and choice error RTs
figure; hold on;

clr = cell(3,1);
clr{1,1} = [0 0.5 1];
clr{2,1} = [0 0.5 0];
clr{3,1} = [1 0 0];

ln = cell(3,1);
ln{1,1} = '-';
ln{2,1} = '-';
ln{3,1} = '-';

mrk = cell(3,1);
mrk{1,1} = 'o';
mrk{2,1} = 'o';
mrk{3,1} = 'o';

% Plot predicted correct Go RTs
cellfun(@(a,b,c) plot(a,mtb_edf(a,a),'Color',b,'LineStyle',c,'LineWidth',2),prd.rtGoCorr,clr,ln,'Uni',0);

% Plot observed correct Go RTs
cellfun(@(a,b,c,d,e) scatter(a,b./c,50,d,'Marker',e),obsOptimData.rtQ(:,1),obsOptimData.pDefect(:,1),num2cell(obsOptimData.P(:,1)),clr,mrk,'Uni',0)

clr       = cell(3,1);
clr{1,1}  = [0 0.5 1];
clr{2,1}  = [0 0.5 0];
clr{3,1}  = [1 0 0];

ln        = cell(3,1);
ln{1,1}   = '--';
ln{2,1}   = '--';
ln{3,1}   = '--';

mrk       = cell(3,1);
mrk{1,1}  = '^';
mrk{2,1}  = '^';
mrk{3,1}  = '^';

% Plot predicted error Go RTs
cellfun(@(a,b,c) plot(a,mtb_edf(a,a),'Color',b,'LineStyle',c,'LineWidth',2),prd.rtGoComm,clr,ln,'Uni',0);

% Plot observed error Go  RTs
cellfun(@(a,b,c,d,e) scatter(a,b./c,50,d,'Marker',e),obsOptimData.rtQ(:,2),obsOptimData.pDefect(:,2),num2cell(obsOptimData.P(:,2)),clr,mrk,'Uni',0);

% Set axes
xlabel('RT (ms)');
ylabel('Cumulative probability');

% Plot correct and inhibition error RTs
% =======================

figure; hold on;

clr = cell(3,1);
clr{1,1} = [0 0.5 1];
clr{2,1} = [0 0.5 0];
clr{3,1} = [1 0 0];

ln = cell(3,1);
ln{1,1} = '-';
ln{2,1} = '-';
ln{3,1} = '-';

mrk = cell(3,1);
mrk{1,1} = 'o';
mrk{2,1} = 'o';
mrk{3,1} = 'o';

% Plot predicted correct Go RTs
cellfun(@(a,b,c) plot(a,mtb_edf(a,a),'Color',b,'LineStyle',c,'LineWidth',2),prd.rtGoCorr,clr,ln,'Uni',0);

% Plot observed correct Go RTs
cellfun(@(a,b,c,d,e) scatter(a,b./c,50,d,'Marker',e),obsOptimData.rtQ(:,1),obsOptimData.pDefect(:,1),num2cell(obsOptimData.P(:,1)),clr,mrk,'Uni',0)

clr       = cell(3,1);
clr{1,1}  = [0 0.5 1];
clr{2,1}  = [0 0.5 0];
clr{3,1}  = [1 0 0];

ln        = cell(3,1);
ln{1,1}   = ':';
ln{2,1}   = ':';
ln{3,1}   = ':';

mrk       = cell(3,1);
mrk{1,1}  = 's';
mrk{2,1}  = 's';
mrk{3,1}  = 's';

iSsd = 3;

% Plot predicted StopFailure RT
cellfun(@(a,b,c) plot(a,mtb_edf(a,a),'Color',b,'LineStyle',c,'LineWidth',2),prd.rtStopFailure(:,iSsd),clr,ln,'Uni',0);

% Plot observed StopFailure RT
cellfun(@(a,b,c,d,e) scatter(a,b./c,50,d,'Marker',e),obsOptimData.rtQ(:,2+iSsd),obsOptimData.pDefect(:,2+iSsd),num2cell(obsOptimData.P(:,2+iSsd)),clr,mrk,'Uni',0)

xlabel('RT (ms)');
ylabel('Cumulative probability');
title('Signal-respond RTs for middle SSD in subject 08');



% Plot inhibition function
% =======================

figure; hold on;

clr = cell(3,1);
clr{1,1} = [0 0.5 1];
clr{2,1} = [0 0.5 0];
clr{3,1} = [1 0 0];

ln = cell(3,1);
ln{1,1} = '-';
ln{2,1} = '-';
ln{3,1} = '-';

mrk = cell(3,1);
mrk{1,1} = 'o';
mrk{2,1} = 'o';
mrk{3,1} = 'o';

ssdCell = mat2cell(obs.ssd,ones(3,1),5);
pStopFailureCell = mat2cell(obs.pStopFailure,ones(3,1),5);

% Plot predicted inhibition functions
cellfun(@(a,b,c,d) plot(a,b,'Color',c,'LineStyle',d,'LineWidth',2),ssdCell,prd.inhibFunc,clr,ln,'Uni',0);

% Plot observed inhibition functions
cellfun(@(a,b,c,d) scatter(a,b,50,c,'Marker',d),ssdCell,pStopFailureCell,clr,mrk,'Uni',0);






















% Compute how many NaNs to add
nansToAdd = num2cell(nSim - cellfun(@(a) numel(a),prdOptimData.rt));

% Plot defective data, y-axis [0,1]
subplot(1,2,1);
cla;hold on;

prdOptimData.rt = cellfun(@(a,b) [a(:);nan(b,1)],prdOptimData.rt,nansToAdd,'Uni',0);

% Plot predictions
cellfun(@(a,b,c) plot(a,mtb_edf(a,a),'Color',b,'LineStyle',c,'LineWidth',2),prdOptimData.rt,clr,ln,'Uni',0);

% Plot observations
cellfun(@(a,b,c,d) scatter(a,b,50,c,'Marker',d),obsOptimData.rtQ,obsOptimData.pDefect,clr,mrk,'Uni',0)

% Set axes
set(gca,'PlotBoxAspectRatio',[1.6 1 1],'XLim',[0 1200],'YLim',[0 1]);

% % Sort data, and make sure it's a column vector
% prdOptimData.rt = cellfun(@(a) sort(a(:)),prdOptimData.rt,'Uni',0);
% 
% % Plot predictions
% cellfun(@(a,b,c) plot(a,mtb_edf(a,a),'Color',b,'LineStyle',c,'LineWidth',2),prdOptimData.rt,clr,ln,'Uni',0);
% 
% % Plot observations
% cellfun(@(a,b,c,d,e) scatter(a,b./c,20,d,'Marker',e),obsOptimData.rtQ,obsOptimData.pDefect,num2cell(obsOptimData.P),clr,mrk,'Uni',0);

% Plot defective data, y-axis [0,0.1]
subplot(1,2,2);
cla;hold on;

prdOptimData.rt = cellfun(@(a,b) [a(:);nan(b,1)],prdOptimData.rt,nansToAdd,'Uni',0);

% Plot predictions
cellfun(@(a,b,c) plot(a,mtb_edf(a,a),'Color',b,'LineStyle',c,'LineWidth',2),prdOptimData.rt,clr,ln,'Uni',0);

% Plot observations
cellfun(@(a,b,c,d) scatter(a,b,50,c,'Marker',d),obsOptimData.rtQ,obsOptimData.pDefect,clr,mrk,'Uni',0)

% Set axes
set(gca,'PlotBoxAspectRatio',[1.6 1 1],'XLim',[0 1200],'YLim',[0 0.1]);


switch size(obsOptimData.rt,2)
  case 7
    
    figure;
    
    for iCnd = 1:3
    
      subplot(1,3,iCnd);
      cla;hold on;

      iStopTr = 3:7;
      
      % Plot predictions
      cellfun(@(a,b,c) plot(a,mtb_edf(a,a),'Color',b,'LineStyle',c,'LineWidth',2),prdOptimDataForStop.rt(iCnd,[1,iStopTr]),clr(iCnd,[1,iStopTr]),ln(iCnd,[1,iStopTr]),'Uni',0);

      % Plot observations
      cellfun(@(a,b,c,d,e) scatter(a,b./c,50,d,'Marker',e),obsOptimData.rtQ(iCnd,[1,iStopTr]),obsOptimData.pDefect(iCnd,[1,iStopTr]),num2cell(obsOptimData.P(iCnd,[1,iStopTr])),clr(iCnd,[1,iStopTr]),mrk(iCnd,[1,iStopTr]),'Uni',0)

      % Set axes
      set(gca,'PlotBoxAspectRatio',[1.6 1 1],'XLim',[0 1200],'YLim',[0 1]);
      
    end
    
    
end



% Draw now
drawnow;
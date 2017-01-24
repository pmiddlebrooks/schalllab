function sam_plot_obs_prd(SAM,obsOptimData,prdOptimData)
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

    clr = cell(3,2);
    clr{1,1} = [0 0.5 1];
    clr{1,2} = [0 0.5 1];
    clr{2,1} = [0 0.5 0];
    clr{2,2} = [0 0.5 0];
    clr{3,1} = [1 0 0];
    clr{3,2} = [1 0 0];

    ln = cell(3,2);
    ln{1,1} = '-';
    ln{1,2} = '--';
    ln{2,1} = '-';
    ln{2,2} = '--';
    ln{3,1} = '-';
    ln{3,2} = '--';

    mrk = cell(3,2);
    mrk{1,1} = 'o';
    mrk{1,2} = '^';
    mrk{2,1} = 'o';
    mrk{2,2} = '^';
    mrk{3,1} = 'o';
    mrk{3,2} = '^';
    
  case 7
    
    qntls = cell(3,7);

    clr = cell(3,7);
    clr{1,1} = [0 0.5 1];
    clr{1,2} = [0 0.5 1];
    clr{1,3} = [0 0.5 1];
    clr{1,4} = [0 0.5 1];
    clr{1,5} = [0 0.5 1];
    clr{1,6} = [0 0.5 1];
    clr{1,7} = [0 0.5 1];
        
    clr{2,1} = [0 0.5 0];
    clr{2,2} = [0 0.5 0];
    clr{2,3} = [0 0.5 0];
    clr{2,4} = [0 0.5 0];
    clr{2,5} = [0 0.5 0];
    clr{2,6} = [0 0.5 0];
    clr{2,7} = [0 0.5 0];
    
    clr{3,1} = [1 0 0];
    clr{3,2} = [1 0 0];
    clr{3,3} = [1 0 0];
    clr{3,4} = [1 0 0];
    clr{3,5} = [1 0 0];
    clr{3,6} = [1 0 0];
    clr{3,7} = [1 0 0];
    
    ln = cell(3,7);
    ln{1,1} = '-';
    ln{1,2} = '--';
    ln{1,3} = '-.';
    ln{1,4} = '-.';
    ln{1,5} = '-.';
    ln{1,6} = '-.';
    ln{1,7} = '-.';
    
    ln{2,1} = '-';
    ln{2,2} = '--';
    ln{2,3} = '-.';
    ln{2,4} = '-.';
    ln{2,5} = '-.';
    ln{2,6} = '-.';
    ln{2,7} = '-.';
    
    ln{3,1} = '-';
    ln{3,2} = '--';
    ln{3,3} = '-.';
    ln{3,4} = '-.';
    ln{3,5} = '-.';
    ln{3,6} = '-.';
    ln{3,7} = '-.';

    mrk = cell(3,7);
    mrk{1,1} = 'o';
    mrk{1,2} = '^';
    mrk{1,3} = 's';
    mrk{1,4} = 's';
    mrk{1,5} = 's';
    mrk{1,6} = 's';
    mrk{1,7} = 's';
    
    mrk{2,1} = 'o';
    mrk{2,2} = '^';
    mrk{2,3} = 's';
    mrk{2,4} = 's';
    mrk{2,5} = 's';
    mrk{2,6} = 's';
    mrk{2,7} = 's';
    
    mrk{3,1} = 'o';
    mrk{3,2} = '^';
    mrk{3,3} = 's';
    mrk{3,4} = 's';
    mrk{3,5} = 's';
    mrk{3,6} = 's';
    mrk{3,7} = 's';
    
    prdOptimDataForStop = prdOptimData;
    
end

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
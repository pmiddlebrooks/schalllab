function cmd_session_data_plot(Data, options)

%%
% Set defaults



targAngleArray = Data(1).targAngleArray;
ssdArray    = Data(1).ssdArray;
sessionID   = Data(1).sessionID;


dataType    = options.dataType;
printPlot   = options.printPlot;
filterData  = options.filterData;
stopHz      = options.stopHz;
figureHandle = options.figureHandle;
targAngle  = options.targAngle;


PLOT_ERROR  = true;

epochArray = {'targOn', 'stopSignalOn', 'responseOnset', 'rewardOn'};

nUnit = length(Data);

if strcmp(options.targAngle, 'collapse')
   nTarg = 1;
elseif strcmp(options.targAngle, 'each')
   nTarg = length(targAngleArray);
else
   nTarg = length(options.targAngle);
end


cMap = cmd_colormap;


targLineW = 2;

for kDataIndex = 1 : nUnit
   nRow = nTarg;
   nEpoch = length(epochArray);
   nColumn = nEpoch;
   figureHandle = figureHandle + 1;
   if printPlot
      [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
   else
      [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
   end
   clf
   
   yLimMax = Data(kDataIndex).yMax * 1.1;
   yLimMin = min(Data(kDataIndex).yMin * 1.1);
   
   %    yLimMin = 0;
   %    yLimMax = 65;
   
   
   
   
   for iTarg = 1 : nTarg
      
      for mEpoch = 1 : nEpoch
         mEpochName = epochArray{mEpoch};
         epochRange = ccm_epoch_range(mEpochName, 'plot');
         
         % _______  Set up axes  ___________
         
         
         % Set up plot axes
         % Left target Go trials
         ax(iTarg, mEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(iTarg, mEpoch) yAxesPosition(iTarg, mEpoch) axisWidth axisHeight]);
         set(ax(iTarg, mEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
         cla
         hold(ax(iTarg, mEpoch), 'on')
         plot(ax(iTarg, mEpoch), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
         title(epochArray{mEpoch})
         
         if mEpoch == 1
             yl = sprintf('Angle: %d', targAngleArray(iTarg));
             ylabel(ax(iTarg,mEpoch), yl)
         end
         if mEpoch > 1
            set(ax(iTarg, mEpoch), 'yticklabel', [], 'ycolor', [1 1 1])            
%             set(ax(iTarg, mEpoch), 'ycolor', [1 1 1])
         end
         
         
         
         
         % Go trials
         switch dataType
            case 'neuron'
               dataSignal = 'sdfMean';
            case 'lfp'
               dataSignal = 'lfpMean';
            case 'erp'
               dataSignal = 'erp';
         end
         
         if ~strcmp(mEpochName, 'stopSignalOn')  % No stop signals on go trials
            
            alignGoTarg = Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).alignTime;
            
            if ~isempty(alignGoTarg)
               sigGoTarg = Data(kDataIndex).angle(iTarg).goTarg.(mEpochName).(dataSignal);
               %                         plot(ax(iTarg, mEpoch), epochRange, sigGoTarg(alignGoTarg + epochRange), 'color', goC(iPropIndex,:), 'linewidth', targLineW)
               plot(ax(iTarg, mEpoch), epochRange, sigGoTarg(alignGoTarg + epochRange), 'color', cMap.goTarg, 'linewidth', targLineW)
            end
         end
         
         
         % Stop signal trials
         switch dataType
            case 'neuron'
               dataSignal = 'sdf';
            case 'lfp'
               dataSignal = 'lfp';
            case 'erp'
               dataSignal = 'eeg';
         end
         stopTargSig = cell(1, length(ssdArray));
         stopTargAlign = cell(1, length(ssdArray));
         stopStopSig = cell(1, length(ssdArray));
         stopStopAlign = cell(1, length(ssdArray));
         for jSSDIndex = 1 : length(ssdArray)
            stopTargSig{jSSDIndex} = Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).(dataSignal);
            stopTargAlign{jSSDIndex} = Data(kDataIndex).angle(iTarg).stopTarg.ssd(jSSDIndex).(mEpochName).alignTime;
            
            
            if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
               stopStopSig{jSSDIndex} = Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).(dataSignal);
               stopStopAlign{jSSDIndex} = Data(kDataIndex).angle(iTarg).stopStop.ssd(jSSDIndex).(mEpochName).alignTime;
            end
            
         end  % jSSDIndex = 1 : length(ssdArray)
         
         [rasStopTarg, alignStopTarg] = align_raster_sets(stopTargSig, stopTargAlign);
         switch dataType
            case 'neuron'
               sigStopTarg = nanmean(rasStopTarg, 1);
            case {'lfp','erp'}
               if filterData
                  sigStopTarg = lowpass(nanmean(rasStopTarg, 1)', stopHz)';
               else
                  sigStopTarg = nanmean(rasStopTarg, 1);
               end
         end
         if size(sigStopTarg, 2) == 1, sigStopTarg = []; end;
         
         if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
            [rasStopCorrect, alignStopCorrect] = align_raster_sets(stopStopSig, stopStopAlign);
            switch dataType
               case 'neuron'
                  sigStopCorrect = nanmean(rasStopCorrect, 1);
               case {'lfp','erp'}
                  if filterData
                     sigStopCorrect = lowpass(nanmean(rasStopCorrect, 1)', stopHz)';
                  else
                     sigStopCorrect = nanmean(rasStopCorrect, 1);
                  end
            end
            if size(sigStopCorrect, 2) == 1, sigStopCorrect = []; end;
         end
         
         
         
         
         if ~isempty(sigStopTarg)
            plot(ax(iTarg, mEpoch), epochRange, sigStopTarg(alignStopTarg + epochRange), 'color', cMap.stopTarg, 'linewidth', targLineW)
         end
         
         if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
            if ~isempty(sigStopCorrect)
               plot(ax(iTarg, mEpoch), epochRange, sigStopCorrect(alignStopCorrect + epochRange), 'color', cMap.stopStop, 'linewidth', targLineW)
            end
         end
         
         
         
%          if mEpoch == 1
%       legend('Go','StopTarg','StopStop')
%          end

         
         
      end % mEpoch
      
   end
   %                             legend(ax(iTarg, 1), {num2cell(targAngleArray'), num2str(targAngleArray')})
   
   %         colorbar('peer', ax(iTarg, 1), 'location', 'west')
   %         colorbar('peer', ax(iTarg, 1), 'location', 'west')
   h=axes('Position', [0 0 1 1], 'Visible', 'Off');
   titleString = sprintf('%s \t %s', sessionID, Data(kDataIndex).name);
   text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top','interpreter','none')
   if printPlot
      print(figureHandle,[local_figure_path, sessionID, '_', Data(kDataIndex).name, '_ccm_',dataType,'.pdf'],'-dpdf', '-r300')
   end
end % kDataIndex



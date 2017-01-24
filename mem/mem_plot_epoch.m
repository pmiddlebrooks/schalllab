function DataOut = mem_plot_epoch(subjectID, sessionID, options)

% If called without options structure, returns default options structure to tailor.
% If called with only DataIn input, assign default option structure.
%
% options:
%
% epochName: 


if nargin < 3
%     options.epochName     	= 'targOn';
%     options.eventMark1   	= 'responseCue';
%     options.eventMark2   	= 'responseOnset';
%     options.epochRange   	= [-1499 : 1500];
    options.yAxis           = [];
    options.targSide           = 'rightTarg'; % or leftTarg
    options.printPlot           = 'false';
end
if nargin == 0
    DataOut = options;
    return
end


opt = mem_session_data;
opt.printPlot = false;
opt.plotFlag = false;
DataIn = mem_session_data(subjectID,sessionID,opt);

   
   
DataOut            = struct();

% epochName       = options.epochName;
% eventMark1      = options.eventMark1;
% eventMark2      = options.eventMark2;
figureHandle    = 20;
% dataType    	= DataIn(1).options.dataType;
% epochRange      = options.epochRange;
epochTarg       = -400 : 2200;
epochSacc       = -2200 : 400;
printPlot       = options.printPlot;
yAxisRange      = options.yAxis;
targSide        = options.targSide; % Option to plot data in only one direction

nUnit           = length(DataIn);
dataType        = DataIn.dataType;


plotFlag        = 1;
rightColor      = [0 0 1];
leftColor       = rightColor;
markerColorOpen    	= [.5 .5 .5];
markerColor    	= [0 0 0];
markerSize    	= 60;


if plotFlag
    targLineW = 2;
    tickWidth = 6;
    
    for k = 1 : nUnit
        
        
        
        %  % Figure out a good y-axis limit to use
        if ~isempty(yAxisRange)
            yMin = yAxisRange(1);
            yMax = yAxisRange(2);
        else
            yMin = min([0, min(DataIn(k).leftTarg.targOn.signalMean), min(DataIn(k).leftTarg.responseOnset.signalMean), min(DataIn(k).rightTarg.targOn.signalMean), min(DataIn(k).rightTarg.responseOnset.signalMean)]);
            yMax = max([max(DataIn(k).leftTarg.targOn.signalMean), max(DataIn(k).leftTarg.responseOnset.signalMean), max(DataIn(k).rightTarg.targOn.signalMean), max(DataIn(k).rightTarg.responseOnset.signalMean)]);
            %          yMax = DataIn(k).yMax;
        end
        yLimRas  = 40; % What's a reasonable y-axis for the raster plots?
        
        
        nLeftTrial = size(DataIn(k).leftTarg.targOn.raster, 1);
        nRightTrial = size(DataIn(k).rightTarg.targOn.raster, 1);
        
        
        
        nRow = 2;
        nColumn = 2;
        
        screenOrSave = 'screen';
        figureHandle = figureHandle + 1;
        if printPlot
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
        else
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
        end
        clf
        
        %                                 neuronTitle = sprintf('%s', DataIn(k).name);
        %         titleH = title(neuronTitle, 'FontWeight', 'Bold', 'color', 'r');
        %         taxH = ax(1, axLeftTarg);
        %         set(titleH,'Position',[taxH(1)*1.2 taxH(4)*.8 0]);
        
        %         boxMargin = .5;
        %         x = xAxesPosition(end, 1);% - boxMargin;
        %         y = yAxesPosition(end, 1);% - boxMargin;
        %         w = axisWidth * nColumn/2;
        %         h = axisHeight * nRow;
        %                 rectangle('Position', [x, y, w, h], 'edgecolor', 'b')
        
        
        
        % ______________  SET UP AXES  __________________
        % axes names
        axSig = 1; % SDFs on top
        axRas = 2; % Rasters on bottom plot
        axLeftTarg = 1;
        axLeftSacc = 2;
        axRightTarg = 3;
        axRightSacc = 4;
        
        % LEFT TARGETS
        % Left SDFs
        ax(axSig, axLeftTarg) = axes('units', 'centimeters', 'position', [xAxesPosition(axSig, axLeftTarg) yAxesPosition(axSig, axLeftTarg) axisWidth axisHeight]);
        set(ax(axSig, axLeftTarg), 'ylim', [yMin yMax * 1.1], 'xlim', [epochTarg(1) epochTarg(end)])
        cla
        hold(ax(axSig, axLeftTarg), 'on')
        plot(ax(axSig, axLeftTarg), [1 1], [yMin yMax * 1.1], '-k', 'linewidth', 2)
        %         ttl = sprintf('%s', 'target On');
        %         title(ttl)
        
        % Left Rasters
        ax(axRas, axLeftTarg) = axes('units', 'centimeters', 'position', [xAxesPosition(axRas, axLeftTarg) yAxesPosition(axRas, axLeftTarg) axisWidth axisHeight]);
                set(ax(axRas, axLeftTarg), 'ylim', [0 max([yLimRas, nLeftTrial])], 'xlim', [epochTarg(1) epochTarg(end)])
%         set(ax(axRas, axLeftTarg), 'ylim', [yMin yMax * 1.1], 'xlim', [epochRange(1) epochRange(end)])
        cla
        hold(ax(axRas, axLeftTarg), 'on')
        plot(ax(axRas, axLeftTarg), [1 1], [0 nLeftTrial], '-k', 'linewidth', 2)
        
        % LEFT SACCADES
        % Left SDFs
        ax(axSig, axLeftSacc) = axes('units', 'centimeters', 'position', [xAxesPosition(axSig, axLeftSacc) yAxesPosition(axSig, axLeftSacc) axisWidth axisHeight]);
        set(ax(axSig, axLeftSacc), 'ylim', [yMin yMax * 1.1], 'xlim', [epochSacc(1) epochSacc(end)])
        cla
        hold(ax(axSig, axLeftSacc), 'on')
        plot(ax(axSig, axLeftSacc), [1 1], [yMin yMax * 1.1], '-k', 'linewidth', 2)
        %         ttl = sprintf('%s', 'target On');
        %         title(ttl)
        
        % Left Rasters
        ax(axRas, axLeftSacc) = axes('units', 'centimeters', 'position', [xAxesPosition(axRas, axLeftSacc) yAxesPosition(axRas, axLeftSacc) axisWidth axisHeight]);
                set(ax(axRas, axLeftSacc), 'ylim', [0 max([yLimRas, nLeftTrial])], 'xlim', [epochSacc(1) epochSacc(end)])
%         set(ax(axRas, axLeftSacc), 'ylim', [yMin yMax * 1.1], 'xlim', [epochRange(1) epochRange(end)])
        cla
        hold(ax(axRas, axLeftSacc), 'on')
        plot(ax(axRas, axLeftSacc), [1 1], [0 nLeftTrial], '-k', 'linewidth', 2)
        
        
        
%         % RIGHT TARGETS
%         % Right SDFs
%         ax(axSig, axRightTarg) = axes('units', 'centimeters', 'position', [xAxesPosition(axSig, axRightTarg) yAxesPosition(axSig, axRightTarg) axisWidth axisHeight]);
%         set(ax(axSig, axRightTarg), 'ylim', [yMin yMax * 1.1], 'xlim', [epochTarg(1) epochTarg(end)])
%         cla
%         hold(ax(axSig, axRightTarg), 'on')
%         plot(ax(axSig, axRightTarg), [1 1], [yMin yMax * 1.1], '-k', 'linewidth', 2)
%         %         ttl = sprintf('%s', epochName);
%         %         title(ttl)
%         
%         % Right Rasters
%         ax(axRas, axRightTarg) = axes('units', 'centimeters', 'position', [xAxesPosition(axRas, axRightTarg) yAxesPosition(axRas, axRightTarg) axisWidth axisHeight]);
%                 set(ax(axRas, axRightTarg), 'ylim', [0 max([yLimRas, nRightTrial])], 'xlim', [epochTarg(1) epochTarg(end)])
% %         set(ax(axRas, axRightTarg), 'ylim', [yMin yMax * 1.1], 'xlim', [epochRange(1) epochRange(end)])
%         cla
%         hold(ax(axRas, axRightTarg), 'on')
%         plot(ax(axRas, axRightTarg), [1 1], [0 nRightTrial], '-k', 'linewidth', 2)
%         
%         % RIGHT SACCADES
%         % Right SDFs
%         ax(axSig, axRightSacc) = axes('units', 'centimeters', 'position', [xAxesPosition(axSig, axRightSacc) yAxesPosition(axSig, axRightSacc) axisWidth axisHeight]);
%         set(ax(axSig, axRightSacc), 'ylim', [yMin yMax * 1.1], 'xlim', [epochSacc(1) epochSacc(end)])
%         cla
%         hold(ax(axSig, axRightSacc), 'on')
%         plot(ax(axSig, axRightSacc), [1 1], [yMin yMax * 1.1], '-k', 'linewidth', 2)
%         %         ttl = sprintf('%s', epochName);
%         %         title(ttl)
%         
%         % Right Rasters
%         ax(axRas, axRightSacc) = axes('units', 'centimeters', 'position', [xAxesPosition(axRas, axRightSacc) yAxesPosition(axRas, axRightSacc) axisWidth axisHeight]);
%                 set(ax(axRas, axRightSacc), 'ylim', [0 max([yLimRas, nRightTrial])], 'xlim', [epochSacc(1) epochSacc(end)])
% %         set(ax(axRas, axRighOtSacc), 'ylim', [yMin yMax * 1.1], 'xlim', [epochRange(1) epochRange(end)])
%         cla
%         hold(ax(axRas, axRightSacc), 'on')
%         plot(ax(axRas, axRightSacc), [1 1], [0 nRightTrial], '-k', 'linewidth', 2)
        
           colormap([1 1 1; 0 0 0])
     
        
        
        
        
        
        
        
        
        % ______________  PLOT DATA  __________________
        
        % Left Trials Target Aligned
        if ~isempty(DataIn(k).leftTarg.targOn.alignTime)
            % Signals
            plot(ax(axSig, axLeftTarg), epochTarg, DataIn(k).(targSide).targOn.signalMean(DataIn(k).(targSide).targOn.alignTime + epochTarg), 'color', leftColor, 'linewidth', targLineW)
            plot(ax(axSig, axLeftSacc), epochTarg, DataIn(k).(targSide).responseOnset.signalMean(DataIn(k).(targSide).responseOnset.alignTime + epochTarg), 'color', leftColor, 'linewidth', targLineW)
            
            % Rasters
            if strcmp(dataType, 'neuron')
                axes(ax(axRas, axLeftTarg))
                
                
                % Events on Rasters - for now these are hard coded
                leftCueLat = get_event_latency(DataIn(k), targSide, 'targOn', 'responseCueOn');
                leftRTLat = get_event_latency(DataIn(k), targSide, 'targOn', 'responseOnset');

                % Sort trials based on RT
                [B,iRTOrder] = sort(leftRTLat);
                
                leftRasTarg = fat_raster(DataIn(k).(targSide).targOn.raster(iRTOrder,:), tickWidth);
                imagesc(epochTarg, 1 : nLeftTrial, leftRasTarg(:, DataIn(k).(targSide).targOn.alignTime + epochTarg))

                scatter(leftCueLat(iRTOrder,:), 1 : nLeftTrial, 'd', 'markeredgecolor', markerColor, 'markerfacecolor', markerColorOpen, 'sizeData', markerSize)
                scatter(leftRTLat(iRTOrder,:), 1 : nLeftTrial, 'o', 'markeredgecolor', markerColor, 'markerfacecolor', markerColor, 'sizeData', markerSize)

                
             
                % Alignment line
                plot(ax(axRas, axLeftTarg), [1 1], [0 nLeftTrial], '-k', 'linewidth', 2)
            end
        end
        
        % Left Trials Saccade Aligned
        if ~isempty(DataIn(k).(targSide).responseOnset.alignTime)
            % Signals
            plot(ax(axSig, axLeftSacc), epochSacc, DataIn(k).(targSide).responseOnset.signalMean(DataIn(k).(targSide).responseOnset.alignTime + epochSacc), 'color', leftColor, 'linewidth', targLineW)
            
            % Rasters
            if strcmp(dataType, 'neuron')
                axes(ax(axRas, axLeftSacc))

                % Events on Rasters - for now these are hard coded
                leftCueLat = get_event_latency(DataIn(k), targSide, 'responseOnset', 'responseCueOn');
                leftTargLat = get_event_latency(DataIn(k), targSide, 'responseOnset', 'targOn');
               	
                [B,iTargOrder] = sort(leftTargLat);
                iTargOrder = wrev(iTargOrder);

               	leftRasSacc = fat_raster(DataIn(k).(targSide).responseOnset.raster(iTargOrder,:), tickWidth);               
                imagesc(epochSacc, 1 : nLeftTrial, leftRasSacc(:, DataIn(k).(targSide).responseOnset.alignTime + epochSacc))
                
                scatter(leftCueLat(iTargOrder,:), 1 : nLeftTrial, 'd', 'markeredgecolor', markerColor, 'markerfacecolor', markerColorOpen, 'sizeData', markerSize)
                scatter(leftTargLat(iTargOrder,:), 1 : nLeftTrial, 'o', 'markeredgecolor', markerColor, 'markerfacecolor', markerColor, 'sizeData', markerSize)

                % Alignment line
                plot(ax(axRas, axLeftSacc), [1 1], [0 nLeftTrial], '-k', 'linewidth', 2)
            end
        end
        
        
        
        
%         % Right Target Trials Target Aligned
%         if ~isempty(DataIn(k).rightTarg.targOn.alignTime)
%             % Signal
%             plot(ax(axSig, axRightTarg), epochTarg, DataIn(k).rightTarg.targOn.signalMean(DataIn(k).rightTarg.targOn.alignTime + epochTarg), 'color', rightColor, 'linewidth', targLineW)
%             plot(ax(axSig, axRightSacc), epochTarg, DataIn(k).rightTarg.responseOnset.signalMean(DataIn(k).rightTarg.responseOnset.alignTime + epochTarg), 'color', rightColor, 'linewidth', targLineW)
%             
%             % Raster
%             if strcmp(dataType, 'neuron')
%                 axes(ax(axRas, axRightTarg))
% 
%                 % Events on Rasters
%                 rightCueLat = get_event_latency(DataIn(k), 'rightTarg', 'targOn', 'responseCueOn');
%                 rightRTLat = get_event_latency(DataIn(k), 'rightTarg', 'targOn', 'responseOnset');
%                 
%                  % Sort trials based on RT
%                 [B,iRTOrder] = sort(rightRTLat);
% 
%                 rightRasTarg = fat_raster(DataIn(k).rightTarg.targOn.raster(iRTOrder,:), tickWidth);                
%                 imagesc(epochTarg, 1 : nRightTrial, rightRasTarg(:, DataIn(k).rightTarg.targOn.alignTime + epochTarg))
% 
%                 scatter(rightCueLat(iRTOrder), 1 : nRightTrial, 'o', 'markeredgecolor', markerColor, 'sizeData', markerSize)
%                 scatter(rightRTLat(iRTOrder), 1 : nRightTrial, 'o', 'markeredgecolor', markerColor, 'markerfacecolor', markerColor, 'sizeData', markerSize)
%                 
%                 % Alignment line
%                 plot(ax(axRas, axRightTarg), [1 1], [0 nRightTrial], '-k', 'linewidth', 2)
%             end
%         end
%         
%         
%         % Right Target Trials Saccade Aligned
%         if ~isempty(DataIn(k).rightTarg.responseOnset.alignTime)
%             % Signal
%             plot(ax(axSig, axRightSacc), epochSacc, DataIn(k).rightTarg.responseOnset.signalMean(DataIn(k).rightTarg.responseOnset.alignTime + epochSacc), 'color', rightColor, 'linewidth', targLineW)
%             
%             % Raster
%             if strcmp(dataType, 'neuron')
%                 axes(ax(axRas, axRightSacc))
%                 
%                 % Events on Rasters
%                 rightCueLat = get_event_latency(DataIn(k), 'rightTarg', 'responseOnset', 'responseCueOn');
%                 rightTargLat = get_event_latency(DataIn(k), 'rightTarg', 'responseOnset', 'targOn');
% 
%                  % Sort trials based on RT
%                 [B,iTargOrder] = sort(rightTargLat);
%                 iTargOrder = wrev(iTargOrder);
% 
%                 rightRasSacc = fat_raster(DataIn(k).rightTarg.responseOnset.raster(iTargOrder,:), tickWidth);                
%                 imagesc(epochSacc, 1 : nRightTrial, rightRasSacc(:, DataIn(k).rightTarg.responseOnset.alignTime + epochSacc))
% 
%                 scatter(rightCueLat(iTargOrder), 1 : nRightTrial, 'o', 'markeredgecolor', markerColor, 'sizeData', markerSize)
%                 scatter(rightTargLat(iTargOrder), 1 : nRightTrial, 'o', 'markeredgecolor', markerColor, 'markerfacecolor', markerColor, 'sizeData', markerSize)
%                 
%                 % Alignment line
%                 plot(ax(axRas, axRightSacc), [1 1], [0 nRightTrial], '-k', 'linewidth', 2)
%             end
%         end
        
      h=axes('Position', [0 0 1 1], 'Visible', 'Off');
      titleString = sprintf('%s \t %s', DataIn(k).sessionID, DataIn(k).name);
      text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top', 'color', 'k')
      if printPlot
         localFigurePath = local_figure_path;
         print(figureHandle,[localFigurePath, DataIn(k).sessionID, '_', DataIn(k).name, '_mem_plot_' dataType],'-dpdf', '-r300')
%          print(figureHandle,[localFigurePath, DataIn.sessionID, '_', DataIn(k).name, '_mem_plot_' dataType],'-djpeg')
      end
        
    end
end


end % function





% ******************************************
%               SUBFUNCTIONS
% ******************************************

function [eventLatency] = get_event_latency(UnitCondition, targSide, epochName, eventMarkName)

eventLatency = [];
alignTime = UnitCondition.(targSide).(epochName).alignTime;
unitAlignList = UnitCondition.(targSide).(epochName).alignTimeList;


if ~isempty(eventMarkName)  % If an eventMarkName was entered, use it
    eventAlignList = UnitCondition.(targSide).(eventMarkName).alignTimeList;
else  % Else use defaults
    switch epochName
        case 'fixWindowEntered'
            % target onset
            eventAlignList = UnitCondition.(targSide).targOn.alignTimeList;
        case 'targOn'
            % target onset
            eventAlignList = UnitCondition.(targSide).responseOnset.alignTimeList;
        case 'responseOnset'
            % target onset
            eventAlignList = UnitCondition.(targSide).targOn.alignTimeList;
        case 'rewardOn'
            % target onset
            eventAlignList = UnitCondition.(targSide).responseOnset.alignTimeList;
    end
end
if ~isempty(alignTime)
    eventLatency = eventAlignList - unitAlignList;
end
end



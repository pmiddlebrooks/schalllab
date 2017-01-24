function data = cmd_plot_epoch(Unit, options)

% If singal strength or ssd index vectors are not input, assume user wants to
% collapse across all of them
if nargin < 2
    options.epochName = 'responseOnset';
    options.eventMarkName = 'targOn';
    options.ssdArray = Unit(1).ssdArray;
   options.targAngle        = 'each';
    if nargin == 0
        data = options;
    end
end
epochName = options.epochName;
eventMarkName = options.eventMarkName;
ssdArray = options.ssdArray;

    switch options.targAngle
        case 'each'
            angleArray = Unit(1).targAngleArray;
            nTarg = length(angleArray);
        case 'collapse'
            nTarg = 1;
    end
        
nUnit = length(Unit);


epochRange = ccm_epoch_range(epochName, 'plot');
% sdfMax = 80;




plotFlag = 1;
goColor         = [0 0 1];
stopGoColor     = [1 0 0];
stopStopColor   = [.5 .5 .5];
markerColor       = [0 0 0];
markerSize        = 40;
if plotFlag
    
    
    figureHandle = 674;
    targLineW = 2;
    tickWidth = 10;
    
    
    
    for k = 1 : nUnit
        
        
        
        
        nRow = 2;
        nColumn = 2;
        screenOrSave = 'screen';
        figureHandle = figureHandle + 1;
        [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = big_figure(nRow, nColumn, figureHandle, screenOrSave);
        clf
        
   for iTarg = 1 : nTarg;
      
    GoTarg      = cmd_concat_neural_conditions(Unit(k), epochName, eventMarkName, {'goTarg'}, angleArray(iTarg), ssdArray);
    StopTarg    = cmd_concat_neural_conditions(Unit(k), epochName, eventMarkName, {'stopTarg'}, angleArray(iTarg), ssdArray);
    StopCorrect = cmd_concat_neural_conditions(Unit(k), epochName, eventMarkName, {'stopStop'}, angleArray(iTarg), ssdArray);

    nGoTargTrial    = size(GoTarg.raster, 1);
        nStopTargTrial  = size(StopTarg.raster, 1);
        nStopCorrectTrial = size(StopCorrect.raster, 1);
        
        nTargTrial = nStopCorrectTrial + nStopTargTrial + nGoTargTrial;
        
        
        %  % Figure out a good y-axis limit to use
        sdfMax = max([1; GoTarg.sdf'; StopTarg.sdf'; StopCorrect.sdf']);
        
      
      if strcmp(options.targAngle, 'collapse')
         selectOpt.targAngle = 'all';
      elseif strcmp(options.targAngle, 'each')
         selectOpt.targAngle = angleArray(iTarg);
      else
         selectOpt.targAngle = options.targAngle;
      end
        
        
        
        % _______  Set up axes  ___________
        % axes names
        axTarg = 1;
        
        % Set up plot axes
        % to target (and stop correct) trials
        ax(1, axTarg) = axes('units', 'centimeters', 'position', [xAxesPosition(1, axTarg) yAxesPosition(1, axTarg) axisWidth axisHeight]);
        set(ax(1, axTarg), 'ylim', [0 sdfMax * 1.1], 'xlim', [epochRange(1) epochRange(end)])
        cla
        hold(ax(1, axTarg), 'on')
        plot(ax(1, axTarg), [1 1], [0 sdfMax * 1.1], '-k', 'linewidth', 2)
        ttl = sprintf('Target:  %s', epochName);
        title(ttl)
        
        
        ax(2, axTarg) = axes('units', 'centimeters', 'position', [xAxesPosition(2, axTarg) yAxesPosition(2, axTarg) axisWidth axisHeight]);
        set(ax(2, axTarg), 'ylim', [0 nTargTrial], 'xlim', [epochRange(1) epochRange(end)])
        cla
        hold(ax(2, axTarg), 'on')
        plot(ax(2, axTarg), [1 1], [0 nTargTrial], '-k', 'linewidth', 2)
        
        
        
        
        % Plot SDFs
        if ~isempty(StopCorrect(k).align)
            plot(ax(1, axTarg), epochRange, StopCorrect(k).sdf(StopCorrect(k).align + epochRange), 'color', stopStopColor, 'linewidth', targLineW)
        end
        if ~isempty(StopTarg(k).align)
            plot(ax(1, axTarg), epochRange, StopTarg(k).sdf(StopTarg(k).align + epochRange), 'color', stopGoColor, 'linewidth', targLineW)
        end
        if ~isempty(GoTarg(k).align)
            plot(ax(1, axTarg), epochRange, GoTarg(k).sdf(GoTarg(k).align + epochRange), 'color', goColor, 'linewidth', targLineW)
        end
        
        
        
        % Plot Rasters
        stopCorrRas = fat_raster(StopCorrect(k).raster, tickWidth);
        stopCorrRas = stopCorrRas .* 3;
        stopTargRas = fat_raster(StopTarg(k).raster, tickWidth);
        stopTargRas = stopTargRas .* 2;
        goTargRas = fat_raster(GoTarg(k).raster, tickWidth);

        
        colormap([1 1 1; goColor; stopGoColor; stopStopColor])
        axes(ax(2, axTarg))
        if ~isempty(StopCorrect(k).align)
            imagesc(epochRange, 1 : nStopCorrectTrial, stopCorrRas(:, StopCorrect(k).align + epochRange))
            if ~isempty(StopCorrect(k).eventLatency)
                scatter(StopCorrect(k).eventLatency, 1 : nStopCorrectTrial, '.', 'markeredgecolor', markerColor, 'markerfacecolor', markerColor, 'sizeData', markerSize)
            end
        end
        
        if ~isempty(StopTarg(k).align)
            imagesc(epochRange, nStopCorrectTrial + 1 : nStopCorrectTrial + nStopTargTrial, stopTargRas(:, StopTarg(k).align + epochRange))
            if ~isempty(StopTarg(k).eventLatency)
                scatter(StopTarg(k).eventLatency, nStopCorrectTrial + 1 : nStopCorrectTrial + nStopTargTrial, '.', 'markeredgecolor', markerColor, 'markerfacecolor', markerColor, 'sizeData', markerSize)
            end
        end
        
        if ~isempty(GoTarg(k).align)
            imagesc(epochRange, nStopCorrectTrial + nStopTargTrial + 1 : nStopCorrectTrial + nStopTargTrial + nGoTargTrial, goTargRas(:, GoTarg(k).align + epochRange))
            if ~isempty(GoTarg(k).eventLatency)
                scatter(GoTarg(k).eventLatency, nStopCorrectTrial + nStopTargTrial + 1 : nStopCorrectTrial + nStopTargTrial + nGoTargTrial, '.', 'markeredgecolor', markerColor, 'markerfacecolor', markerColor, 'sizeData', markerSize)
            end
        end
        
        plot(ax(2, axTarg), [1 1], [0 nTargTrial], '-k', 'linewidth', 2)
        
        
        
   end % kUnit
    end % iTarg
end



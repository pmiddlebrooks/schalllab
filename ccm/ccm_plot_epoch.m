function Data = ccm_plot_epoch(Unit, options)

% If called without options structure, returns default options structure to tailor.
% If called with only Unit input, assign default option structure.
%

if nargin == 0
    options.epochName       	= 'checkerOn';
    options.eventMarkName   	= 'responseOnset';
    options.epochRange   	= [];
    Data = options;
    return
end
if nargin == 1
    options.epochName       	= 'checkerOn';
    options.eventMarkName   	= 'responseOnset';
    options.epochRange   	= [];    
end

epochName       = options.epochName;
eventMarkName   = options.eventMarkName;
printPlot       = false;
figureHandle    = 674;
ssdArray        = Unit(1).ssdArray;
colorCohArray   = Unit(1).pSignalArray;
dataType        = Unit(1).options.dataType;
epochRange      = options.epochRange;


nUnit = length(Unit);

if isempty(epochRange)
rangeFactor = 2;
epochRange = ccm_epoch_range(epochName, 'plot', rangeFactor);
end
% yMax = 80;




plotFlag = 1;
goColor         = [0 0 1];
stopGoColor     = [1 0 0];
stopStopColor   = [.5 .5 .5];
markerColor       = [0 0 0];
markerSize        = 40;
if plotFlag
    
    
    targLineW = 2;
    distLineW = 1;
    tickWidth = 6;
    
    opt = ccm_concat_neural_conditions(Unit); % Get default options structure
    
    opt.epochName = epochName;
    opt.eventMarkName = eventMarkName;
    opt.conditionArray = {'goTarg'};
    opt.colorCohArray = colorCohArray;
    opt.ssdArray = ssdArray;
    opt.dataType = dataType;
    GoTarg      = ccm_concat_neural_conditions(Unit, opt);
    opt.conditionArray = {'goDist'};
    GoDist      = ccm_concat_neural_conditions(Unit, opt);
    opt.conditionArray = {'stopTarg'};
    StopTarg    = ccm_concat_neural_conditions(Unit, opt);
    opt.conditionArray = {'stopDist'};
    StopDist    = ccm_concat_neural_conditions(Unit, opt);
    opt.conditionArray = {'stopStop'};
    StopCorrect = ccm_concat_neural_conditions(Unit, opt);
    
    
    for k = 1 : nUnit
        
        nGoTargTrial    = size(GoTarg(k).signal, 1);
        nGoDistTrial    = size(GoDist(k).signal, 1);
        nStopTargTrial  = size(StopTarg(k).signal, 1);
        nStopDistTrial  = size(StopDist(k).signal, 1);
        nStopCorrectTrial = size(StopCorrect(k).signal, 1);
        
        nTargTrial = nStopCorrectTrial + nStopTargTrial + nGoTargTrial;
        
        
        %  % Figure out a good y-axis limit to use
        if strcmp(dataType, 'neuron')
        yMax = max([1; GoTarg(k).signalFn'; StopTarg(k).signalFn'; StopCorrect(k).signalFn']);
        else
                    yMax = max([GoTarg(k).signalFn'; StopTarg(k).signalFn'; StopCorrect(k).signalFn']);
        end
        yMin = min([0; GoTarg(k).signalFn'; StopTarg(k).signalFn'; StopCorrect(k).signalFn']);
        
        
        
        
        nRow = 2;
        nColumn = 2;
        screenOrSave = 'screen';
        figureHandle = figureHandle + 1;
        [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = big_figure(nRow, nColumn, figureHandle, screenOrSave);
        clf
        
        %                                 neuronTitle = sprintf('%s', Unit(k).name);
        %         titleH = title(neuronTitle, 'FontWeight', 'Bold', 'color', 'r');
        %         taxH = ax(1, axTarg);
        %         set(titleH,'Position',[taxH(1)*1.2 taxH(4)*.8 0]);
        
        %         boxMargin = .5;
        %         x = xAxesPosition(end, 1);% - boxMargin;
        %         y = yAxesPosition(end, 1);% - boxMargin;
        %         w = axisWidth * nColumn/2;
        %         h = axisHeight * nRow;
        %                 rectangle('Position', [x, y, w, h], 'edgecolor', 'b')
        
        
        
        % _______  Set up axes  ___________
        % axes names
        axTarg = 1;
        axDist = 2;
        
        % Set up plot axes
        % to target (and stop correct) trials
        ax(1, axTarg) = axes('units', 'centimeters', 'position', [xAxesPosition(1, axTarg) yAxesPosition(1, axTarg) axisWidth axisHeight]);
        set(ax(1, axTarg), 'ylim', [yMin yMax * 1.1], 'xlim', [epochRange(1) epochRange(end)])
        cla
        hold(ax(1, axTarg), 'on')
        plot(ax(1, axTarg), [1 1], [yMin yMax * 1.1], '-k', 'linewidth', 2)
        ttl = sprintf('Target:  %s', epochName);
        title(ttl)
        
        % to distractor (and stop correct) trials
        ax(1, axDist) = axes('units', 'centimeters', 'position', [xAxesPosition(1, axDist) yAxesPosition(1, axDist) axisWidth axisHeight]);
        set(ax(1, axDist), 'ylim', [yMin yMax * 1.1], 'xlim', [epochRange(1) epochRange(end)])
        cla
        hold(ax(1, axDist), 'on')
        plot(ax(1, axDist), [1 1], [yMin yMax * 1.1], '-k', 'linewidth', 2)
        ttl = sprintf('Distractor:  %s', epochName);
        title(ttl)
        %             set(ax(axRight, mEpoch), 'Xtick', [0 : 100 : epochRange(end) - epochRange(1)], 'XtickLabel', [epochRange(1) : 100: epochRange(end)])
        
        
        ax(2, axTarg) = axes('units', 'centimeters', 'position', [xAxesPosition(2, axTarg) yAxesPosition(2, axTarg) axisWidth axisHeight]);
        set(ax(2, axTarg), 'ylim', [0 nTargTrial], 'xlim', [epochRange(1) epochRange(end)])
        cla
        hold(ax(2, axTarg), 'on')
        plot(ax(2, axTarg), [1 1], [0 nTargTrial], '-k', 'linewidth', 2)
        
        % to distractor (and stop correct) trials
        ax(2, axDist) = axes('units', 'centimeters', 'position', [xAxesPosition(2, axDist) yAxesPosition(2, axDist) axisWidth axisHeight]);
        set(ax(2, axDist), 'ylim', [0 nTargTrial], 'xlim', [epochRange(1) epochRange(end)])
        cla
        hold(ax(2, axDist), 'on')
        plot(ax(2, axDist), [1 1], [0 nTargTrial], '-k', 'linewidth', 2)
        
        
        
        
        % Plot signals
        if ~isempty(StopCorrect(k).align)
            plot(ax(1, axTarg), epochRange, StopCorrect(k).signalFn(StopCorrect(k).align + epochRange), 'color', stopStopColor, 'linewidth', targLineW)
        end
        if ~isempty(StopTarg(k).align)
            plot(ax(1, axTarg), epochRange, StopTarg(k).signalFn(StopTarg(k).align + epochRange), 'color', stopGoColor, 'linewidth', targLineW)
        end
        if ~isempty(GoTarg(k).align)
            plot(ax(1, axTarg), epochRange, GoTarg(k).signalFn(GoTarg(k).align + epochRange), 'color', goColor, 'linewidth', targLineW)
        end
        
        if ~isempty(StopCorrect(k).align)
            plot(ax(1, axDist), epochRange, StopCorrect(k).signalFn(StopCorrect(k).align + epochRange), 'color', stopStopColor, 'linewidth', distLineW)
        end
        if ~isempty(StopDist(k).align)
            plot(ax(1, axDist), epochRange, StopDist(k).signalFn(StopDist(k).align + epochRange), '--', 'color', stopGoColor, 'linewidth', distLineW)
        end
        if ~isempty(GoDist(k).align)
            plot(ax(1, axDist), epochRange, GoDist(k).signalFn(GoDist(k).align + epochRange), '--', 'color', goColor, 'linewidth', distLineW)
        end
        
        if strcmp(dataType, 'neuron')
        % Plot Rasters
        stopCorrRas = fat_raster(StopCorrect(k).signal, tickWidth);
        stopCorrRas = stopCorrRas .* 3;
        stopTargRas = fat_raster(StopTarg(k).signal, tickWidth);
        %         stopTargRas = stopTargRas .* (2/3);
        stopTargRas = stopTargRas .* 2;
        goTargRas = fat_raster(GoTarg(k).signal, tickWidth);
        %         goTargRas = goTargRas .* (1/3);
        
        stopDistRas = fat_raster(StopDist(k).signal, tickWidth);
        stopDistRas = stopDistRas .* 2;
        goDistRas = fat_raster(GoDist(k).signal, tickWidth);
        goDistRas = goDistRas;
        
        colormap([1 1 1; goColor; stopGoColor; stopStopColor])
        %         if isempty(StopCorrect(k).align)
        %         colormap([1 1 1; goColor; stopGoColor])
        %         elseif isempty(StopTarg(k).align)
        %         colormap([1 1 1; goColor; stopStopColor])
        %         elseif isempty(GoTarg(k).align)
        %         colormap([1 1 1; stopGoColor; stopStopColor])
        %         end
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
        
        
        
        axes(ax(2, axDist))
        if ~isempty(StopCorrect(k).align)
            imagesc(epochRange, 1 : nStopCorrectTrial, stopCorrRas(:, StopCorrect(k).align + epochRange))
            if ~isempty(StopCorrect(k).eventLatency)
                scatter(StopCorrect(k).eventLatency, 1 : nStopCorrectTrial, '.', 'markeredgecolor', markerColor, 'markerfacecolor', markerColor, 'sizeData', markerSize)
            end
        end
        if ~isempty(StopDist(k).align)
            imagesc(epochRange, nStopCorrectTrial + 1 : nStopCorrectTrial + nStopDistTrial, stopDistRas(:, StopDist(k).align + epochRange))
            if ~isempty(StopDist(k).eventLatency)
                scatter(StopDist(k).eventLatency, nStopCorrectTrial + 1 : nStopCorrectTrial + nStopDistTrial, '.', 'markeredgecolor', markerColor, 'markerfacecolor', markerColor, 'sizeData', markerSize)
            end
        end
        if ~isempty(GoDist(k).align)
            imagesc(epochRange, nStopCorrectTrial + nStopDistTrial + 1 : nStopCorrectTrial + nStopDistTrial + nGoDistTrial, goDistRas(:, GoDist(k).align + epochRange))
            if ~isempty(GoDist(k).eventLatency)
                scatter(GoDist(k).eventLatency, nStopCorrectTrial + nStopDistTrial + 1 : nStopCorrectTrial + nStopDistTrial + nGoDistTrial, '.', 'markeredgecolor', markerColor, 'markerfacecolor', markerColor, 'sizeData', markerSize)
            end
        end
        plot(ax(2, axDist), [1 1], [0 nTargTrial], '-k', 'linewidth', 2)
        end
    end
end



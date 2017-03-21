function data = ccm_inhibition_ssd_metrics(subject, sessionID, options)
%%
data = [];

% Set default options or return a default options structure
if nargin < 3
    options.collapseSignal   	= false;
    options.collapseTarg        = true;
    options.include50           = false;
    options.USE_PRE_SSD         = true;
    options.USE_TWO_COLORS         = false;
    
    options.plotFlag            = true;
    options.printPlot           = true;
    options.figureHandle      	= 333;
    
    % Return just the default options struct if no input
    if nargin == 0
        Data           = options;
        return
    end
end

plotFlag        = options.plotFlag;
printPlot       = options.printPlot;




opt = ccm_options;
opt.plotFlag = 0;
opt.printPlot = 0;
DataInh = ccm_inhibition_rt(subject, sessionID, opt);


nSignal = length(DataInh.pSignalArray);
nSsd = length(DataInh.ssdArray);





% Set up plot
if plotFlag
    cMap = ccm_colormap(DataInh.pSignalArray);
    nCol = nSignal;
    nRow = 3;
    figureSSD = 27;
    ssdSpace = 30;  % Padding for the x-axis limits when plotting SSDs
    makerSize = 30;
    
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nCol, figureSSD);
    clf
    
    nCol = nSsd;
    nRow = 3;
    figureSSRT = 28;
    cohSpace = .02;  % Padding for the x-axis limits when plotting SSDs
    makerSizeSSRT = 30;
    
    [axisWidthSSRT, axisHeightSSRT, xAxesPositionSSRT, yAxesPositionSSRT] = standard_landscape(nRow, nCol, figureSSRT);
    clf
end




for iCoh = 1 : nSignal
    figure(figureSSD)
    % SSD vs. SSRT
    ax(1, iCoh) = axes('units', 'centimeters', 'position', [xAxesPosition(1, iCoh) yAxesPosition(1, iCoh) axisWidth axisHeight]);
    hold(ax(1, iCoh), 'on')
    set(ax(1, iCoh), 'ylim', [-100 400], 'xlim', [DataInh.ssdArray(1)-ssdSpace DataInh.ssdArray(end)+ssdSpace])
    grid on
    
    plot(ax(1, iCoh), DataInh.ssd{iCoh}, DataInh.ssrtIntegration{iCoh}, '.', 'color', cMap(iCoh,:), 'markerSize', makerSize)
    plot(ax(1, iCoh), [DataInh.ssdArray(1) DataInh.ssdArray(end)], [0 0], '--k')
    
    
    
    % SSD vs. Predicted-Observed noncanceled RTs
    ax(2, iCoh) = axes('units', 'centimeters', 'position', [xAxesPosition(2, iCoh) yAxesPosition(2, iCoh) axisWidth axisHeight]);
    hold(ax(2, iCoh), 'on')
    set(ax(2, iCoh), 'ylim', [-500 500], 'xlim', [DataInh.ssdArray(1)-ssdSpace DataInh.ssdArray(end)+ssdSpace])
    grid on
    
    plot(ax(2, iCoh), DataInh.ssdArray, DataInh.stopRespondRT(iCoh,:) - DataInh.stopRespondRTPredict(iCoh,:), '.', 'color', cMap(iCoh,:), 'markerSize', makerSize)
    plot(ax(2, iCoh), [DataInh.ssdArray(1) DataInh.ssdArray(end)], [0 0], '--k')
    
    
    
    % SSD vs. mean(noncanceled RTs) - SSDs
    ax(3, iCoh) = axes('units', 'centimeters', 'position', [xAxesPosition(3, iCoh) yAxesPosition(3, iCoh) axisWidth axisHeight]);
    hold(ax(3, iCoh), 'on')
    set(ax(3, iCoh), 'ylim', [-400 400], 'xlim', [DataInh.ssdArray(1)-ssdSpace DataInh.ssdArray(end)+ssdSpace])
    grid on
    
    plot(ax(3, iCoh), DataInh.ssdArray, DataInh.stopRespondRT(iCoh,:) - DataInh.ssdArray', '.', 'color', cMap(iCoh,:), 'markerSize', makerSize)
    plot(ax(3, iCoh), [DataInh.ssdArray(1) DataInh.ssdArray(end)], [0 0], '--k')
    
    
    if iCoh == 1
        ylabel(ax(1, iCoh), 'SSRT')
        ylabel(ax(2, iCoh), 'Obs - Predict: StopRTs')
        ylabel(ax(3, iCoh), 'StopRTs - SSD')
    end
    if iCoh > 1
        set(ax(1, iCoh), 'yticklabel', [])
        set(ax(2, iCoh), 'yticklabel', [])
        set(ax(3, iCoh), 'yticklabel', [])
    end
    
end

















for iSsd = 1 : nSsd
    figure(figureSSRT)
    
    % SSD vs. SSRT
    ax(1, iSsd) = axes('units', 'centimeters', 'position', [xAxesPositionSSRT(1, iSsd) yAxesPositionSSRT(1, iSsd) axisWidthSSRT axisHeightSSRT]);
    hold(ax(1, iSsd), 'on')
    set(ax(1, iSsd), 'ylim', [-100 400], 'xlim', [DataInh.pSignalArray(1)-cohSpace DataInh.pSignalArray(end)+cohSpace], 'xtick', DataInh.pSignalArray, 'xtickLabel', 100*DataInh.pSignalArray)
    grid on
    
    
    for jCoh = 1 : length(DataInh.pSignalArray)
        if ismember(DataInh.ssdArray(iSsd), DataInh.ssd{jCoh})
            jSsdInd = ismember(DataInh.ssd{jCoh}, DataInh.ssdArray(iSsd));
            plot(ax(1, iSsd), DataInh.pSignalArray(jCoh), DataInh.ssrtIntegration{jCoh}(jSsdInd), '.', 'color', cMap(jCoh,:), 'markerSize', makerSizeSSRT)
        end
    end
    
    if iSsd == 1
        ylabel(ax(1, iSsd), 'SSRT')
%         ylabel(ax(2, iSsd), 'Obs - Predict: StopRTs')
%         ylabel(ax(3, iSsd), 'StopRTs - SSD')
    end
    if iSsd > 1
        set(ax(1, iSsd), 'yticklabel', [])
%         set(ax(2, iSsd), 'yticklabel', [])
%         set(ax(3, iSsd), 'yticklabel', [])
    end
end












figure(figureSSD)
h=axes('Position', [0 0 1 1], 'Visible', 'Off');
titleString = sprintf('\n\n%s \t %s', subject,sessionID);
text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')

figure(figureSSRT)
h=axes('Position', [0 0 1 1], 'Visible', 'Off');
titleString = sprintf('\n\n%s \t %s', subject,sessionID);
text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')


if printPlot
    print(figureSSD,fullfile(local_figure_path, subject, [sessionID,'_ccm_ssd_metrics.pdf']),'-dpdf', '-r300')
    print(figureSSRT,fullfile(local_figure_path, subject, [sessionID,'_ccm_ssrt_metrics.pdf']),'-dpdf', '-r300')
end

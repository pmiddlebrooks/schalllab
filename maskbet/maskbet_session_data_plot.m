function maskbet_session_data_plot(Data, options)


nUnit = length(Data);
epochArray = options.epochArray;
figureHandle = options.figureHandle;
printPlot= options.printPlot;

% Which SOAs to analyze
if strcmp(options.soa, 'collapse')
    nSOA = 1;
    soaArray = 'collapse';
elseif strcmp(options.soa, 'each')
    nSOA = 4;
    soaArray = options.soaArray;
elseif strcmp(options.soa, 'mem')
    nSOA = 1;
    soaArray = options.soaArray(end);
else
    nSOA = 1;
    %     nSOA = length(options.soa);
    soaArray = options.soa;
end

% If collapsing into all left and all right (for decision stage),
% need to note here that there are "2" angles to deal with
% (important for calling maskbet_trial_selection.m)
switch options.decCollapseDir
    case {'leftRight','upDown'}
        nDecAngle = 2;
    case {'all'}
        nDecAngle = 1;
    otherwise
        nDecAngle = length(options.maskAngle);
end
switch options.betCollapseDir
    case {'leftRight','upDown'}
        nBetAngle = 2;
    case {'all'}
        nBetAngle = 1;
    otherwise
        nBetAngle = length(options.betAngle);
end




nEpoch = length(epochArray);


targLineW = 2;
distLineW = 1;
for kData = 1 : nUnit
    nRow = 4;
    nColumn = nEpoch;
    figureHandle = figureHandle + 1;
    if printPlot
        [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
    else
        [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
    end
    clf
    
%     kMax = options.signalMax(kData,:,:,:);
%     yLimMax = max(kMax(:));
    yLimMax = max(options.signalMax(kData,:));
    yLimMin = 0;
%     yLimMax = Data(kData).yMax * 1.1;
%     yLimMin = min(Data(kData).yMin * 1.1);
    
    %    yLimMin = 0;
    %    yLimMax = 65;
    for jEpoch = 1 : nEpoch
        
        jEpochName = epochArray{jEpoch};
        epochRange = maskbet_epoch_range(jEpochName, 'plot');
        
        
        
        if jEpoch <=4
            
            % **************************************************
            % Decision Stage data
            % **************************************************
            for m = 1 : nDecAngle
                
                
                
                % _______  Set up axes  ___________
                % Set up plot axes
                % Left target Go trials
                ax(m, jEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(m, jEpoch) yAxesPosition(m, jEpoch) axisWidth axisHeight]);
                set(ax(m, jEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
                cla
                hold(ax(m, jEpoch), 'on')
                plot(ax(m, jEpoch), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
                
                if m == 1
                title(epochArray{jEpoch})
                end
                if jEpoch > 1
                    set(ax(m, jEpoch), 'yticklabel', [])
                    set(ax(m, jEpoch), 'ycolor', [1 1 1])
                end
                
                
                
                
                
                for iSOA = 1 : nSOA
                    
                    plot_data
                    
                end  % SOA
                
                if jEpoch == 1
                    hAngle = sprintf('Angle: %s',num2str(options.maskAngle(m)));
%                     set(ax(m, jEpoch), 'ylabel', hAngle)
                    ylabel(hAngle)
                end
            end % nDecAnlge
            
            
        else
            % **************************************************
            % Bet Stage data
            % **************************************************
            for m = 1 : nBetAngle
                
                % _______  Set up axes  ___________
                % Set up plot axes
                % Left target Go trials
                ax(m, jEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(m, jEpoch) yAxesPosition(m, jEpoch) axisWidth axisHeight]);
                set(ax(m, jEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
                cla
                hold(ax(m, jEpoch), 'on')
                plot(ax(m, jEpoch), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
                
                if m == 1
                title(epochArray{jEpoch})
                end
                if jEpoch > 1
                    set(ax(m, jEpoch), 'yticklabel', [])
                    set(ax(m, jEpoch), 'ycolor', [1 1 1])
                end
                
                
                
                
                
                for iSOA = 1 : nSOA
                    
                    
                    plot_data
                    
                end  % SOA
                
            end % nBetAngle
        end % if m <= 4
    end % Epoch
    
end % Unit


    function plot_data
        % Combine the data from various outcomes if needed
        if strcmp(options.decOutcome, 'collapse') && strcmp(options.betOutcome, 'collapse')
            
            align = Data(kData).soa(iSOA).(jEpochName).angle(m).all.alignTime;
            if ~isempty(align)
                signal = nanmean(Data(kData).soa(iSOA).(jEpochName).angle(m).all.sdf, 1);
                plot(ax(m, jEpoch), epochRange, signal(align + epochRange), 'color', 'k', 'linewidth', targLineW)
            end
            
        elseif strcmp(options.decOutcome, 'each') && strcmp(options.betOutcome, 'collapse')
            
            % Target
            align = Data(kData).soa(iSOA).(jEpochName).angle(m).targ.alignTime;
            if ~isempty(align)
                signal = nanmean(Data(kData).soa(iSOA).(jEpochName).angle(m).targ.sdf, 1);
                plot(ax(m, jEpoch), epochRange, signal(align + epochRange), 'color', 'b', 'linewidth', targLineW)
            end
            
            % Distractor
            align = Data(kData).soa(iSOA).(jEpochName).angle(m).dist.alignTime;
            if ~isempty(align)
                signal = nanmean(Data(kData).soa(iSOA).(jEpochName).angle(m).dist.sdf, 1);
                plot(ax(m, jEpoch), epochRange, signal(align + epochRange), 'color', 'r', 'linewidth', targLineW)
            end
            
        elseif strcmp(options.decOutcome, 'collapse') && strcmp(options.betOutcome, 'each')
            
            % High
            align = Data(kData).soa(iSOA).(jEpochName).angle(m).high.alignTime;
            if ~isempty(align)
                signal = nanmean(Data(kData).soa(iSOA).(jEpochName).angle(m).high.sdf, 1);
                plot(ax(m, jEpoch), epochRange, signal(align + epochRange), 'color', 'b', 'linewidth', targLineW)
            end
            
            % Low
            align = Data(kData).soa(iSOA).(jEpochName).angle(m).low.alignTime;
            if ~isempty(align)
                signal = nanmean(Data(kData).soa(iSOA).(jEpochName).angle(m).low.sdf, 1);
                plot(ax(m, jEpoch), epochRange, signal(align + epochRange), 'color', 'r', 'linewidth', targLineW)
            end
            
        elseif strcmp(options.decOutcome, 'each') && strcmp(options.betOutcome, 'each')
            
            % Target - High
            align = Data(kData).soa(iSOA).(jEpochName).angle(m).targHigh.alignTime;
            if ~isempty(align)
                signal = nanmean(Data(kData).soa(iSOA).(jEpochName).angle(m).targHigh.sdf, 1);
                plot(ax(m, jEpoch), epochRange, signal(align + epochRange), 'color', 'b', 'linewidth', targLineW)
            end
            
            % Target - Low
            align = Data(kData).soa(iSOA).(jEpochName).angle(m).targLow.alignTime;
            if ~isempty(align)
                signal = nanmean(Data(kData).soa(iSOA).(jEpochName).angle(m).targLow.sdf, 1);
                plot(ax(m, jEpoch), epochRange, signal(align + epochRange), 'color', 'g', 'linewidth', targLineW)
            end
            
            
            % Distractor - High
            align = Data(kData).soa(iSOA).(jEpochName).angle(m).distHigh.alignTime;
            if ~isempty(align)
                signal = nanmean(Data(kData).soa(iSOA).(jEpochName).angle(m).distHigh.sdf, 1);
                plot(ax(m, jEpoch), epochRange, signal(align + epochRange), 'color', 'r', 'linewidth', targLineW)
            end
            
            % Distractor - Low
            align = Data(kData).soa(iSOA).(jEpochName).angle(m).distLow.alignTime;
            if ~isempty(align)
                signal = nanmean(Data(kData).soa(iSOA).(jEpochName).angle(m).distLow.sdf, 1);
                plot(ax(m, jEpoch), epochRange, signal(align + epochRange), 'color', 'y', 'linewidth', targLineW)
            end
            
        end
    end

end





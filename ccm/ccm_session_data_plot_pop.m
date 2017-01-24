function ccm_session_data_plot_pop(Data, Opt)

%%
% Set defaults


pSignalArray = Data(1).pSignalArray;
ssdArray    = Data(1).ssdArray;
sessionID   = Data(1).sessionID;
subjectID   = Data(1).subjectID;

dataType    = Opt.dataType;
printPlot   = Opt.printPlot;
filterData  = Opt.filterData;
stopHz      = Opt.stopHz;
figureHandle = Opt.figureHandle;
collapseSignal  = Opt.collapseSignal;
doStops     = Opt.doStops;

PLOT_ERROR  = false;

% epochArray = {'targOn', 'checkerOn', 'stopSignalOn', 'responseOnset', 'rewardOn'};
epochArray = {'targOn', 'checkerOn', 'stopSignalOn', 'responseOnset', 'toneOn'};


[nUnit, nTargPair] = size(Data);

nSignal = length(pSignalArray);
% If collapsing data across signal strength, adjust the pSignalArray here
if collapseSignal
    nSignal = 2;
end
easyLeftInd = 1;
hardLeftInd = 2;
hardRightInd = 3;
easyRightInd = 4;


nRow = 6;
nEpoch = length(epochArray);
nColumn = nEpoch * 2 + 1;
if ~strcmp(Opt.howProcess, 'step') && ~strcmp(Opt.howProcess,'print')
    figureHandle = figureHandle + 1;
end
figureHandle = 102;
% if printPlot
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
% else
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
% end
clf

yLimMin = 0;
yLimMax = nUnit;


for mEpoch = 1 : nEpoch
    mEpochName = epochArray{mEpoch};
    epochRange = ccm_epoch_range(mEpochName, 'plot');
    
    % _______  Set up axes  ___________
    % axes names
    axGoEasy = 1;
    axGoHard = 2;
    axStopGoEasy = 3;
    axStopGoHard = 4;
    axStopStopEasy = 5;
    axStopStopHard = 6;
    nRow = 6;
    
    alignLineX = [-epochRange(1) -epochRange(1)];
    % Set up plot axes
    
    % Left axes
    for i = 1 : nRow
        ax(i, mEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(i, mEpoch) yAxesPosition(i, mEpoch) axisWidth axisHeight]);
        set(ax(i, mEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [1  length(epochRange)], 'xticklabel', [])
        cla
        hold(ax(i, mEpoch), 'on')
        if i == 1
            title(epochArray{mEpoch})
        end
        
        % Right axes
        ax(i, mEpoch+nEpoch) = axes('units', 'centimeters', 'position', [xAxesPosition(i, mEpoch+nEpoch+1) yAxesPosition(i, mEpoch+nEpoch+1) axisWidth axisHeight]);
        set(ax(i, mEpoch+nEpoch), 'ylim', [yLimMin yLimMax], 'xlim', [1  length(epochRange)], 'xticklabel', [])
        cla
        hold(ax(i, mEpoch+nEpoch), 'on')
        if i == 1
            title(epochArray{mEpoch})
        end
        
        if mEpoch > 1
            set(ax(i, mEpoch), 'yticklabel', [])
            set(ax(i, mEpoch+nEpoch), 'yticklabel', [])
            set(ax(i, mEpoch), 'ycolor', [1 1 1])
            set(ax(i, mEpoch+nEpoch), 'ycolor', [1 1 1])
        end
    end
    
    
    
    
    
    % __________ Loop signal strengths and plot  _________
    
    
    
    
    
    
    % PLOT LEFT TARGET TRIALS
    
    
    
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
        for iCoh = 1 : 4
            iSdf = [];
            %             alignGoTarg = Data(1).(mEpochName).colorCoh(iCoh).goTarg.alignTime;
            %                alignGoDist = Data(kDataIndex, jTarg).(mEpochName).colorCoh(iPropIndexL).goDist.alignTime;
            
            for kUnit = 1:nUnit
                kNorm = Data(kUnit).yMax;
                kSdf = Data(kUnit).(mEpochName).colorCoh(iCoh).goTarg.(dataSignal) / kNorm;
                iSdf = [iSdf; kSdf];
            end
            if iCoh < 3
                iAx = iCoh;
                iEpoch = mEpoch;
            else
                iAx = iCoh - 2;
                iEpoch = mEpoch + nEpoch;
            end
            imagesc(ax(iAx, iEpoch), iSdf)
            plot(ax(iAx, iEpoch), alignLineX, [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
caxis(ax(iAx, iEpoch), [0 1])
            %                if ~isempty(alignGoDist) && PLOT_ERROR
            %                   sigGoDist = Data(kDataIndex, jTarg).(mEpochName).colorCoh(iPropIndexL).goDist.(dataSignal);
            %                   plot(ax(axGoEasy, mEpoch), epochRange, sigGoDist(alignGoDist + epochRange), '--', 'color', cMap(iPropIndexL,:), 'linewidth', distLineW)
            %                end
        end
    end
    
    
    
    
    
    
    % Stop trials
    if doStops
        switch dataType
            case 'neuron'
                dataSignal = 'sdf';
            case 'lfp'
                dataSignal = 'lfp';
            case 'erp'
                dataSignal = 'eeg';
        end
        
        for iCoh = 1 : 4
            iSdfNC = [];
            iSdfC = [];
            
            for kUnit = 1:nUnit
                kNorm = Data(kUnit).yMax;
                
                stopTargSig = cell(1, length(ssdArray));
                stopTargAlign = cell(1, length(ssdArray));
                %                             stopDistSig = cell(1, length(ssdArray));
                %                             stopDistAlign = cell(1, length(ssdArray));
                stopStopSig = cell(1, length(ssdArray));
                stopStopAlign = cell(1, length(ssdArray));
                for jSSDIndex = 1 : length(ssdArray)
                    stopTargSig{jSSDIndex} = Data(kUnit).(mEpochName).colorCoh(iCoh).stopTarg.ssd(jSSDIndex).(dataSignal);
                    stopTargAlign{jSSDIndex} = 1;
                    
                    %                 stopDistSig{jSSDIndex} = Data(kUnit).(mEpochName).colorCoh(iCoh).stopDist.ssd(jSSDIndex).(dataSignal);
                    %                 stopDistAlign{jSSDIndex} = 1;
                    %
                    if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
                        stopStopSig{jSSDIndex} = Data(kUnit).(mEpochName).colorCoh(iCoh).stopStop.ssd(jSSDIndex).(dataSignal);
                        stopStopAlign{jSSDIndex} = 1;
                    end
                    
                end  % jSSDIndex = 1 : length(ssdArray)
                
                [rasStopTarg, alignStopTarg] = align_raster_sets(stopTargSig, stopTargAlign);
                %             [rasStopDist, alignStopDist] = align_raster_sets(stopDistSig, stopDistAlign);
                switch dataType
                    case 'neuron'
                        sigStopTarg = nanmean(rasStopTarg, 1);
                        kSdfNC = sigStopTarg / kNorm;
                        iSdfNC = [iSdfNC; kSdfNC];
                        %                     sigStopDist = nanmean(rasStopDist, 1);
                    case {'lfp','erp'}
                        if filterData
                            sigStopTarg = lowpass(nanmean(rasStopTarg, 1)', stopHz)';
                            %                         sigStopDist = lowpass(nanmean(rasStopDist, 1)', stopHz)';
                        else
                            sigStopTarg = nanmean(rasStopTarg, 1);
                            %                         sigStopDist = nanmean(rasStopDist, 1);
                        end
                end
                if size(sigStopTarg, 2) == 1, sigStopTarg = []; end;
                %             if size(sigStopDist, 2) == 1, sigStopDist = []; end;
                
                if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
                    [rasStopCorrect, alignStopCorrect] = align_raster_sets(stopStopSig, stopStopAlign);
                    switch dataType
                        case 'neuron'
                            sigStopCorrect = nanmean(rasStopCorrect, 1);
                            kSdfC = sigStopCorrect / kNorm;
                            iSdfC = [iSdfC; kSdfC];
                        case {'lfp','erp'}
                            if filterData
                                sigStopCorrect = lowpass(nanmean(rasStopCorrect, 1)', stopHz)';
                            else
                                sigStopCorrect = nanmean(rasStopCorrect, 1);
                            end
                    end
                    if size(sigStopCorrect, 2) == 1, sigStopCorrect = []; end;
                end
                
            end
            
            if iCoh < 3
                iAxNC = iCoh + 2; % Noncanceled plots
                iAxC = iCoh + 4; % Canceled plots
                iEpoch = mEpoch;
            else
                iAxNC = iCoh;
                iAxC = iCoh + 2;
                iEpoch = mEpoch + nEpoch;
            end
            
            imagesc(ax(iAxNC, iEpoch), iSdfNC)
            plot(ax(iAxNC, iEpoch), alignLineX, [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
caxis(ax(iAxNC, iEpoch), [0 1])
                if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
            imagesc(ax(iAxC, iEpoch), iSdfC)
            plot(ax(iAxC, iEpoch), alignLineX, [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
caxis(ax(iAxC, iEpoch), [0 1])
                end
            
            %             if ~isempty(sigStopTarg)
            %                 plot(ax(axStopGoEasy, mEpoch), epochRange, sigStopTarg(alignStopTarg + epochRange), 'color', cMap(iPropIndexL,:), 'linewidth', targLineW)
            %             end
            %             if PLOT_ERROR && ~isempty(sigStopDist)
            %                 plot(ax(axStopGoEasy, mEpoch), epochRange, sigStopDist(alignStopDist + epochRange), '--', 'color', cMap(iPropIndexL,:), 'linewidth', distLineW)
            %             end
            %
            %             if ~strcmp(mEpochName, 'responseOnset')  % No stop signals on go trials
            %                 if ~isempty(sigStopCorrect)
            %                     plot(ax(axStopStopEasy, mEpoch), epochRange, sigStopCorrect(alignStopCorrect + epochRange), 'color', cMap(iPropIndexL,:), 'linewidth', targLineW)
            %                 end
            %             end
            
        end
    end % if doStops
    
    
    
    
end % mEpoch
colormap('jet')

%                             legend(ax(axGo, 1), {num2cell(pSignalArray'), num2str(pSignalArray')})

%         colorbar('peer', ax(axGo, 1), 'location', 'west')
%         colorbar('peer', ax(axStopGo, 1), 'location', 'west')
% h=axes('Position', [0 0 1 1], 'Visible', 'Off');
titleString = sprintf('%s', sessionID);
text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')
% if printPlot && ~collapseSignal
%     print(figureHandle,fullfile(local_figure_path, subjectID, [sessionID, '_ccm_population.pdf']),'-dpdf', '-r300')
    print(figureHandle,fullfile(local_figure_path, subjectID, [sessionID, '_ccm_population']),'-depsc', '-r300')
% elseif printPlot && collapseSignal
%     print(figureHandle,fullfile(local_figure_path, subjectID, [sessionID, '_ccm_', Data(kDataIndex, jTarg).name, '_', dataType, '_collapse.pdf']),'-dpdf', '-r300')
%     % micalaFolder = '/Volumes/SchallLab/Users/Paul/micala/';
%     % print(figureHandle,[micalaFolder, sessionID, '_ccm_', Data(kDataIndex, jTarg).name, '_',dataType,'_collapse.pdf'],'-dpdf', '-r300')
% end




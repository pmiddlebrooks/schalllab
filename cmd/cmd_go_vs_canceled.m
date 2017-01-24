function Data = cmd_go_vs_canceled(subjectID, sessionID, options)

% ccm_neuron_stop_vs_go(subjectID, sessionID, options)
%
% Compares noncanceled stops trials vs. latency matched (fast) go trials and canceled stop trials vs. latency matched (slower) go trials.
%
% If called without any arguments, returns a default options structure.
% If options are input but one is not specified, it assumes default.
%
% input:
%   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%
% Possible options are (default listed first):
%     options.dataType        	= which recorded signal to analyze:
%           'neuron','lfp',eeg'
%     options.unitArray       	= which units to analyze:
%           'each', <list of units in a cell array: e.g.:{'spikeUnit17a', spikeUnit17b'}>
%     options.collapseTarg 	= collapse angle/directions of the CORRECT TARGET
%           'none','leftRight','upDown','all'. 'None' analyzes each target
%           direction separately. E.g. 'leftRight' combines all left
%           targets and compares it with all right targets.
%     options.latencyMatchMethod  = wich method to use to match go trial latencies with canceled and noncanceled stop trials:
%           'ssrt','match','mean';
%     options.minTrialPerCond  	= how many trials must a condition have to
%           include in the analyses?
%
%     options.plotFlag       = true, false;
%     options.printPlot       = false, true;
%     options.figureHandle  = optional way to assign the figure to a handle


% Set default options or return a default options structure
if nargin < 3
    % Data type to collect/analyze
    options.dataType        	= 'neuron';
    
    options.unitArray          	= 'each';
    options.collapseTarg        = 'none';
    options.latencyMatchMethod 	= 'ssrt';
    options.minTrialPerCond     = 10;
    options.cellType            =      'move';
    options.ssrt                =      [];
    options.Unit                =      [];
   options.filterData       = false;
    
    options.plotFlag            = true;
    options.printPlot           = false;
    options.figureHandle      	= 420;
    
    % Return just the default options struct if no input
    if nargin == 0
        Data           = options;
        return
    end
end


filterData          = options.filterData;
collapseTarg    	= options.collapseTarg;
latencyMatchMethod  = options.latencyMatchMethod;
minTrialPerCond     = options.minTrialPerCond;

plotFlag            = options.plotFlag;
printPlot           = options.printPlot;
figureHandle        = options.figureHandle;

Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;

spikeWindow = -20 : 20;

epochName       = 'targOn';
eventMarkName   = 'responseOnset';
encodeTime = 10;

epochRangeTarg      = 1 : 400;
epochRangeSacc      = -199 : 200;


%% Get neural data from the session/unit:
if isempty(options.Unit)
    optSess = cmd_session_data;
    optSess.plotFlag = false;
    Unit = cmd_session_data(subjectID, sessionID, optSess);
else
    Unit = options.Unit;
end

if isempty(Unit)
    fprintf('Session %s does not contain spike data \n', sessionID)
    return
end

% How many units were recorded?
nUnit           = size(Unit, 2);
targAngleArray  = Unit(1).targAngleArray;
nAngle          = length(targAngleArray);
ssdArray        = Unit(1).ssdArray;
nSSD            = length(ssdArray);
unitArray     	= Unit(1).unitArray;





%% Get inhibition data from the session:
        if isempty(options.ssrt)
dataInh = cmd_inhibition(subjectID, sessionID, 'plotFlag', false);
        end




%%  Loop through Units and target pairs to collect and plot data


for kUnitIndex = 1 : nUnit
    
    
    for iTarg = 1 : nAngle
        iTarget = targAngleArray(iTarg);
        
        
        stopStopSpike           	= cell(nSSD, 1);
        stopStopTargRas          = cell(nSSD, 1);
        stopStopTargSDF          = cell(nSSD, 1);
        stopStopTargEventLat  	= cell(nSSD, 1);
        stopStopTargAlign      	= cell(nSSD, 1);
%         stopStopSaccRas             = cell(nSSD, 1);
%         stopStopSaccSDF             = cell(nSSD, 1);
%         stopStopSaccEventLat        = cell(nSSD, 1);
%         stopStopSaccAlign           = cell(nSSD, 1);
        
        goTargSlowSpike             = cell(nSSD, 1);
        goTargSlowTargRas        = cell(nSSD, 1);
        goTargSlowTargSDF        = cell(nSSD, 1);
        goTargSlowTargEventLat  	= cell(nSSD, 1);
        goTargSlowTargAlign    	= cell(nSSD, 1);
        goTargSlowSaccRas        	= cell(nSSD, 1);
        goTargSlowSaccSDF           = cell(nSSD, 1);
        goTargSlowSaccEventLat      = cell(nSSD, 1);
        goTargSlowSaccAlign         = cell(nSSD, 1);
        
        
        % For now, use the grand SSRT via integratin method
        if isempty(options.ssrt)
            ssrt = round(mean(dataInh.ssrtIntegrationWeightedAngle(iTarg)));
        else
            ssrt = options.ssrt;
        end
        
        
        usableCondition = zeros(nSSD, 1);
        
        
        
        
        % Get the go trial data: these need to be split to latency-match with
        % the stop trial data
        jGoTargTarg      = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'goTarg'}, iTarget);
        jGoTargSacc      = cmd_concat_neural_conditions(Unit(kUnitIndex), 'responseOnset', 'targOn', {'goTarg'}, iTarget);
        
        
        
        for k = 1 : nSSD
            kSSD = ssdArray(k);
            
            
            
            % Get the stop trial data
            jStopStopTarg    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, 'fixOn', {'stopStop'}, iTarget, kSSD);
%             jStopTargSacc    = cmd_concat_neural_conditions(Unit(kUnitIndex), 'responseOnset', 'targOn', {'stopStop'}, iTarget, kSSD);
            jStopTargTarg    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'stopStop'}, iTarget, kSSD);
            
            if size(jStopStopTarg.eventLatency, 1) > 0
                switch latencyMatchMethod
                    case 'ssrt'
                        jStopLatency = kSSD + ssrt;
                        jGoSlowTrial = jGoTargTarg.eventLatency > jStopLatency;
                    case 'mean'
                        jGoTargRT = sort(jGoTargTarg.eventLatency);
                        meanStopTargRT = nanmean(jStopTargTarg.eventLatency);
                        while nanmean(jGoTargRT) > meanStopTargRT
                            jGoTargRT(end) = [];
                        end
                        jStopLatency = jGoTargRT(end);
                        jGoSlowTrial = jGoTargTarg.eventLatency <= jStopLatency;
                    case 'match'
                        % Use nearest neighbor method to get rt-matched
                        % trials
                        nStopCorrect = size(jStopStopTarg.raster, 1);
                        data = ccm_match_rt(jGoTargTarg.eventLatency, jStopTargTarg.eventLatency(jStopTargTrial), nStopCorrect);
                        jGoSlowTrial = data.goFastTrial;
                end
                
                
                if sum(jGoSlowTrial)
                    
                    
                    stopStopTargRas{k}          = jStopStopTarg.raster;
                    jStopStopTargSDF            = nanmean(spike_density_function(stopStopTargRas{k}, Kernel), 1);
                    stopStopTargSDF{k}          = jStopStopTargSDF;
                    stopStopTargAlign{k}        = jStopStopTarg.align;
                    stopStopTargEventLat{k}     = jStopStopTarg.eventLatency;
%                     stopStopSaccRas{k}          = jStopTargSacc.raster(jStopTargTrial,:);
%                     jStopTargSaccSDF            = nanmean(spike_density_function(stopStopSaccRas{k}, Kernel), 1);
%                     stopStopSaccSDF{k}          = jStopTargSaccSDF;
%                     stopStopSaccAlign{k}        = jStopTargSacc.align;
%                     stopStopSaccEventLat{k} 	= jStopTargSacc.eventLatency(jStopTargTrial,:);
                    
                    goTargSlowTargRas{k}        = jGoTargTarg.raster(jGoSlowTrial,:);
                    jGoSlowTargSDF          	= nanmean(spike_density_function(goTargSlowTargRas{k}, Kernel), 1);
                    goTargSlowTargSDF{k}        = jGoSlowTargSDF;
                    goTargSlowTargAlign{k}      = jGoTargTarg.align;
                    goTargSlowTargEventLat{k}	= jGoTargTarg.eventLatency(jGoSlowTrial,:);
                    goTargSlowSaccRas{k}      	= jGoTargSacc.raster(jGoSlowTrial,:);
                    jGoSlowSaccSDF                              = nanmean(spike_density_function(goTargSlowSaccRas{k}, Kernel), 1);
                    goTargSlowSaccSDF{k}         = jGoSlowSaccSDF;
                    goTargSlowSaccAlign{k}       = jGoTargSacc.align;
                    goTargSlowSaccEventLat{k}    = jGoTargSacc.eventLatency(jGoSlowTrial,:);
                    
                    
                         % TEST #1:
                   % Hanes et al 1998 t-test of spike rates 40 ms surrounding estimated ssrt
                    stopStopSpike{k}     = sum(jStopStopTarg.raster(:, spikeWindow + jStopStopTarg.align + kSSD + ssrt), 2);
                    goTargSlowSpike{k}   = sum(jGoTargTarg.raster(jGoSlowTrial, spikeWindow + jGoTargTarg.align + kSSD + ssrt), 2);
                    [h,p,ci,sts]                        = ttest2(stopStopSpike{k}, goTargSlowSpike{k});
                    
                    pValue{k}    = p;
                    stats{k}     = sts;
                    
                    
                    
                        
                        % TEST #2:
                        % Hanes et al 1998 (p.822) differential sdf test
                        % ------------------------------------------------------------------------
                        cancelTime = nan; % Initialize to NaN;
                        sdfDiff = jGoSlowTargSDF(jGoTargTarg.align + (-599 : epochRangeTarg(end)))' - jStopStopTargSDF(jStopStopTarg.align + (-599 : epochRangeTarg(end)))';
                        % If user thinks it's a fixation/stopping type cell, flip
                        % the sign of the differential sdf
                        if strcmp(options.cellType, 'fix')
                            sdfDiff = -sdfDiff;
                        end
                        sdfDiffCheckerOn = 600;
                        stdDiff = std(sdfDiff(1:sdfDiffCheckerOn));
                        
                        
                        % are there times at which the difference between sdfs is
                        % greater than 2 standard deviations of the difference 600
                        % ms before checkerboard onset?
                        std2Ind = sdfDiff(sdfDiffCheckerOn : end) > 2*stdDiff;
                        
                        % Look for a sequence of 50 ms for which the go sdf is 2
                        % std greater than the stop sdf.
                        % First whether the differential sdf was > 2*Std for the
                        % first 50 ms
                        pass50msTestInd = [];
                        if sum(std2Ind(1:50)) == 50
                            pass50msTestInd = 1;
                        else
                            % If it wasn't, determein whether there was a time
                            % after the checkerboard onset that the differential
                            % sdf was > 2*Std for at least 50 ms.
                            riseAbove2Std = find([0; diff(std2Ind)] == 1);
                            sinkBelow2Std = find([0; diff(std2Ind)] == -1);
                            if ~isempty(riseAbove2Std)
                                % Get rid of occasions for which the signals differ
                                % going into the epoch (and therefore they will
                                % cease to differ before they begin again to
                                % differ)
                                if ~isempty(sinkBelow2Std)
                                    sinkBelow2Std(sinkBelow2Std < riseAbove2Std(1)) = [];
                                end
                                % Now riseAbove2Std length shouldb equal to or 1
                                % more than sinkBelow2Std. If they're equal, see if
                                % any of the riseAbove2Std streaks go longer than
                                % 50ms
                                if length(riseAbove2Std) == length(sinkBelow2Std)
                                    ind = find(sinkBelow2Std - riseAbove2Std >= 50, 1);
                                    if ~isempty(ind)
                                        pass50msTestInd = riseAbove2Std(ind);
                                    end
                                    % If they're not equal, the last riseAbove2Std
                                    % will last until the end of the sdf: see if
                                    % that is at least 50 ms
                                elseif length(riseAbove2Std) > length(sinkBelow2Std)
                                    if riseAbove2Std(end) < length(std2Ind) - 49
                                        pass50msTestInd = riseAbove2Std(end);
                                    end
                                end
                            end
                        end
                        
                        
                        % If there werwe 50 consecutive ms of sdfDiff > 2*Std,
                        % check whether the difference ever reached 6*Std
                        if ~isempty(pass50msTestInd)
                            std6Ind = sdfDiff(pass50msTestInd : end) > 6*stdDiff;
                            if sum(std6Ind)
                                cancelTime = kSSD + ssrt - pass50msTestInd;
                            end
                        end
                        
                        
                        
                        
                        % Mark the data for later display if there were enough trials in both conditions
                    if sum(jGoSlowTrial) >= options.minTrialPerCond && size(jStopStopTarg.eventLatency, 1) >= options.minTrialPerCond
                        usableCondition(k) = 1;
                        
                            fprintf('ssd: %d  \tgo v. stop: %.2f  %.2f sp, p = %.2f\t canceltime = %d\n',...
                                kSSD, mean(goTargSlowSpike{k}), mean(stopStopSpike{k}), p, cancelTime);
                    end
                end % if sum(jGoSlowTrial)
            end % if sum(jStopTargTrial)
        end
        
        
        
        % ****************************  PLOTTING  ****************************
        
        % Figure out how many conditions to plot/analyze
        usableCondition = logical(usableCondition(:));
        ssdVector = reshape(ssdArray, nSSD, 1);
        ssdVector = ssdVector(usableCondition);
        
        
        nUsable = sum(usableCondition);
        
        
        
        
        
        
        
        
        
        
        
        
        
        % Creat a new figure if needed (for a new pair of targets and/or
        % unit
        if plotFlag && nUsable > 0
            
            cMap = cmd_colormap;
            
            goLineW = 2;
            stopLineW = 2;
            markSize = 20;
            nRow = max(3, nUsable);
            figureHandle = figureHandle + 1;
            nColumn = 2;
            if printPlot
                [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nRow, nColumn, 'portrait', figureHandle);
            else
                [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
            end
            axisHeight = axisHeight * .9;
            clf
            
            
            % Figure out y-axis limits (to be consistent across graphs)
            dataMx    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'goTarg', 'stopStop'}, targAngleArray(iTarg), ssdArray);
            sdfAll   = nanmean(spike_density_function(dataMx.raster(:,:), Kernel), 1);
            sdfMax = max(sdfAll .* 1.7);
            %             epochRange = ccm_epoch_range(epochName, 'plot');
            
            % Title for the figure
            h=axes('Position', [0 0 1 1], 'Visible', 'Off');
            titleString = sprintf('%s \tUnit %s:\tAngle: %d', sessionID, Unit(kUnitIndex).name, targAngleArray(iTarg));
            text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top', 'interpreter', 'none')
            
            
            
            for i = 1 : nUsable
                
                % which pSignalArray index and which ssdArray index?
                iSsdInd = find(ssdArray == ssdVector(i));
                iSSD = ssdArray(iSsdInd);
                
                
                % _______  Set up axes  ___________
                iRow = i;
                colTarg = 1;
                colSacc = 2;
                
                % Data aligned on checkerboard onset
                ax(iRow, colTarg) = axes('units', 'centimeters', 'position', [xAxesPosition(iRow, colTarg) yAxesPosition(iRow, colTarg) axisWidth axisHeight]);
                set(ax(iRow, colTarg), 'ylim', [0 sdfMax * 1.1], 'xlim', [epochRangeTarg(1) epochRangeTarg(end)])
                cla
                hold(ax(iRow, colTarg), 'on')
                plot(ax(iRow, colTarg), [1 1], [0 sdfMax], '-k', 'linewidth', 2)
                ttl = sprintf('SSD: %d ', iSSD);
                title(ttl)
                
                % Data aligned on response onset
                ax(iRow, colSacc) = axes('units', 'centimeters', 'position', [xAxesPosition(iRow, colSacc) yAxesPosition(iRow, colSacc) axisWidth axisHeight]);
                set(ax(iRow, colSacc), 'ylim', [0 sdfMax * 1.1], 'xlim', [epochRangeSacc(1) epochRangeSacc(end)])
                set(ax(iRow, colSacc), 'yticklabel', [], 'ycolor', [1 1 1])
                cla
                hold(ax(iRow, colSacc), 'on')
                plot(ax(iRow, colSacc), [1 1], [0 sdfMax], '-k', 'linewidth', 2)
                
                
                %             if k > 1
                %                 set(ax(kRow, k), 'yticklabel', [])
                %                 set(ax(kRow, k), 'ycolor', [1 1 1])
                %             end
                
                
                
                iGoTargSlowTargSDF = goTargSlowTargSDF{iSsdInd}(goTargSlowTargAlign{iSsdInd} + epochRangeTarg);
                iStopStopTargSDF = stopStopTargSDF{iSsdInd}(stopStopTargAlign{iSsdInd} + epochRangeTarg);
                plot(ax(iRow, colTarg), [iSSD, iSSD], [0 sdfMax], 'color', [.2 .2 .2], 'linewidth', 1)
                plot(ax(iRow, colTarg), [iSSD + ssrt, iSSD + ssrt], [0 sdfMax], '--', 'color', [0 0 0], 'linewidth', 1)
                plot(ax(iRow, colTarg), epochRangeTarg, iGoTargSlowTargSDF, 'color', cMap.goTarg, 'linewidth', goLineW)
                iGoTargRTMean = round(mean(goTargSlowTargEventLat{iSsdInd}));
                plot(ax(iRow, colTarg), iGoTargRTMean, iGoTargSlowTargSDF(iGoTargRTMean), '.b','markersize', markSize)
                plot(ax(iRow, colTarg), epochRangeTarg, iStopStopTargSDF, 'color', cMap.stopTarg, 'linewidth', stopLineW)
                
                
                iGoTargSlowSaccSDF = goTargSlowSaccSDF{iSsdInd}(goTargSlowSaccAlign{iSsdInd} + epochRangeSacc);
                goTargSlowRTMean = round(mean(goTargSlowTargEventLat{iSsdInd}));
                iStopTargSaccSDF = stopStopTargSDF{iSsdInd}(stopStopTargAlign{iSsdInd} + goTargSlowRTMean + epochRangeSacc);
                plot(ax(iRow, colSacc), epochRangeSacc, iGoTargSlowSaccSDF, 'color', cMap.goTarg, 'linewidth', goLineW)
                plot(ax(iRow, colSacc), epochRangeSacc, iStopTargSaccSDF, 'color', cMap.stopTarg, 'linewidth', stopLineW)
                
                
                
                % plot(ax(kRow, k), [kSSD, kSSD], [0 sdfMax * .15], '--', 'color', [.2 .2 .2], 'linewidth', 1)
                %             plot(ax(kRow, k), [kSSD + ssrt, kSSD + ssrt], [0 sdfMax], 'color', [0 0 0], 'linewidth', 1)
                %             plot(ax(kRow, k), epochRangeTarg, jGoSlowSDF(jGoTargTarg.align + epochRangeTarg), 'color', cMap.goTarg, 'linewidth', goLineW)
                %             plot(ax(kRow, k), epochRangeTarg, jStopStopTarg.sdf(jStopStopTarg.align + epochRangeTarg), 'color', cMap.stopStop, 'linewidth', goLineW)
            end % for i = 1 : nUsable
            
            
        end % if plotFlag
            
            
            if printPlot
                localFigurePath = local_figure_path;
                print(figureHandle,[localFigurePath, sessionID, '_',Unit(kUnitIndex).name, '_', iTarg, '_cmd_go_vs_canceled.pdf'],'-dpdf', '-r300')
            end
        
        
        
        
                % Collect the data for later analyses
        Data(kUnitIndex).targ(iTarg).stopStopSpike      = stopStopSpike;
        Data(kUnitIndex).targ(iTarg).stopStopTargRas        = stopStopTargRas;
        Data(kUnitIndex).targ(iTarg).stopStopTargAlign        = stopStopTargAlign;
        Data(kUnitIndex).targ(iTarg).stopStopTargEventLat        = stopStopTargEventLat;
%         Data(kUnitIndex).targ(iTarg).stopStopSaccRas  	= stopStopSaccRas;
%         Data(kUnitIndex).targ(iTarg).stopStopSaccAlign        = stopStopSaccAlign;
%         Data(kUnitIndex).targ(iTarg).stopStopSaccEventLat        = stopStopSaccEventLat;
        
        Data(kUnitIndex).targ(iTarg).goTargSlowSpike    = goTargSlowSpike;
        Data(kUnitIndex).targ(iTarg).goTargSlowTargRas      = goTargSlowTargRas;
        Data(kUnitIndex).targ(iTarg).goTargSlowTargAlign      = goTargSlowTargAlign;
        Data(kUnitIndex).targ(iTarg).goTargSlowTargEventLat     = goTargSlowTargEventLat;
        Data(kUnitIndex).targ(iTarg).goTargSlowSaccRas 	= goTargSlowSaccRas;
        Data(kUnitIndex).targ(iTarg).goTargSlowSaccAlign      = goTargSlowSaccAlign;
        Data(kUnitIndex).targ(iTarg).goTargSlowSaccEventLat     = goTargSlowSaccEventLat;

    end % iTargInd
end % for kUnitIndex = 1 : nUnit


return
end

%         if usableStopStop(k)
%             ax(kRow, k) = axes('units', 'centimeters', 'position', [xAxesPosition(kRow, k) yAxesPosition(kRow, k) axisWidth axisHeight]);
%             set(ax(kRow, k), 'ylim', [0 sdfMax * 1.1], 'xlim', [epochRangeTarg(1) epochRangeTarg(end)])
%             cla
%             hold(ax(kRow, k), 'on')
%             plot(ax(kRow, k), [0 0], [0 sdfMax], '-k', 'linewidth', 2)
%             if ~usableStopTarg(k)
%                 ttl = sprintf('SSD:  %d', kSSD);
%                 title(ttl)
%             end
%
%             if k > 1
%                 set(ax(kRow, k), 'yticklabel', [])
%                 set(ax(kRow, k), 'ycolor', [1 1 1])
%             end
%
%             % Get the stop trial data
%             jStopStopTarg    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'stopStop'}, iTarget, kSSD);
%             jStopTargTarg    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'stopStop'}, iTarget, kSSD);
%
%             if isempty(jStopTargTarg.sdf) && ~strcmp(latencyMatchMethod, 'ssrt')
%                 continue
%             end
%             switch latencyMatchMethod
%                 case 'ssrt'
%                     jStopLatency = kSSD + ssrt;
%                     jGoSlowTrial = jGoTargTarg.eventLatency > jStopLatency;
%                 case 'mean'
%                     if isempty(jStopTargTarg.sdf)
%                         continue
%                     end
%                     jGoTargRT = sort(jGoTargTarg.eventLatency);
%                     meanStopTargRT = nanmean(jStopTargTarg.eventLatency);
%                     while nanmean(jGoTargRT) > meanStopTargRT
%                         jGoTargRT(end) = [];
%                     end
%                     jStopLatency = jGoTargRT(end);
%                     jGoSlowTrial = jGoTargTarg.eventLatency > jStopLatency;
%                 case 'match'
%                     if isempty(jStopTargTarg.sdf)
%                         continue
%                     end
%                     % Use nearest neighbor method to get rt-matched
%                     % trials
%                     nStopCorrect = size(jStopStopTarg.raster, 1);
%                     data = ccm_match_rt(jGoTargTarg.eventLatency, jStopTargTarg.eventLatency, nStopCorrect);
%                     jGoSlowTrial = data.goSlowTrial;
%             end
%
%             jGoSlowSDF   = nanmean(spike_density_function(jGoTargTarg.raster(jGoSlowTrial,:), Kernel), 1);
%
%             plot(ax(kRow, k), [kSSD, kSSD], [0 sdfMax * .15], '--', 'color', [.2 .2 .2], 'linewidth', 1)
%             plot(ax(kRow, k), [kSSD + ssrt, kSSD + ssrt], [0 sdfMax], 'color', [0 0 0], 'linewidth', 1)
%             plot(ax(kRow, k), epochRangeTarg, jGoSlowSDF(jGoTargTarg.align + epochRangeTarg), 'color', cMap.goTarg, 'linewidth', goLineW)
%             plot(ax(kRow, k), epochRangeTarg, jStopStopTarg.sdf(jStopStopTarg.align + epochRangeTarg), 'color', cMap.stopStop, 'linewidth', goLineW)
%         end





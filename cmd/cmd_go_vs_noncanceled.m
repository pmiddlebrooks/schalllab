function Data = cmd_go_vs_noncanceled(subjectID, sessionID, options)

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
    options.figureHandle      	= 400;
    
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
        
        
        stopTargSpike           	= cell(nSSD, 1);
        stopTargTargRas          = cell(nSSD, 1);
        stopTargTargSDF          = cell(nSSD, 1);
        stopTargTargEventLat  	= cell(nSSD, 1);
        stopTargTargAlign      	= cell(nSSD, 1);
        stopTargSaccRas             = cell(nSSD, 1);
        stopTargSaccSDF             = cell(nSSD, 1);
        stopTargSaccEventLat        = cell(nSSD, 1);
        stopTargSaccAlign           = cell(nSSD, 1);
        
        goTargFastSpike             = cell(nSSD, 1);
        goTargFastTargRas        = cell(nSSD, 1);
        goTargFastTargSDF        = cell(nSSD, 1);
        goTargFastTargEventLat  	= cell(nSSD, 1);
        goTargFastTargAlign    	= cell(nSSD, 1);
        goTargFastSaccRas        	= cell(nSSD, 1);
        goTargFastSaccSDF           = cell(nSSD, 1);
        goTargFastSaccEventLat      = cell(nSSD, 1);
        goTargFastSaccAlign         = cell(nSSD, 1);
        
        
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
            jStopTargTarg    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'stopTarg'}, iTarget, kSSD);
            jStopTargSacc    = cmd_concat_neural_conditions(Unit(kUnitIndex), 'responseOnset', 'targOn', {'stopTarg'}, iTarget, kSSD);
            jStopStopTarg    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'stopStop'}, iTarget, kSSD);
            
            % Use a subsample of the noncanceled stop RTs that are
            % later than the SSD plus some time to encode the stimulus
            jStopTargTrial = jStopTargTarg.eventLatency > kSSD + encodeTime;
            if sum(jStopTargTrial)
                switch latencyMatchMethod
                    case 'ssrt'
                        jStopLatency = kSSD + ssrt;
                        jGoFastTrial = jGoTargTarg.eventLatency <= jStopLatency & jGoTargTarg.eventLatency >  kSSD + encodeTime;
                    case 'mean'
                        jGoTargRT = sort(jGoTargTarg.eventLatency);
                        meanStopTargRT = nanmean(jStopTargTarg.eventLatency);
                        while nanmean(jGoTargRT) > meanStopTargRT
                            jGoTargRT(end) = [];
                        end
                        jStopLatency = jGoTargRT(end);
                        jGoFastTrial = jGoTargTarg.eventLatency <= jStopLatency;
                    case 'match'
                        % Use nearest neighbor method to get rt-matched
                        % trials
                        nStopCorrect = size(jStopStopTarg.raster, 1);
                        data = ccm_match_rt(jGoTargTarg.eventLatency, jStopTargTarg.eventLatency(jStopTargTrial), nStopCorrect);
                        jGoFastTrial = data.goFastTrial;
                end
                
                
                if sum(jGoFastTrial)
                    
                    
                    stopTargTargRas{k}          = jStopTargTarg.raster(jStopTargTrial,:);
                    jStopTargTargSDF            = nanmean(spike_density_function(stopTargTargRas{k}, Kernel), 1);
                    stopTargTargSDF{k}          = jStopTargTargSDF;
                    stopTargTargAlign{k}        = jStopTargTarg.align;
                    stopTargTargEventLat{k}     = jStopTargTarg.eventLatency(jStopTargTrial,:);
                    stopTargSaccRas{k}          = jStopTargSacc.raster(jStopTargTrial,:);
                    jStopTargSaccSDF            = nanmean(spike_density_function(stopTargSaccRas{k}, Kernel), 1);
                    stopTargSaccSDF{k}          = jStopTargSaccSDF;
                    stopTargSaccAlign{k}        = jStopTargSacc.align;
                    stopTargSaccEventLat{k} 	= jStopTargSacc.eventLatency(jStopTargTrial,:);
                    
                    goTargFastTargRas{k}        = jGoTargTarg.raster(jGoFastTrial,:);
                    jGoFastTargSDF          	= nanmean(spike_density_function(goTargFastTargRas{k}, Kernel), 1);
                    goTargFastTargSDF{k}        = jGoFastTargSDF;
                    goTargFastTargAlign{k}      = jGoTargTarg.align;
                    goTargFastTargEventLat{k}	= jGoTargTarg.eventLatency(jGoFastTrial,:);
                    goTargFastSaccRas{k}      	= jGoTargSacc.raster(jGoFastTrial,:);
                    jGoFastSaccSDF                              = nanmean(spike_density_function(goTargFastSaccRas{k}, Kernel), 1);
                    goTargFastSaccSDF{k}         = jGoFastSaccSDF;
                    goTargFastSaccAlign{k}       = jGoTargSacc.align;
                    goTargFastSaccEventLat{k}    = jGoTargSacc.eventLatency(jGoFastTrial,:);
                    
                    % Hanes et al 1998 t-test of spike rates 40 ms surrounding estimated ssrt
                    stopTargSpike{k}     = sum(jStopTargTarg.raster(jStopTargTrial, spikeWindow + jStopTargTarg.align + kSSD + ssrt), 2);
                    goTargFastSpike{k}   = sum(jGoTargTarg.raster(jGoFastTrial, spikeWindow + jGoTargTarg.align + kSSD + ssrt), 2);
                    [h,p,ci,sts]                        = ttest2(stopTargSpike{k}, goTargFastSpike{k});
                    
                    pValue{k}    = p;
                    stats{k}     = sts;
                    
                    
                    
                    % Mark the data for later display if there were enough trials in both conditions
                    if sum(jGoFastTrial) >= options.minTrialPerCond && sum(jStopTargTrial) >= options.minTrialPerCond
                        usableCondition(k) = 1;
                        
                    end
                end % if sum(jGoFastTrial)
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
            dataMx    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'goTarg', 'stopTarg'}, targAngleArray(iTarg), ssdArray);
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
                
                
                
                iGoTargFastTargSDF = goTargFastTargSDF{iSsdInd}(goTargFastTargAlign{iSsdInd} + epochRangeTarg);
                iStopTargTargSDF = stopTargTargSDF{iSsdInd}(stopTargTargAlign{iSsdInd} + epochRangeTarg);
                plot(ax(iRow, colTarg), [iSSD, iSSD], [0 sdfMax], 'color', [.2 .2 .2], 'linewidth', 1)
                plot(ax(iRow, colTarg), [iSSD + ssrt, iSSD + ssrt], [0 sdfMax], '--', 'color', [0 0 0], 'linewidth', 1)
                plot(ax(iRow, colTarg), epochRangeTarg, iGoTargFastTargSDF, 'color', cMap.goTarg, 'linewidth', goLineW)
                iGoTargRTMean = round(mean(goTargFastTargEventLat{ iSsdInd}));
                plot(ax(iRow, colTarg), iGoTargRTMean, iGoTargFastTargSDF(iGoTargRTMean), '.b','markersize', markSize)
                plot(ax(iRow, colTarg), epochRangeTarg, iStopTargTargSDF, 'color', cMap.stopTarg, 'linewidth', stopLineW)
                iStopTargRTMean = round(mean(stopTargTargEventLat{ iSsdInd}));
                plot(ax(iRow, colTarg), iStopTargRTMean, iStopTargTargSDF(iStopTargRTMean), '.b','markersize', markSize)
                
                
                iGoTargFastSaccSDF = goTargFastSaccSDF{iSsdInd}(goTargFastSaccAlign{iSsdInd} + epochRangeSacc);
                iStopTargSaccSDF = stopTargSaccSDF{iSsdInd}(stopTargSaccAlign{iSsdInd} + epochRangeSacc);
                plot(ax(iRow, colSacc), epochRangeSacc, iGoTargFastSaccSDF, 'color', cMap.goTarg, 'linewidth', goLineW)
                %                 iGoTargTargMean = round(mean(goTargFastSaccEventLat{ iSsdInd}));
                %                 plot(ax(iRow, colSacc), iGoTargTargMean, iStopTargSaccSDF(iGoTargTargMean), 'ok','markersize', 10)
                plot(ax(iRow, colSacc), epochRangeSacc, iStopTargSaccSDF, 'color', cMap.stopTarg, 'linewidth', stopLineW)
                %                 iGoTargTargMean = round(mean(stopTargSaccEventLat{ iSsdInd}));
                %                 plot(ax(iRow, colSacc), iGoTargTargMean, iStopTargSaccSDF(iGoTargTargMean), 'ok','markersize', 10)
                
                
                
                % plot(ax(kRow, k), [kSSD, kSSD], [0 sdfMax * .15], '--', 'color', [.2 .2 .2], 'linewidth', 1)
                %             plot(ax(kRow, k), [kSSD + ssrt, kSSD + ssrt], [0 sdfMax], 'color', [0 0 0], 'linewidth', 1)
                %             plot(ax(kRow, k), epochRangeTarg, jGoFastSDF(jGoTargTarg.align + epochRangeTarg), 'color', cMap.goTarg, 'linewidth', targLineW)
                %             plot(ax(kRow, k), epochRangeTarg, jStopTargTarg.sdf(jStopTargTarg.align + epochRangeTarg), 'color', cMap.stopTarg, 'linewidth', targLineW)
            end % for i = 1 : nUsable
            
            
        end % if plotFlag
            
            
            if printPlot
                localFigurePath = local_figure_path;
                print(figureHandle,[localFigurePath, sessionID, '_',Unit(kUnitIndex).name, '_', iTarg, '_cmd_go_vs_noncanceled.pdf'],'-dpdf', '-r300')
            end
        
        
        
        
                % Collect the data for later analyses
        Data(kUnitIndex).targ(iTarg).stopTargSpike      = stopTargSpike;
        Data(kUnitIndex).targ(iTarg).stopTargTargRas        = stopTargTargRas;
        Data(kUnitIndex).targ(iTarg).stopTargTargAlign        = stopTargTargAlign;
        Data(kUnitIndex).targ(iTarg).stopTargTargEventLat        = stopTargTargEventLat;
        Data(kUnitIndex).targ(iTarg).stopTargSaccRas  	= stopTargSaccRas;
        Data(kUnitIndex).targ(iTarg).stopTargSaccAlign        = stopTargSaccAlign;
        Data(kUnitIndex).targ(iTarg).stopTargSaccEventLat        = stopTargSaccEventLat;
        
        Data(kUnitIndex).targ(iTarg).goTargFastSpike    = goTargFastSpike;
        Data(kUnitIndex).targ(iTarg).goTargFastTargRas      = goTargFastTargRas;
        Data(kUnitIndex).targ(iTarg).goTargFastTargAlign      = goTargFastTargAlign;
        Data(kUnitIndex).targ(iTarg).goTargFastTargEventLat     = goTargFastTargEventLat;
        Data(kUnitIndex).targ(iTarg).goTargFastSaccRas 	= goTargFastSaccRas;
        Data(kUnitIndex).targ(iTarg).goTargFastSaccAlign      = goTargFastSaccAlign;
        Data(kUnitIndex).targ(iTarg).goTargFastSaccEventLat     = goTargFastSaccEventLat;

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
%             jStopTargTarg    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'stopTarg'}, iTarget, kSSD);
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
%             plot(ax(kRow, k), epochRangeTarg, jGoSlowSDF(jGoTargTarg.align + epochRangeTarg), 'color', cMap.goTarg, 'linewidth', targLineW)
%             plot(ax(kRow, k), epochRangeTarg, jStopStopTarg.sdf(jStopStopTarg.align + epochRangeTarg), 'color', cMap.stopStop, 'linewidth', targLineW)
%         end





function Data = cmd_neuron_stop_vs_go(subjectID, sessionID, options)

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
    options.minTrialPerCond     = 5;
    options.ssrt                = [];
    
    options.plotFlag            = true;
    options.printPlot           = false;
    options.figureHandle      	= 400;
    
    % Return just the default options struct if no input
    if nargin == 0
        Data           = options;
        return
    end
end


unitArray           = options.unitArray;
collapseTarg    	= options.collapseTarg;
latencyMatchMethod  = options.latencyMatchMethod;
minTrialPerCond     = options.minTrialPerCond;

plotFlag            = options.plotFlag;
printPlot           = options.printPlot;
figureHandle        = options.figureHandle;

Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;


%%
% Get neural data from the session/unit:
optSess = cmd_session_data;
optSess.plotFlag = false;
Unit = cmd_session_data(subjectID, sessionID, optSess);
if isempty(Unit)
    fprintf('Session %s does not contain spike data \n', sessionID)
    return
end

% How many units were recorded?
nUnit = size(Unit, 2);
targAngleArray = Unit(1).targAngleArray;
nAngle = length(targAngleArray);
ssdArray = Unit(1).ssdArray;
%%
% Get inhibition data from the session:
dataInh = cmd_inhibition(subjectID, sessionID, 'plotFlag', false);
% For now, use the grand SSRT via integratin method
if isempty(options.ssrt)
ssrt = dataInh.ssrt.grand;
else
    ssrt = options.ssrt;
end


%%
% Set up the plot with appropriate number of axes
if plotFlag
    
    epochName       = 'targOn';
    eventMarkName   = 'responseOnset';
    cMap = cmd_colormap;
    targLineW = 2;
    tickWidth = 10;
    nRow = 2;
    for kUnitIndex = 1 : nUnit
        
        
        for iTarg = 1 : nAngle
            iTarget = targAngleArray(iTarg);
            if nAngle > 1
                usableStopStop = dataInh.nStopStopAngle(iTarg,:) >= minTrialPerCond;
                usableStopTarg = dataInh.nStopTargAngle(iTarg,:) >= minTrialPerCond;
            else
                usableStopStop = dataInh.nStopStop >= minTrialPerCond;
                usableStopTarg = dataInh.nStopTarg >= minTrialPerCond;
            end
            usableStop = usableStopStop | usableStopTarg;
            if ~any(usableStop)
                continue
            end
            usableStopTarg = usableStopTarg(usableStop);
            usableStopStop = usableStopStop(usableStop);
            ssdArrayTarg = ssdArray(usableStop);
            
            % Find a suitable y-axis scale based on a few
            % conditions/outcomes
            dataMx    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'goTarg', 'stopTarg'}, targAngleArray(iTarg), ssdArray);
            sdfAll   = nanmean(spike_density_function(dataMx.raster(:,:), Kernel), 1);
            sdfMax = max(sdfAll .* 1.7);
            epochRange = ccm_epoch_range(epochName, 'plot');
            epochRange = 0 : 400;
            
            
            figureHandle = figureHandle + 1;
            nColumn = max(length(ssdArrayTarg), 4);
            if printPlot
                [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
            else
                [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
            end
            axisHeight = axisHeight * .9;
            clf
            
            
            
            
            
            
            
            
            
            
            
            for k = 1 : length(ssdArrayTarg)
                kSSD = ssdArrayTarg(k);
                % _______  Set up axes  ___________
                % axes names
                rStopTarg = 1;
                rStopStop = 2;
                
                % Get the go trial data: these need to be split to latency-match with
                % the stop trial data
                jGoTarg      = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'goTarg'}, iTarget);
                
                
                
                
                if usableStopTarg(k)
                    % Set up plot axes
                    ax(rStopTarg, k) = axes('units', 'centimeters', 'position', [xAxesPosition(rStopTarg, k) yAxesPosition(rStopTarg, k) axisWidth axisHeight]);
                    set(ax(rStopTarg, k), 'ylim', [0 sdfMax * 1.1], 'xlim', [epochRange(1) epochRange(end)])
                    cla
                    hold(ax(rStopTarg, k), 'on')
                    plot(ax(rStopTarg, k), [0 0], [0 sdfMax], '-k', 'linewidth', 2)
                    ttl = sprintf('SSD:  %d', kSSD);
                    title(ttl)
                    
                    if k > 1
                        set(ax(rStopTarg, k), 'yticklabel', [])
                        set(ax(rStopTarg, k), 'ycolor', [1 1 1])
                    end
                    
                    
                    % Get the stop trial data
                    jStopTarg    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'stopTarg'}, iTarget, kSSD);
                    
                    switch latencyMatchMethod
                        case 'ssrt'
                            jStopLatency = kSSD + ssrt;
                            jGoFastTrial = jGoTarg.eventLatency <= jStopLatency;
                        case 'mean'
                            jGoTargRT = sort(jGoTarg.eventLatency);
                            meanStopTargRT = nanmean(jStopTarg.eventLatency);
                            while nanmean(jGoTargRT) > meanStopTargRT
                                jGoTargRT(end) = [];
                            end
                            jStopLatency = jGoTargRT(end);
                            jGoFastTrial = jGoTarg.eventLatency <= jStopLatency;
                        case 'match'
                            % Use nearest neighbor method to get rt-matched
                            % trials
                            nStopCorrect = size(jStopStop.raster, 1);
                            data = ccm_match_rt(jGoTarg.eventLatency, jStopTarg.eventLatency, nStopCorrect);
                            jGoFastTrial = data.goFastTrial;
                    end
                    
                    jGoFastSDF   = nanmean(spike_density_function(jGoTarg.raster(jGoFastTrial,:), Kernel), 1);
                    
                    plot(ax(rStopTarg, k), [kSSD, kSSD], [0 sdfMax * .15], '--', 'color', [.2 .2 .2], 'linewidth', 1)
                    plot(ax(rStopTarg, k), [kSSD + ssrt, kSSD + ssrt], [0 sdfMax], 'color', [0 0 0], 'linewidth', 1)
                    plot(ax(rStopTarg, k), epochRange, jGoFastSDF(jGoTarg.align + epochRange), 'color', cMap.goTarg, 'linewidth', targLineW)
                    plot(ax(rStopTarg, k), epochRange, jStopTarg.sdf(jStopTarg.align + epochRange), 'color', cMap.stopTarg, 'linewidth', targLineW)
                    
                end
                
                
                
                
                
                if usableStopStop(k)
                    ax(rStopStop, k) = axes('units', 'centimeters', 'position', [xAxesPosition(rStopStop, k) yAxesPosition(rStopStop, k) axisWidth axisHeight]);
                    set(ax(rStopStop, k), 'ylim', [0 sdfMax * 1.1], 'xlim', [epochRange(1) epochRange(end)])
                    cla
                    hold(ax(rStopStop, k), 'on')
                    plot(ax(rStopStop, k), [0 0], [0 sdfMax], '-k', 'linewidth', 2)
                    if ~usableStopTarg(k)
                        ttl = sprintf('SSD:  %d', kSSD);
                        title(ttl)
                    end
                    
                    if k > 1
                        set(ax(rStopStop, k), 'yticklabel', [])
                        set(ax(rStopStop, k), 'ycolor', [1 1 1])
                    end
                    
                    % Get the stop trial data
                    jStopStop    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'stopStop'}, iTarget, kSSD);
                    jStopTarg    = cmd_concat_neural_conditions(Unit(kUnitIndex), epochName, eventMarkName, {'stopTarg'}, iTarget, kSSD);
                    
                    if isempty(jStopTarg.sdf) && ~strcmp(latencyMatchMethod, 'ssrt')
                        continue
                    end
                    switch latencyMatchMethod
                        case 'ssrt'
                            jStopLatency = kSSD + ssrt;
                            jGoSlowTrial = jGoTarg.eventLatency > jStopLatency;
                        case 'mean'
                            if isempty(jStopTarg.sdf)
                                continue
                            end
                            jGoTargRT = sort(jGoTarg.eventLatency);
                            meanStopTargRT = nanmean(jStopTarg.eventLatency);
                            while nanmean(jGoTargRT) > meanStopTargRT
                                jGoTargRT(end) = [];
                            end
                            jStopLatency = jGoTargRT(end);
                            jGoSlowTrial = jGoTarg.eventLatency > jStopLatency;
                        case 'match'
                             if isempty(jStopTarg.sdf)
                                continue
                            end
                           % Use nearest neighbor method to get rt-matched
                            % trials
                            nStopCorrect = size(jStopStop.raster, 1);
                            data = ccm_match_rt(jGoTarg.eventLatency, jStopTarg.eventLatency, nStopCorrect);
                            jGoSlowTrial = data.goSlowTrial;
                    end
                    
                    jGoSlowSDF   = nanmean(spike_density_function(jGoTarg.raster(jGoSlowTrial,:), Kernel), 1);
                    
                    plot(ax(rStopStop, k), [kSSD, kSSD], [0 sdfMax * .15], '--', 'color', [.2 .2 .2], 'linewidth', 1)
                    plot(ax(rStopStop, k), [kSSD + ssrt, kSSD + ssrt], [0 sdfMax], 'color', [0 0 0], 'linewidth', 1)
                    plot(ax(rStopStop, k), epochRange, jGoSlowSDF(jGoTarg.align + epochRange), 'color', cMap.goTarg, 'linewidth', targLineW)
                    plot(ax(rStopStop, k), epochRange, jStopStop.sdf(jStopStop.align + epochRange), 'color', cMap.stopStop, 'linewidth', targLineW)
                end
                
                
                
                
                
                
                
                
            end
            
            h=axes('Position', [0 0 1 1], 'Visible', 'Off');
            titleString = sprintf('%s \tUnit %s:\tAngle: %d', sessionID, Unit(kUnitIndex).name, targAngleArray(iTarg));
            text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top', 'interpreter', 'none')
            
            if printPlot
                localFigurePath = local_figure_path;
                print(figureHandle,[localFigurePath, sessionID, '_',Unit(kUnitIndex).name, '_', iTarg, '_ccm_neuron_stop_vs_go.pdf'],'-dpdf', '-r300')
            end
        end % iTarg
    end % for kUnitIndex = 1 : nUnit
end % if plotFlag



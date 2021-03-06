function Data = ccm_neuron_stop_vs_go(subjectID, sessionID, unitArray, options)

% ccm_neuron_stop_vs_go(subjectID, sessionID, options)
%
% Compares noncanceled stops trials vs. latency matched (fast) go trials and canceled stop trials vs. latency matched (slower) go trials.
%
% CANCEL TIMES ARE RELATIVE TO CHECKERBOARD ONSET, SO YOU NEED
% TO SUBTRACT SSD AND SSRT TO CALCULATE CLASSIC "CANCEL TIME" (RELATIVE TO
% SSRT)
%
%
% If called without any arguments, returns a default options structure.
% If options are input but one is not specified, it assumes default.
%
% input:
%   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%     unitArray       	= which units to analyze:
%           'each', <list of units in a cell array: e.g.:{'spikeUnit17a', spikeUnit17b'}>
%
% Possible options are (default listed first):
%     options.dataType        	= which recorded signal to analyze:
%           'neuron','lfp',eeg'
%     options.collapseSignal    = Collapse across signal strength (difficulty conditions)?
%            false, true
%     options.collapseTarg 	= collapse angle/directions of the CORRECT
%           TARGET within each hemifield
%           false, true
%     options.latencyMatchMethod  = wich method to use to match go trial latencies with canceled and noncanceled stop trials:
%           'ssrt','match','mean';
%     options.minTrialPerCond  	= how many trials must a condition have to
%           include in the analyses?
%     options.cellType  	= Are we treating it as a movement or fixaiton
%     cell? A movement cell with have higher firiring rate for go vs. stop,
%     and a fixaiton with be reversed
%           'move','fix'.
%
%     options.plotFlag       = true, false;
%     options.printPlot       = false, true;
%     options.figureHandle  = optional way to assign the figure to a handle

%%

% Set default options or return a default options structure
if nargin < 4
    %Data type to collect/analyze
    options.dataType        	= 'neuron';
    options.multiUnit        	= false;
    options.ANALYZE_NONCANCELED = true;
    options.ANALYZE_CANCELED    = true;
    %     options.ANALYZE_CANCELED    = [];
    
    options.normalize   	= false;
    options.collapseSignal   	= false;
    options.collapseTarg        = true;
    options.latencyMatchMethod 	= 'ssrt';
    options.minTrialPerCond     = 10;
    options.cellType            =      'move';
    options.ssrt            =      'intWeightPerSession'; % 'intWeightPerSession, inttPerSession, intWeightPerCoherence, intPerCoherence, intPerSsd
    options.USE_PRE_SSD = true;
    options.ms2Std = 75;  % How many consecutive ms should neural signals be separated to "cancel"?
    
    options.plotFlag            = true;
    options.plotSem            = true;
    options.printPlot           = true;
    options.plotSingle          = false;
    
    options.category = [];
    
    % Return just the default options struct if no input
    if nargin == 0
        Data           = options;
        return
    end
end

if ~isempty(options.category)
    appendNorm = [];
    if options.normalize
        appendNorm = '_norm';
        options.category = [options.category, appendNorm];
    end
end

% subjectID = 'broca';
% sessionID = 'bp093n02';
% unitArray = {'spikeUnit17a'};

ANALYZE_CANCELED = options.ANALYZE_CANCELED;
ANALYZE_NONCANCELED = options.ANALYZE_NONCANCELED;

maxSdfProportion = .67;
alphaVal = 0.05;  % p-value significance threshold.


dataType            = options.dataType;
collapseSignal      = options.collapseSignal;
latencyMatchMethod  = options.latencyMatchMethod;
minTrialPerCond     = options.minTrialPerCond;

plotFlag            = options.plotFlag;
printPlot           = options.printPlot;
stopStopFig        = 100;
stopTargFig        = 150;
nPlot               = 6; % plot only 3 conditions- those with highest number of stop trials
faceAlpha           = .25; % how opaque to make the error ribbons.

% latencyMatchMethod = 'match';
% latencyMatchMethod = 'rt';


Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;


spikeWindow = -20 : 20;


epochName       = 'checkerOn';
markEvent   = 'responseOnset';
encodeTime      = 10;

epochRangeChecker      = 1 : 600;  % Make this big enough to catch really late cancel times.
epochRangeSacc      = -199 : 100;

normBaseEpoch = -100:-1; % How many ms prior target onset to use to calculate mean and Std Sdfs for normalization;

%%  Loop through Units and target pairs to collect and plot data
[nUnit, nTargPair]  = size(unitArray);
for kUnitIndex = 1 : nUnit
    
    fprintf('%s: %s\n',sessionID, (unitArray{kUnitIndex}))
    
    
    %  Get neural data from the session/unit:
    tic
    optSess             = ccm_options;
    optSess.multiUnit    = options.multiUnit;
    optSess.plotFlag    = false;
    optSess.printPlot    = false;
    optSess.collapseTarg = options.collapseTarg;
    optSess.unitArray   = unitArray(kUnitIndex);
    optSess.USE_PRE_SSD   = options.USE_PRE_SSD;
    Unit                = ccm_session_data(subjectID, sessionID, optSess);
    dataInh = Unit(1).dataInh;
    rf = ccm_find_saccade_rf(Unit);
    
    fprintf('ccm_session_data time: %.2f\n', toc)
    if isempty(Unit)
        fprintf('Session %s does not contain spike data \n', sessionID)
        return
    end
    
    
    
    pSignalArray        = Unit(1).pSignalArray;
    ssdArray            = Unit(1).ssdArray;
    
    nSSD             = length(ssdArray);
    
    
    
%     %   Get the saccade receptive field
%     %   Get neural data from the session/unit:
%     optCollapse             = ccm_options;
%     optCollapse.multiUnit    = options.multiUnit;
%     optCollapse.plotFlag    = false;
%     optCollapse.printFlag    = false;
%     optCollapse.collapseTarg = true;
%     optCollapse.collapseSignal = true;
%     optCollapse.unitArray   = unitArray(kUnitIndex);
%     optCollapse.doStops   = false;
%     UnitCollapse                = ccm_session_data(subjectID, sessionID, optCollapse);
%     rf = ccm_find_saccade_rf(UnitCollapse);
    
    
    
    %   Use rf to determine which color coherence conditions to analyze
    switch rf
        case 'left'
            colorInd = pSignalArray < .5;
        case 'right'
            colorInd = pSignalArray > .5;
        case 'none'
            switch Unit(1).hemisphere
                case 'left'
                    colorInd = pSignalArray > .5;
                case 'right'
                    colorInd = pSignalArray < .5;
            end
    end
    rfColorCohArray = pSignalArray(colorInd);
    if rfColorCohArray(1) < .5
        cohEasy = rfColorCohArray(1);
        cohHard = rfColorCohArray(end);
    else
        cohEasy = rfColorCohArray(end);
        cohHard = rfColorCohArray(1);
    end
    
    
    
    %%
    
    % Use only the color coherence conditions for the response field figured
    % out above
    nStopStop = dataInh.nStopStop(colorInd,:);
    [nStopStop, indSS] = sort(nStopStop(:));
    nStopTarg = dataInh.nStopTarg(colorInd,:);
    [nStopTarg, indST] = sort(nStopTarg(:));
    
    % Create arrays of the color coherence and ssd values that are
    % usable (have enough data). Will loop through them below for each
    % set of comparisons
    usableStopStop = nStopStop >= minTrialPerCond;
    usableStopTarg = nStopTarg >= minTrialPerCond;
    
    % Vectorize color coherence and ssd matrices
    colorCohMat = repmat(rfColorCohArray(:),1,nSSD);
    colorCohMat = colorCohMat(:);
    
    ssdMat = repmat(reshape(ssdArray,1,nSSD),length(rfColorCohArray),1);
    ssdMat = ssdMat(:);
    
    % Use sorted trial numbers to make lists of color coherence
    % and ssd values to loop through for each condition
    stopStopCoh = colorCohMat(indSS);
    stopStopCoh = flipud(stopStopCoh(usableStopStop));
    stopStopSsd = ssdMat(indSS);
    stopStopSsd = flipud(stopStopSsd(usableStopStop));
    stopStopSsrt = nan(length(stopStopSsd), 1);
    
    stopTargCoh = colorCohMat(indST);
    stopTargCoh = flipud(stopTargCoh(usableStopTarg));
    stopTargSsd = ssdMat(indST);
    stopTargSsd = flipud(stopTargSsd(usableStopTarg));
    stopTargSsrt = nan(length(stopTargSsd), 1);
    
    % Code the easy and hard color coherences for later analyses
    stopStopCond = nan(length(stopStopCoh), 1);
    stopStopCond(stopStopCoh == cohEasy) = 1;
    stopStopCond(stopStopCoh == cohHard) = 2;
    
    stopTargCond = nan(length(stopTargCoh), 1);
    stopTargCond(stopTargCoh == cohEasy) = 1;
    stopTargCond(stopTargCoh == cohHard) = 2;
    
    
    
    
    
    
    
    
    % Get all GO RTs and all STOP RTs to use as baseline differential SDF for
    % standard deviations to determine neural separation time
    opt                 = options; % Get default options structure
    
    opt.epochName       = 'checkerOn';
    opt.eventMarkName   = 'targOn';
    opt.colorCohArray   = rfColorCohArray;
    opt.conditionArray  = {'goTarg', 'goDist'};
    opt.ssdArray        = [];
    goChecker           = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
    
    opt.conditionArray  = {'stopStop'};
    opt.ssdArray        = ssdArray;
    stopChecker         = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
    
    
    if options.normalize
        opt.epochName       = 'targOn';
        opt.eventMarkName   = 'fixWindowEntered';
        opt.colorCohArray   = pSignalArray;
        opt.conditionArray  = {'goTarg', 'goDist'};
        opt.ssdArray        = [];
        goTargOnNorm           = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
        
        goTargOnNormSdf = nanmean(spike_density_function(goTargOnNorm.signal(:, :), Kernel), 1);
        goTargOnSdfMean    = nanmean(goTargOnNormSdf(goTargOnNorm.align + normBaseEpoch));
        goTargOnSdfSD      = nanstd(goTargOnNormSdf(goTargOnNorm.align + normBaseEpoch));
        
    end
    
    
    % Pre-allocate cell arrays
    cancelTimeDist = cell(length(stopStopCoh), 1);
    cancelTimeSdf = nan(length(stopStopCoh), 1);
    cancelTime6Std = nan(length(stopStopCoh), 1);
    cancelTime4Std = nan(length(stopStopCoh), 1);
    cancelTime2Std = nan(length(stopStopCoh), 1);  % The index in the sdf after checker onset at wich the differential sdf passes the 50ms test
    pValue40msStopStop = nan(length(stopStopCoh), 1);  % The index in the sdf after checker onset at wich the differential sdf passes the 50ms test
    
    pValue40msStopTarg = nan(length(stopTargCoh), 1);  % The index in the sdf after checker onset at wich the differential sdf passes the 50ms test
    
    
    pValuePostSsrtRt = nan(length(stopStopCoh), 1);  % for jeff's new proposed analysis
    pValuePostSsrt100 = nan(length(stopStopCoh), 1);  % for jeff's new proposed analysis
    targToStopDiffSD = nan(length(stopStopCoh), 1);  % for jeff's new proposed analysis
    preSaccDiffSD = nan(length(stopStopCoh), 1);  % for jeff's new proposed analysis
    
    
    
    % Canceled Stop and Slow Go
    stopStopCheckerData          = cell(length(stopStopCoh), 1);
    stopStopCheckerFn          = cell(length(stopStopCoh), 1);
    stopStopCheckerFnSem          = cell(length(stopStopCoh), 1);
    stopStopCheckerEventLat  	= cell(length(stopStopCoh), 1);
    stopStopCheckerAlign      	= nan(length(stopStopCoh), 1);
    stopStopCheckerSdf      	= cell(length(stopStopCoh), 1);
    
    goTargSlowCheckerData        = cell(length(stopStopCoh), 1);
    goTargSlowCheckerFn        = cell(length(stopStopCoh), 1);
    goTargSlowCheckerFnSem        = cell(length(stopStopCoh), 1);
    goTargSlowCheckerEventLat  	= cell(length(stopStopCoh), 1);
    goTargSlowCheckerAlign    	= nan(length(stopStopCoh), 1);
    goTargSlowCheckerSdf        = cell(length(stopStopCoh), 1);
    
    goTargSlowSaccData        	= cell(length(stopStopCoh), 1);
    goTargSlowSaccFn           = cell(length(stopStopCoh), 1);
    goTargSlowSaccFnSem           = cell(length(stopStopCoh), 1);
    goTargSlowSaccEventLat      = cell(length(stopStopCoh), 1);
    goTargSlowSaccAlign         = nan(length(stopStopCoh), 1);
    goTargSlowSaccSdf           = cell(length(stopStopCoh), 1);
    
    
    stopStopSpike               = cell(length(stopStopCoh), 1);
    goTargSlowSpike               = cell(length(stopStopCoh), 1);
    stopStopSpikePostSsrtRt               = cell(length(stopStopCoh), 1);
    goTargSlowSpikePostSsrtRt               = cell(length(stopStopCoh), 1);
    stopStopSpikePostSsrt100               = cell(length(stopStopCoh), 1);
    goTargSlowSpikePostSsrt100               = cell(length(stopStopCoh), 1);
    
    
    % Nonanceled Stop and Fast Go
    stopTargCheckerData          = cell(length(stopTargCoh), 1);
    stopTargCheckerFn          = cell(length(stopTargCoh), 1);
    stopTargCheckerFnSem          = cell(length(stopTargCoh), 1);
    stopTargCheckerEventLat  	= cell(length(stopTargCoh), 1);
    stopTargCheckerAlign      	= nan(length(stopTargCoh), 1);
    stopTargCheckerSdf          = cell(length(stopTargCoh), 1);
    
    stopTargSaccData             = cell(length(stopTargCoh), 1);
    stopTargSaccFn             = cell(length(stopTargCoh), 1);
    stopTargSaccFnSem             = cell(length(stopTargCoh), 1);
    stopTargSaccEventLat        = cell(length(stopTargCoh), 1);
    stopTargSaccAlign           = nan(length(stopTargCoh), 1);
    stopTargSaccSdf             = cell(length(stopTargCoh), 1);
    
    goTargFastCheckerData        = cell(length(stopTargCoh), 1);
    goTargFastCheckerFn        = cell(length(stopTargCoh), 1);
    goTargFastCheckerFnSem        = cell(length(stopTargCoh), 1);
    goTargFastCheckerEventLat  	= cell(length(stopTargCoh), 1);
    goTargFastCheckerAlign    	= nan(length(stopTargCoh), 1);
    goTargFastCheckerSdf        = cell(length(stopTargCoh), 1);
    
    goTargFastSaccData        	= cell(length(stopTargCoh), 1);
    goTargFastSaccFn           = cell(length(stopTargCoh), 1);
    goTargFastSaccFnSem           = cell(length(stopTargCoh), 1);
    goTargFastSaccEventLat      = cell(length(stopTargCoh), 1);
    goTargFastSaccAlign         = nan(length(stopTargCoh), 1);
    goTargFastSaccSdf           = cell(length(stopTargCoh), 1);
    
    stopTargSpike               = cell(length(stopTargCoh), 1);
    goTargFastSpike               = cell(length(stopTargCoh), 1);
    
    
    
    
    
    
    
    % Loop through each (set of) usable (enough trials as per
    % minTrialPerCond) signal strength/ssd condition.
    
    
    
    
    % Creat a new figure if needed (for a new pair of targets and/or
    % left/right checker signal
    if plotFlag
        grayScale = 1;
        cMap = ccm_colormap(pSignalArray, grayScale);
        if length(pSignalArray) > 4
            stopColorMap = [1 0 0; .7 0 0; .4 0 0; .4 0 0; .7 0 0; 1 0 0];
        elseif length(pSignalArray) == 4
            stopColorMap = [1 0 0; .7 0 0; .7 0 0; 1 0 0];
        end
        
        opt = ccm_concat_neural_conditions; % Get default options structure
        opt.colorCohArray = rfColorCohArray;
        opt.ssdArray = ssdArray;
        
        opt.epochName = 'targOn';
        opt.markEvent = markEvent;
        opt.conditionArray = {'goTarg'};
        dataTarg           = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
        
        opt.epochName = 'responseOnset';
        opt.markEvent = markEvent;
        opt.conditionArray = {'goTarg'};
        dataSacc           = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
        
        
        
        if options.normalize
            range = -500:300;
            sdfSacc = nanmean(spike_density_function(dataSacc.signal(:, :), Kernel), 1);
            sdfSacc = (sdfSacc(dataSacc.align+range) - goTargOnSdfMean) / goTargOnSdfSD;
            max(sdfSacc)
            yLimMax         = max(sdfSacc) + .5 * max(sdfSacc)
            yLimMax         = max(sdfSacc) * 1.2
            yLimMax         = max(sdfSacc) + 6 * goTargOnSdfSD
            %             yLimMin         = min(sdfSacc) - .2 * min(sdfSacc);
            yLimMin         = -20;
        else
            yLimMax         = max(nanmean(spike_density_function(dataSacc.signal(:,:), Kernel), 1)) * 1.2;
            yLimMin = 0;
            
        end
        
        
        
        goLineW = 2;
        stopLineW = 2;
        markSize = 30;
        tickWidth = 10;
        nRow = nPlot;
        nColumn = 2;
        colChkr = 1;
        colSacc = 2;
        
        stopStopFig = stopStopFig + 1;
        stopTargFig = stopTargFig + 1;
        %         if printPlot
        [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nRow, nColumn, 'portrait', stopStopFig);
        clf
        h=axes('Position', [0 0 1 1], 'Visible', 'Off');
        set(gcf, 'Name','Go v Canceled','NumberTitle','off')
        titleString = sprintf('%s \t %s', sessionID, Unit(kUnitIndex).name);
        text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')
        
        if options.plotSingle
            stopStopFigSingle = 2039;
            [axisW, axisH, xAxesPos, yAxesPos] = standard_figure(2      , 2, 'portrait', stopStopFigSingle);
        end
        
        if ANALYZE_NONCANCELED
            [~, ~, ~, ~] = standard_figure(nRow, nColumn, 'portrait', stopTargFig);
            clf
            h=axes('Position', [0 0 1 1], 'Visible', 'Off');
            set(gcf, 'Name','Go v Noncanceled','NumberTitle','off')
            titleString = sprintf('%s \t %s', sessionID, Unit(kUnitIndex).name);
            text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')
            axisHeight = axisHeight * .9;
        end
    end
    
    
    
    
    
    
    
    %             $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    %             $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    %             $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    %%
    if ANALYZE_CANCELED
        %               Canceled Stop vs latency-matched (Slow) Go:
        if isempty(stopStopCoh)
            fprintf('Session %s doesn"t have enough canceled trials in any condiiton to analyze\n', sessionID)
        else
            fprintf('\nGo vs. Canceled:\n')
            for i = 1 : length(stopStopCoh)
                
                switch options.ssrt
                    case 'intWeightPerSession'
                        iSsrt = round(nanmean(dataInh.ssrtIntegrationWeighted));
                    case 'intPerSsd'
                        % Find the ssrt to use
                        iCohSsrtInd = dataInh.pSignalArray == stopStopCoh(i);
                        iSsdSsrtInd = dataInh.ssd{iCohSsrtInd} == stopStopSsd(i);
                        iSsrt = dataInh.ssrtIntegration{iCohSsrtInd}(iSsdSsrtInd);
                end
                if isnan(iSsrt)
                    continue
                end
                stopStopSsrt(i) = iSsrt;
                
                
                % Get the go trial data: these need to be split to latency-match with
                % the stop trial data
                opt                 = options; % Get default options structure
                
                opt.epochName       = epochName;
                opt.eventMarkName   = markEvent;
                opt.conditionArray  = {'goTarg'};
                opt.colorCohArray   = stopStopCoh(i);
                opt.ssdArray        = [];
                iGoTargChecker      = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
                
                opt.epochName       = 'responseOnset';
                opt.eventMarkName   = 'checkerOn';
                iGoTargSacc         = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
                
                opt.epochName       = 'targOn';
                opt.eventMarkName   = 'fixWindowEntered';
                iGoTargTargOn         = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
                
                
                
                % Get the stop trial data
                opt.epochName       = epochName;
                opt.eventMarkName   = 'targOn';
                opt.conditionArray  = {'stopStop'};
                opt.ssdArray        = stopStopSsd(i);
                iStopStopChecker    = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
                
                opt.epochName       = epochName;
                opt.eventMarkName   = markEvent;
                opt.conditionArray  = {'stopTarg'};
                iStopTargChecker       = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
                
                
                % Figure out the minimum length and SDF/raster can be to
                % analyze:
                checkerEpochEnd = min([epochRangeChecker(end), length(iStopStopChecker.signalFn) - iStopStopChecker.align - stopStopSsd(i)]);
                
                % Use a latency-matching method to get a (slow)
                % subsample of go trials that match the canceled
                % stop trials
                switch latencyMatchMethod
                    
                    case 'ssrt'
                        iStopLatency = stopStopSsd(i) + iSsrt;
                        iGoSlowTrial = iGoTargChecker.eventLatency > iStopLatency;
                    case 'mean'
                        iGoTargRT = sort(iGoTargChecker.eventLatency);
                        meanStopTargRT = nanmean(iStopTargChecker.eventLatency);
                        while nanmean(iGoTargRT) > meanStopTargRT
                            iGoTargRT(end) = [];
                        end
                        iStopLatency = iGoTargRT(end);
                        iGoSlowTrial = iGoTargChecker.eventLatency > iStopLatency;
                    case 'match'
                        % Use nearest neighbor method to get rt-matched
                        % trials
                        nStopCorrect = size(iStopStopChecker.signal, 1);
                        data = ccm_match_rt(iGoTargChecker.eventLatency, iStopTargChecker.eventLatency(iStopTargTrial), nStopCorrect);
                        iGoSlowTrial = data.goSlowTrial;
                end
                
                iGoSlowRt = iGoTargChecker.eventLatency(iGoSlowTrial);
                
                
                
                
                
                % ****************************************************************
                % COLLECT (AND ANALYZE) THE RELEVANT DATA
                
                stopStopCheckerData{i}        = iStopStopChecker.signal;
                stopStopCheckerAlign(i)     	= iStopStopChecker.align;
                stopStopCheckerEventLat{i}	= iStopStopChecker.eventLatency;
                stopStopCheckerSdf{i} = spike_density_function(iStopStopChecker.signal, Kernel);
                
                goTargSlowCheckerData{i}      = iGoTargChecker.signal(iGoSlowTrial,:);
                goTargSlowCheckerAlign(i) 	= iGoTargChecker.align;
                goTargSlowCheckerEventLat{i}	= iGoTargChecker.eventLatency(iGoSlowTrial,:);
                goTargSlowCheckerSdf{i} = spike_density_function(goTargSlowCheckerData{i}, Kernel);
                
                goTargSlowSaccData{i}         = iGoTargSacc.signal(iGoSlowTrial,:);
                goTargSlowSaccAlign(i)       = iGoTargSacc.align;
                goTargSlowSaccEventLat{i}    = iGoTargSacc.eventLatency(iGoSlowTrial,:);
                goTargSlowSaccSdf{i} = spike_density_function(goTargSlowSaccData{i}, Kernel);
                
                
                % If we want to normalize (z-score-ish) SDFs,
                if options.normalize
                    stopStopCheckerSdf{i} = (stopStopCheckerSdf{i} - goTargOnSdfMean) / goTargOnSdfSD;
                    goTargSlowCheckerSdf{i} = (goTargSlowCheckerSdf{i} - goTargOnSdfMean) / goTargOnSdfSD;
                    goTargSlowSaccSdf{i} = (goTargSlowSaccSdf{i} - goTargOnSdfMean) / goTargOnSdfSD;
                end
                
                
                stopStopCheckerFn{i} 	= nanmean(stopStopCheckerSdf{i});
                %                 stopStopCheckerFn{i} 	= iStopStopChecker.signalFn;
                goTargSlowCheckerFn{i}      = nanmean(goTargSlowCheckerSdf{i});
                goTargSlowSaccFn{i}         = nanmean(goTargSlowSaccSdf{i});
                
                stopStopCheckerFnSem{i} 	= std(stopStopCheckerSdf{i}) / sqrt(size(stopStopCheckerSdf{i}, 1));
                goTargSlowCheckerFnSem{i}      = std(goTargSlowCheckerSdf{i}) / sqrt(size(goTargSlowCheckerSdf{i}, 1));
                goTargSlowSaccFnSem{i}         = std(goTargSlowSaccSdf{i}) / sqrt(size(goTargSlowSaccSdf{i}, 1));
                
                
                
                
                % TEST #1:
                % Hanes et al 1998 (p.822) t-test of spike rates 40 ms surrounding estimated ssrt
                stopStopSpike{i}     = sum(iStopStopChecker.signal(:, spikeWindow + iStopStopChecker.align + stopStopSsd(i) + iSsrt), 2) / length(spikeWindow) * 1000;
                goTargSlowSpike{i}   = sum(iGoTargChecker.signal(iGoSlowTrial, spikeWindow + iGoTargChecker.align + stopStopSsd(i) + iSsrt), 2) / length(spikeWindow) * 1000;
                [h,p,ci,sts]       	= ttest2(stopStopSpike{i}, goTargSlowSpike{i});
                
                pValue40msStopStop(i)    = p;
                stats{i}     = sts;
                
                
                % TEST #1 ALTERED BY J SCHALL
                % t-test of spike rates from ssrt to median RT
                stopStopSpikePostSsrtRt{i}     = sum(iStopStopChecker.signal(:, iStopStopChecker.align + stopStopSsd(i) + iSsrt : iStopStopChecker.align + nanmedian(iGoSlowRt)), 2) / (nanmedian(iGoSlowRt)-(iStopStopChecker.align + stopStopSsd(i) + iSsrt)) * 1000;
                goTargSlowSpikePostSsrtRt{i}   = sum(iGoTargChecker.signal(iGoSlowTrial, iGoTargChecker.align + stopStopSsd(i) + iSsrt : iGoTargChecker.align + nanmedian(iGoSlowRt)), 2) / (nanmedian(iGoSlowRt)-(iStopStopChecker.align + stopStopSsd(i) + iSsrt)) * 1000;
                [h,p,ci,~]       	= ttest2(stopStopSpikePostSsrtRt{i}, goTargSlowSpikePostSsrtRt{i});
                
                pValuePostSsrtRt(i)    = p;
                
                % t-test of spike rates from ssrt to ssrt + 100
                stopStopSpikePostSsrt100{i}     = sum(iStopStopChecker.signal(:, iStopStopChecker.align + stopStopSsd(i) + iSsrt : iStopStopChecker.align + stopStopSsd(i) + iSsrt + 100), 2) / (nanmedian(iGoSlowRt)-(iStopStopChecker.align + stopStopSsd(i) + iSsrt)) * 1000;
                goTargSlowSpikePostSsrt100{i}   = sum(iGoTargChecker.signal(iGoSlowTrial, iGoTargChecker.align + stopStopSsd(i) + iSsrt : iGoTargChecker.align + stopStopSsd(i) + iSsrt + 100), 2) / (nanmedian(iGoSlowRt)-(iStopStopChecker.align + stopStopSsd(i) + iSsrt)) * 1000;
                [h,p,ci,~]       	= ttest2(stopStopSpikePostSsrt100{i}, goTargSlowSpikePostSsrt100{i});
                
                pValuePostSsrt100(i)    = p;
                
                
                % TEST #2:
                % Hanes et al 1998 (p.822) differential sdf test
                % ------------------------------------------------------------------------
                sdfDiffCheckerOn = 500;
                % Differential sdf go - stopStop
                
                % Baseline sdf differential to caclulate SDF
                sdfDiffBase = goChecker.signalFn(goChecker.align + (-sdfDiffCheckerOn : 1))' - stopChecker.signalFn(stopChecker.align + (-sdfDiffCheckerOn : 1))';
                
                % Get the standard deviation of the differential sdf during the
                % 500 ms before checkerboard onset
                stdDiff = std(sdfDiffBase);
                
                % GoTarg - stopStop SDF from checker Onset to end of checker epoch.
                sdfDiff = goTargSlowCheckerFn{i}(iGoTargChecker.align+stopStopSsd(i) : iGoTargChecker.align+stopStopSsd(i) + checkerEpochEnd)' - stopStopCheckerFn{i}(iStopStopChecker.align+stopStopSsd(i) : iStopStopChecker.align+stopStopSsd(i) + checkerEpochEnd)';
                
                % If user thinks it's a fixation/stopping type cell, flip
                % the sign of the differential sdf
                if strcmp(options.cellType, 'fix')
                    sdfDiff = -sdfDiff;
                end
                
                
                
                % are there times at which the difference between sdfs is
                % greater than 2 standard deviations of the difference 500
                % ms before checkerboard onset? Check from checkerboard onset
                % to end of the checkerboard epoch.
                % So std2Ind starts at checkerboard onset.
                std2Ind = sdfDiff > 2*stdDiff;
                
                % Look for a sequence of options.ms2Std ms for which the go sdf is 2
                % std greater than the stop sdf.
                % First whether the differential sdf was > 2*Std for the
                % first options.ms2Std ms
                if sum(std2Ind(1:options.ms2Std)) == options.ms2Std
                    cancelTime2Std(i) = -iSsrt;
                else
                    % If it wasn't, determein whether there was a time
                    % after the checkerboard onset that the differential
                    % sdf was > 2*Std for at least options.ms2Std ms.
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
                        
                        % If there's one more riseAbove2Std than sinkBelow2Std, the last riseAbove2Std
                        % will last until the end of the sdf: Add to
                        % singkBelowStd the end of the epoch
                        if length(riseAbove2Std) > length(sinkBelow2Std)
                            sinkBelow2Std = [sinkBelow2Std; checkerEpochEnd];
                        end
                        
                        % Now riseAbove2Std length should be equal. See if
                        % any of the riseAbove2Std streaks go longer than
                        % 50ms
                        ind = find(sinkBelow2Std - riseAbove2Std >= options.ms2Std, 1);
                        if ~isempty(ind)
                            cancelTime2Std(i) = riseAbove2Std(ind) - iSsrt;
                        end
                        % If they're not equal, the last riseAbove2Std
                        % will last until the end of the sdf: Add to
                        % singkBelowStd the end of the epoch
                    end
                end
                
                
                % If there werwe options.ms2Std consecutive ms of sdfDiff > 2*Std,
                % check whether the difference ever reached 4*Std
                if ~isnan(cancelTime2Std(i))
                    std4Ind = sdfDiff(cancelTime2Std(i)+stopStopSsd(i)+iSsrt : end) > 4*stdDiff;
                    if sum(std4Ind)
                        
                        % Cancel time will be Negative if
                        % go SDF > stop SDF BEFORE SSRT (which is a good thing-
                        % means the neuron can contribute to controlling saccade
                        % initiation/inhibition.
                        cancelTime4Std(i) = find(std4Ind, 1) + cancelTime2Std(i);
                    end
                    
                    % If it reached 4*Std, check whether and when it reached
                    % 6*Std:
                    if ~isnan(cancelTime4Std(i))
                        std6Ind = sdfDiff(cancelTime4Std(i)+stopStopSsd(i)+iSsrt : end) > 6*stdDiff;
                        if sum(std6Ind)
                            cancelTime6Std(i) = find(std6Ind, 1) + cancelTime4Std(i);
                        end
                    end
                end
                
                
                
                
                
                % TEST #2: ALTERED BY J SCHALL
                % 1. Compare SDFs from checkerOn to stopSignalOn
                % 2. Compare SDFs from SSRT to median goSlowRT (Canceled trials) - 20s pre-saccade (goSlowTrials)
                % ------------------------------------------------------------------------
                
                % 1.
                targToStopGoSlowMean = mean(goTargSlowCheckerFn{i}(iGoTargChecker.align : iGoTargChecker.align + stopStopSsd(i)));
                targToStopGoSlowStd = std(goTargSlowCheckerFn{i}(iGoTargChecker.align : iGoTargChecker.align + stopStopSsd(i)));
                targToStopStopStopMean = mean(stopChecker.signalFn(stopChecker.align : stopChecker.align + stopStopSsd(i)));
                targToStopDiffSD(i) = (targToStopStopStopMean - targToStopGoSlowMean) / targToStopGoSlowStd;
                
                % 2.
                preSaccGoSlowMean = mean(goTargSlowSaccFn{i}(iGoTargSacc.align-20 : iGoTargSacc.align));
                preSaccGoSlowStd = std(goTargSlowSaccFn{i}(iGoTargSacc.align-20 : iGoTargSacc.align));
                ssrtToRtStopStopMean = mean(stopChecker.signalFn(iStopStopChecker.align + stopStopSsd(i) + iSsrt : iStopStopChecker.align + nanmedian(iGoSlowRt)));
                preSaccDiffSD(i) = (ssrtToRtStopStopMean - preSaccGoSlowMean) / preSaccGoSlowStd;
                
                
                
                
                
                
                
                
                
                
                % Cancel time trial-by-trial distribution and average SDF deflection analyses: no point in latency matching, since latency-matching depends on mean activation function/SDF.
                % Slow GO is commented out for this reason: see
                % also plot_inhibition.m
                % #############################################################################
                % #############################################################################
                % #############################################################################
                
                % Find maximum of SDFs after SSD
                % May need to cut off search for max before
                % the end of the signal, so later modulations
                % don't interfere with the analysis.
                stopSignalInd = iStopStopChecker.align + stopStopSsd(i);
                % Use the mean of the latency-matched Go
                % trials as the last index to look for maximum
                % SDF value during canceled stop trials
                meanRtInd = ceil(nanmean(iGoTargChecker.eventLatency(iGoSlowTrial,:)));
                [maxStopStopSdf, maxInd] = max(stopStopCheckerSdf{i}(:, stopSignalInd : stopSignalInd + meanRtInd), [], 2);
                %                 maxStopStopSdf(maxInd >= meanRtInd) = nan;
                %                 maxInd(maxInd >= meanRtInd) = nan;
                
                
                % Average SDF deflection analysis
                % ------------------------------------------------------------------------
                [maxSdfDiff, maxSdfDiffInd] = max(sdfDiff(1 : meanRtInd-stopStopSsd(i)));
                %                 [maxStopStopFn, maxFnInd] = max(stopStopCheckerFn{i}(stopSignalInd : stopSignalInd + meanRtInd-stopStopSsd(i)-50))
                [maxStopStopFn, maxFnInd] = max(stopStopCheckerFn{i}(stopSignalInd : stopSignalInd + maxSdfDiffInd));
                
                
                % Cancel time only gets calculated if max SDF occurs before mean Go slow Rts
                iCancelTimeSdf = nan;
                if maxFnInd < meanRtInd
                    
                    % Find the minimum after the maximum, to
                    % determine relative half-max index
                    [minStopStopFn, minFnInd] = min(stopStopCheckerFn{i}(stopSignalInd + maxFnInd : stopSignalInd + meanRtInd-stopStopSsd(i)));
                    
                    % Use a proportion of the maximum relative to minimum to determine cancel time
                    halfMaxStopStopFn = minStopStopFn + (maxStopStopFn - minStopStopFn) * maxSdfProportion;
                    
                    halfMaxFnInd = find(stopStopCheckerFn{i}(stopSignalInd + maxFnInd:end) < halfMaxStopStopFn, 1);
                    
                    % To calculate cancel time, use half the time between max and "half
                    % max", relative to ssrt.
                    % Also, must pass the pValuePostSsrt100(i) test
                    
                    if ~isempty(halfMaxFnInd) && pValuePostSsrt100(i) < alphaVal
                        iCancelTimeSdf = maxFnInd + halfMaxFnInd/2 - iSsrt;
                        %                                                     % optional- plot individual trials for
                        %                                                     % troubleshooting
                        %                                                     figure(27)
                        %                                                     clf
                        %                                                     hold 'all'
                        %                                                     % Canceled GO goTargSlowCheckerFn
                        %                                                     plot(goTargSlowCheckerFn{i}(goTargSlowCheckerAlign(i):goTargSlowCheckerAlign(i)+700), 'g')
                        %                                                     plot(stopStopCheckerFn{i}(iStopStopChecker.align:iStopStopChecker.align+700), 'k')
                        %                                                     plot([stopStopSsd(i) stopStopSsd(i)], [0 maxStopStopFn*.8], '-r', 'linewidth', 3)
                        %                                                     plot([stopStopSsd(i)+iSsrt stopStopSsd(i)+iSsrt], [0 maxStopStopFn*.8], '--r', 'linewidth', 3)
                        %                                                     plot([stopStopSsd(i) + maxFnInd, stopStopSsd(i) + maxFnInd], [0 maxStopStopFn], '-k', 'linewidth', 3)
                        %                                                     plot([stopStopSsd(i) + maxFnInd + halfMaxFnInd, stopStopSsd(i) + maxFnInd + halfMaxFnInd], [0 halfMaxStopStopFn], '-k', 'linewidth', 3)
                        %                                                     plot([stopStopSsd(i)+(maxFnInd + halfMaxFnInd/2), stopStopSsd(i)+(maxFnInd + halfMaxFnInd/2)], [0 maxStopStopFn*.8], '--k', 'linewidth', 3)
                        %                         printName = sprintf('%s_canceled_Coh%s_Ssd%d_nTrial%d.pdf',Unit(kUnitIndex).name, num2str(stopStopCoh(i)*100), stopStopSsd(i), size(stopStopCheckerSdf{i}, 1));
                        %                                     print(gcf,fullfile(local_figure_path, subjectID, 'go_vs_canceled', options.ssrt, sessionID,  printName),'-dpdf', '-r300')
                    end
                end
                
                % ------------------------------------------------------------------------
                
                % Trial-by-trial analysis: For each trial, determine first index of SDF
                % that falls below half-max
                % ------------------------------------------------------------------------
                iCancelTime = nan(size(iStopStopChecker.signal, 1), 1);
                for a = 1 : size(iStopStopChecker.signal, 1)
                    
                    % Cancel time only gets calculated if max SDF occurs before mean Go slow Rts
                    if maxInd(a) < meanRtInd
                        
                        % Find the minimum after the maximum, to
                        % determine relative half-max index
                        [aMinStopStopSdf, minInd] = min(stopStopCheckerSdf{i}(a, stopSignalInd + maxInd(a) : stopSignalInd + meanRtInd), [], 2);
                        
                        % Use a proportion of the maximum relative to minimum to determine cancel time
                        aHalfMaxStopStopSdf = aMinStopStopSdf + (maxStopStopSdf(a) - aMinStopStopSdf) * maxSdfProportion;
                        
                        aHalfMaxInd = find(stopStopCheckerSdf{i}(a,stopSignalInd + maxInd(a):end) < aHalfMaxStopStopSdf, 1);
                        
                        if ~isempty(aHalfMaxInd)
                            %                             % optional- plot individual trials for
                            %                             % troubleshooting
                            %                             figure(24)
                            %                             clf
                            %                             hold 'all'
                            %                             % Canceled GO
                            %                             plot(stopStopCheckerSdf{i}(a,iStopStopChecker.align:iStopStopChecker.align+700), 'k')
                            %                             plot([stopStopSsd(i) stopStopSsd(i)], [0 maxStopStopSdf(a)*.8], '-r', 'linewidth', 3)
                            %                             plot([stopStopSsd(i)+iSsrt stopStopSsd(i)+iSsrt], [0 maxStopStopSdf(a)*.8], '--r', 'linewidth', 3)
                            %                             plot([stopStopSsd(i) + maxInd(a), stopStopSsd(i) + maxInd(a)], [0 maxStopStopSdf(a)], '-k', 'linewidth', 3)
                            %                             plot([stopStopSsd(i) + maxInd(a) + aHalfMaxInd, stopStopSsd(i) + maxInd(a) + aHalfMaxInd], [0 aHalfMaxStopStopSdf], '-k', 'linewidth', 3)
                            %                             plot([stopStopSsd(i)+(maxInd(a) + aHalfMaxInd/2), stopStopSsd(i)+(maxInd(a) + aHalfMaxInd/2)], [0 maxStopStopSdf(a)*.8], '--k', 'linewidth', 3)
                            % printName = sprintf('%s_canceled_Coh%s_Ssd%d_Trial%d.pdf',Unit(kUnitIndex).name, num2str(stopStopCoh(i)*100), stopStopSsd(i), a);
                            %             print(gcf,fullfile(local_figure_path, subjectID, 'go_vs_canceled', options.ssrt, sessionID,  printName),'-dpdf', '-r300')
                            
                            
                            % To calculate cancel time, use half the time between max and "half
                            % max", relative to ssrt
                            iCancelTime(a) = maxInd(a) + aHalfMaxInd/2 - iSsrt;
                        end
                        
                    end
                    
                    
                end
                % ------------------------------------------------------------------------
                
                
                cancelTimeDist{i} = floor(iCancelTime);
                cancelTimeSdf(i) = floor(iCancelTimeSdf);
                
                
                fprintf('ssd: %d  \tcolor: %.2f  \tgo v. stop: %.1f  %.1f sp, p = %.2f\t canceltime2Std = %d\n',...
                    stopStopSsd(i), stopStopCoh(i), mean(goTargSlowSpike{i}), mean(stopStopSpike{i}), p, cancelTime2Std(i));
                fprintf('ssd: %d  \tcolor: %.2f  \tgo v. stop: %.1f  %.1f sp, p = %.2f\t canceltimeSdf = %d\n',...
                    stopStopSsd(i), stopStopCoh(i), mean(goTargSlowSpike{i}), mean(stopStopSpike{i}), p, cancelTimeSdf(i));
                
                
                if plotFlag && i <= nPlot
                    figure(stopStopFig)
                    
                    % Data aligned on checkerboard onset
                    ax(i, colChkr) = axes('units', 'centimeters', 'position', [xAxesPosition(i, colChkr) yAxesPosition(i, colChkr) axisWidth axisHeight]);
                    set(ax(i, colChkr), 'ylim', [yLimMin yLimMax], 'xlim', [epochRangeChecker(1) epochRangeChecker(end)])
                    cla
                    hold(ax(i, colChkr), 'on')
                    plot(ax(i, colChkr), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
                    ttl = sprintf('SSD: %d  pMag: %.2f  nStop: %d', stopStopSsd(i), stopStopCoh(i), size(iStopStopChecker.signal, 1));
                    title(ttl)
                    
                    % Data aligned on response onset
                    ax(i, colSacc) = axes('units', 'centimeters', 'position', [xAxesPosition(i, colSacc) yAxesPosition(i, colSacc) axisWidth/2 axisHeight]);
                    set(ax(i, colSacc), 'ylim', [yLimMin yLimMax], 'xlim', [epochRangeSacc(1) epochRangeSacc(end)])
                    set(ax(i, colSacc), 'yticklabel', [], 'ycolor', [1 1 1])
                    cla
                    hold(ax(i, colSacc), 'on')
                    plot(ax(i, colSacc), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
                    
                    
                    
                    plot(ax(i, colChkr), [stopStopSsd(i), stopStopSsd(i)], [yLimMin yLimMax], 'color', [.2 .2 .2], 'linewidth', 1)
                    plot(ax(i, colChkr), [stopStopSsd(i) + iSsrt, stopStopSsd(i) + iSsrt], [yLimMin yLimMax], '--', 'color', [0 0 0], 'linewidth', 1)
                    
                    iGoTargRTMean = round(nanmean(goTargSlowCheckerEventLat{i}));
                    
                    % Checkerboard onset aligned
                    iGoTargSlowCheckerFn = goTargSlowCheckerFn{i}(goTargSlowCheckerAlign(i) + epochRangeChecker);
                    iStopStopCheckerFn = stopStopCheckerFn{i}(stopStopCheckerAlign(i) + epochRangeChecker);
                    
                    % Saccade-aligned
                    iGoTargSlowSaccFn = goTargSlowSaccFn{i}(goTargSlowSaccAlign(i) + epochRangeSacc);
                    
                    
                    if options.plotSem
                        % GO Checker aligned
                        iGoTargSlowCheckerFnSem = goTargSlowCheckerFnSem{i}(goTargSlowCheckerAlign(i) + epochRangeChecker);
                        fillX = [epochRangeChecker(1) : epochRangeChecker(end) epochRangeChecker(end) : -1 : epochRangeChecker(1)];
                        fillY = [iGoTargSlowCheckerFn+iGoTargSlowCheckerFnSem fliplr(iGoTargSlowCheckerFn-iGoTargSlowCheckerFnSem)];
                        fill(ax(i, colChkr), fillX, fillY, cMap(pSignalArray == stopStopCoh(i),:), 'facealpha', faceAlpha, 'linestyle' , 'none')
                        
                        % STOP STOP Checker aligned
                        iStopStopCheckerFnSem = stopStopCheckerFnSem{i}(stopStopCheckerAlign(i) + epochRangeChecker);
                        fillX = [epochRangeChecker(1) : epochRangeChecker(end) epochRangeChecker(end) : -1 : epochRangeChecker(1)];
                        fillY = [iStopStopCheckerFn+iStopStopCheckerFnSem fliplr(iStopStopCheckerFn-iStopStopCheckerFnSem)];
                        fill(ax(i, colChkr), fillX, fillY, stopColorMap(pSignalArray == stopStopCoh(i),:), 'facealpha', faceAlpha, 'linestyle' , 'none')
                        
                        % GO Sacc aligned
                        iGoTargSlowSaccFnSem = goTargSlowSaccFnSem{i}(goTargSlowSaccAlign(i) + epochRangeSacc);
                        fillX = [epochRangeSacc(1) : epochRangeSacc(end) epochRangeSacc(end) : -1 : epochRangeSacc(1)];
                        fillY = [iGoTargSlowSaccFn+iGoTargSlowSaccFnSem fliplr(iGoTargSlowSaccFn-iGoTargSlowSaccFnSem)];
                        fill(ax(i, colSacc), fillX, fillY, cMap(pSignalArray == stopStopCoh(i),:), 'facealpha', faceAlpha, 'linestyle' , 'none')
                    end
                    
                    
                    % Checkerboard onset aligned
                    plot(ax(i, colChkr), epochRangeChecker, iGoTargSlowCheckerFn, 'color', cMap(pSignalArray == stopStopCoh(i),:), 'linewidth', goLineW)
                    if iGoTargRTMean < length(iGoTargSlowCheckerFn)
                        plot(ax(i, colChkr), iGoTargRTMean, yLimMin, '.k','markersize', markSize)
                    end
                    
                    plot(ax(i, colChkr), epochRangeChecker, iStopStopCheckerFn, '-', 'color', stopColorMap(pSignalArray == stopStopCoh(i),:), 'linewidth', stopLineW)
                    %                 plot(ax(i, colChkr), epochRangeChecker, sdfDiff(sdfDiffCheckerOn-1:end), 'color', 'g', 'linewidth', stopLineW)
                    
                    % Cancel time
                    if ~isnan(cancelTime2Std(i)) && cancelTime2Std(i)+stopStopSsd(i)+iSsrt < epochRangeChecker(end)
                        plot(ax(i, colChkr), [cancelTime2Std(i)+stopStopSsd(i)+iSsrt cancelTime2Std(i)+stopStopSsd(i)+iSsrt], [yLimMin .95*iStopStopCheckerFn(cancelTime2Std(i)+stopStopSsd(i)+iSsrt)], 'color', 'b', 'linewidth', stopLineW)
                    end
                    if ~isnan(cancelTimeSdf(i)) && cancelTimeSdf(i)+stopStopSsd(i)+iSsrt < epochRangeChecker(end)
                        plot(ax(i, colChkr), [cancelTimeSdf(i)+stopStopSsd(i)+iSsrt cancelTimeSdf(i)+stopStopSsd(i)+iSsrt], [yLimMin .95*iStopStopCheckerFn(cancelTimeSdf(i)+stopStopSsd(i)+iSsrt)], '--', 'color', 'b', 'linewidth', stopLineW)
                    end
                    
                    % Saccade-aligned
                    plot(ax(i, colSacc), epochRangeSacc, iGoTargSlowSaccFn, 'color', cMap(pSignalArray == stopStopCoh(i),:), 'linewidth', goLineW)
                    
                    
                    
                    
                    % Plot  each single valid condition
                    if plotFlag && options.plotSingle
                        figure(stopStopFigSingle)
                        
                        clf
                        
                        % Data aligned on checkerboard onset
                        axSingle(colChkr) = axes('units', 'centimeters', 'position', [xAxesPos(colChkr) yAxesPos(colChkr) axisW axisH]);
                        set(axSingle(colChkr), 'ylim', [yLimMin yLimMax], 'xlim', [epochRangeChecker(1) epochRangeChecker(end)])
                        cla
                        hold(axSingle(colChkr), 'on')
                        plot(axSingle(colChkr), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
                        ttl = sprintf('SSD: %d  pMag: %.2f  nStop: %d', stopStopSsd(i), stopStopCoh(i), size(iStopStopChecker.signal, 1));
                        title(ttl)
                        
                        % Data aligned on response onset
                        axSingle(1, colSacc) = axes('units', 'centimeters', 'position', [xAxesPos(1, colSacc) yAxesPos(1, colSacc) axisW/2 axisH]);
                        set(axSingle(1, colSacc), 'ylim', [yLimMin yLimMax], 'xlim', [epochRangeSacc(1) epochRangeSacc(end)])
                        cla
                        hold(axSingle(1, colSacc), 'on')
                        plot(axSingle(1, colSacc), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
                        ttl = sprintf('SSD: %d  pMag: %.2f  nStop: %d', stopStopSsd(i), stopStopCoh(i), size(iStopStopChecker.signal, 1));
                        title(ttl)
                        
                        
                        plot(axSingle(colChkr), [stopStopSsd(i), stopStopSsd(i)], [yLimMin yLimMax], 'color', [.2 .2 .2], 'linewidth', 1)
                        plot(axSingle(colChkr), [stopStopSsd(i) + iSsrt, stopStopSsd(i) + iSsrt], [yLimMin yLimMax], '--', 'color', [0 0 0], 'linewidth', 1)
                        
                        if options.plotSem
                            % GO Checker aligned
                            %                         iGoTargSlowCheckerFnSem = goTargSlowCheckerFnSem{i}(goTargSlowCheckerAlign(i) + epochRangeChecker);
                            fillX = [epochRangeChecker(1) : epochRangeChecker(end) epochRangeChecker(end) : -1 : epochRangeChecker(1)];
                            fillY = [iGoTargSlowCheckerFn+iGoTargSlowCheckerFnSem fliplr(iGoTargSlowCheckerFn-iGoTargSlowCheckerFnSem)];
                            fill(axSingle(colChkr), fillX, fillY, cMap(pSignalArray == stopStopCoh(i),:), 'facealpha', faceAlpha, 'linestyle' , 'none')
                            
                            % STOP STOP Checker aligned
                            %                         iStopStopCheckerFnSem = stopStopCheckerFnSem{i}(stopStopCheckerAlign(i) + epochRangeChecker);
                            fillX = [epochRangeChecker(1) : epochRangeChecker(end) epochRangeChecker(end) : -1 : epochRangeChecker(1)];
                            fillY = [iStopStopCheckerFn+iStopStopCheckerFnSem fliplr(iStopStopCheckerFn-iStopStopCheckerFnSem)];
                            fill(axSingle(colChkr), fillX, fillY, stopColorMap(pSignalArray == stopStopCoh(i),:), 'facealpha', faceAlpha, 'linestyle' , 'none')
                            
                            % GO Sacc aligned
                            %                         iGoTargSlowSaccFnSem = goTargSlowSaccFnSem{i}(goTargSlowSaccAlign(i) + epochRangeSacc);
                            fillX = [epochRangeSacc(1) : epochRangeSacc(end) epochRangeSacc(end) : -1 : epochRangeSacc(1)];
                            fillY = [iGoTargSlowSaccFn+iGoTargSlowSaccFnSem fliplr(iGoTargSlowSaccFn-iGoTargSlowSaccFnSem)];
                            fill(axSingle(1, colSacc), fillX, fillY, cMap(pSignalArray == stopStopCoh(i),:), 'facealpha', faceAlpha, 'linestyle' , 'none')
                        end
                        
                        
                        % Checkerboard onset aligned
                        plot(axSingle(colChkr), epochRangeChecker, iGoTargSlowCheckerFn, 'color', cMap(pSignalArray == stopStopCoh(i),:), 'linewidth', goLineW)
                        if iGoTargRTMean < length(iGoTargSlowCheckerFn)
                            plot(axSingle(colChkr), iGoTargRTMean, yLimMin, '.k','markersize', markSize)
                        end
                        
                        plot(axSingle(colChkr), epochRangeChecker, iStopStopCheckerFn, '-', 'color', stopColorMap(pSignalArray == stopStopCoh(i),:), 'linewidth', stopLineW)
                        %                 plot(ax(i, colChkr), epochRangeChecker, sdfDiff(sdfDiffCheckerOn-1:end), 'color', 'g', 'linewidth', stopLineW)
                        
                        % Cancel time
                        if ~isnan(cancelTime2Std(i)) && cancelTime2Std(i)+stopStopSsd(i)+iSsrt < epochRangeChecker(end)
                            plot(axSingle(colChkr), [cancelTime2Std(i)+stopStopSsd(i)+iSsrt cancelTime2Std(i)+stopStopSsd(i)+iSsrt], [yLimMin .95*iStopStopCheckerFn(cancelTime2Std(i)+stopStopSsd(i)+iSsrt)], 'color', 'b', 'linewidth', stopLineW)
                        end
                        if ~isnan(cancelTimeSdf(i)) && cancelTimeSdf(i)+stopStopSsd(i)+iSsrt < epochRangeChecker(end)
                            plot(axSingle(colChkr), [cancelTimeSdf(i)+stopStopSsd(i)+iSsrt cancelTimeSdf(i)+stopStopSsd(i)+iSsrt], [yLimMin .95*iStopStopCheckerFn(cancelTimeSdf(i)+stopStopSsd(i)+iSsrt)], '--', 'color', 'b', 'linewidth', stopLineW)
                        end
                        
                        % Saccade-aligned
                        iGoTargSlowSaccFn = goTargSlowSaccFn{i}(goTargSlowSaccAlign(i) + epochRangeSacc);
                        plot(axSingle(1, colSacc), epochRangeSacc, iGoTargSlowSaccFn, 'color', cMap(pSignalArray == stopStopCoh(i),:), 'linewidth', goLineW)
                        
                        
                        if printPlot
                            if ~isdir(fullfile(local_figure_path, subjectID, 'go_vs_canceled', options.ssrt, options.category, sessionID))
                                mkdir(fullfile(local_figure_path, subjectID, 'go_vs_canceled', options.ssrt, options.category, sessionID))
                            end
                            printName = sprintf('%s_canceled_Coh%s_Ssd%d.pdf',Unit(kUnitIndex).name, num2str(stopStopCoh(i)*100), stopStopSsd(i));
                            
                            print(gcf,fullfile(local_figure_path, subjectID, 'go_vs_canceled', options.ssrt, options.category, sessionID,  printName),'-dpdf', '-r300')
                        end
                        
                        
                    end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                end
                
                
            end
        end
        
        %             $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        %             $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        %             $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        %%
        %               Noncanceled Stop vs latency-matched (fast) Go:
        if ANALYZE_NONCANCELED
            if isempty(stopTargCoh)
                fprintf('Session %s doesn"t have enough noncanceled trials in any condiiton to analyze\n', sessionID)
            else
                fprintf('\nGo vs. Noncanceled:\n')
                for i = 1 : length(stopTargCoh)
                    
                    
                    
                    switch options.ssrt
                        case 'intWeightPerSession'
                            iSsrt = round(nanmean(dataInh.ssrtIntegrationWeighted));
                        case 'intPerSsd'
                            % Find the ssrt to use
                            iCohSsrtInd = dataInh.pSignalArray == stopTargCoh(i);
                            iSsdSsrtInd = dataInh.ssd{iCohSsrtInd} == stopTargSsd(i);
                            iSsrt = dataInh.ssrtIntegration{iCohSsrtInd}(iSsdSsrtInd);
                    end
                    if isnan(iSsrt)
                        continue
                    end
                    stopTargSsrt(i) = iSsrt;
                    
                    
                    % Get the go trial data: these need to be split to latency-match with
                    % the stop trial data
                    opt                 = options; % Get default options structure
                    
                    opt.epochName       = epochName;
                    opt.eventMarkName   = markEvent;
                    opt.conditionArray  = {'goTarg'};
                    opt.colorCohArray   = stopTargCoh(i);
                    opt.ssdArray        = [];
                    iGoTargChecker      = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
                    
                    opt.epochName       = 'responseOnset';
                    opt.eventMarkName   = 'checkerOn';
                    iGoTargSacc         = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
                    
                    opt.epochName       = 'targOn';
                    opt.eventMarkName   = 'fixWindowEntered';
                    iGoTargTargOn         = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
                    
                    
                    
                    
                    % Get the stop trial data
                    opt.epochName       = epochName;
                    opt.eventMarkName   = markEvent;
                    opt.conditionArray  = {'stopTarg'};
                    opt.ssdArray        = stopTargSsd(i);
                    iStopTargChecker    = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
                    
                    opt.conditionArray  = {'stopStop'};
                    iStopStopChecker       = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
                    
                    opt.epochName       = 'responseOnset';
                    opt.eventMarkName   = 'checkerOn';
                    opt.conditionArray  = {'stopTarg'};
                    iStopTargSacc       = ccm_concat_neural_conditions(Unit(kUnitIndex), opt);
                    
                    
                    
                    
                    % Use a latency-matching method to get a (fast)
                    % subsample of go trials that match the noncanceled
                    % stop trials
                    switch latencyMatchMethod
                        case 'ssrt'
                            iStopLatency = stopTargSsd(i) + iSsrt;
                            iGoFastTrial = iGoTargChecker.eventLatency <= iStopLatency;
                        case 'mean'
                            iGoTargRT = sort(iGoTargChecker.eventLatency);
                            meanStopTargRT = nanmean(iStopTargChecker.eventLatency);
                            while nanmean(iGoTargRT) > meanStopTargRT
                                iGoTargRT(end) = [];
                            end
                            iStopLatency = iGoTargRT(end);
                            iGoFastTrial = iGoTargChecker.eventLatency <= iStopLatency;
                        case 'match'
                            % Use nearest neighbor method to get rt-matched
                            % trials
                            nStopCorrect = size(iStopStopChecker.signal, 1);
                            data = ccm_match_rt(iGoTargChecker.eventLatency, iStopTargChecker.eventLatency(iStopTargTrial), nStopCorrect);
                            iGoFastTrial = data.goFastTrial;
                    end
                    
                    
                    % Use a subsample of the noncanceled stop RTs that are
                    % later than the SSD plus some time to encode the stimulus
                    if ~options.USE_PRE_SSD
                        iGoFastTrial = iGoFastTrial & iGoTargChecker.eventLatency >  stopTargSsd(i) + encodeTime;
                        iStopTargTrial = iStopTargChecker.eventLatency > stopTargSsd(i) + encodeTime;
                    else
                        % Do nothing for Go Trials, already got 'em above
                        iStopTargTrial = 1 : length(iStopTargChecker.eventLatency);
                    end
                    
                    
                    
                    % ****************************************************************
                    % COLLECT (AND ANALYZE) THE RELEVANT DATA
                    
                    stopTargCheckerData{i}        = iStopTargChecker.signal(iStopTargTrial,:);
                    stopTargCheckerAlign(i)     	= iStopTargChecker.align;
                    stopTargCheckerEventLat{i}	= iStopTargChecker.eventLatency(iStopTargTrial,:);
                    stopTargCheckerSdf{i}        = spike_density_function(stopTargCheckerData{i}, Kernel);
                    
                    stopTargSaccData{i}           = iStopTargSacc.signal(iStopTargTrial,:);
                    stopTargSaccAlign(i)         = iStopTargSacc.align;
                    stopTargSaccEventLat{i}      = iStopTargSacc.eventLatency(iStopTargTrial,:);
                    stopTargSaccSdf{i}        = spike_density_function(stopTargSaccData{i}, Kernel);
                    
                    goTargFastCheckerData{i}      = iGoTargChecker.signal(iGoFastTrial,:);
                    goTargFastCheckerAlign(i) 	= iGoTargChecker.align;
                    goTargFastCheckerEventLat{i}	= iGoTargChecker.eventLatency(iGoFastTrial,:);
                    goTargFastCheckerSdf{i}        = spike_density_function(goTargFastCheckerData{i}, Kernel);
                    
                    goTargFastSaccData{i}         = iGoTargSacc.signal(iGoFastTrial,:);
                    goTargFastSaccAlign(i)       = iGoTargSacc.align;
                    goTargFastSaccEventLat{i}    = iGoTargSacc.eventLatency(iGoFastTrial,:);
                    goTargFastSaccSdf{i}        = spike_density_function(goTargFastSaccData{i}, Kernel);
                    
                    
                    % If we want to normalize (z-score-ish) SDFs,
                    if options.normalize
                        stopTargCheckerSdf{i} = (stopTargCheckerSdf{i} - goTargOnSdfMean) / goTargOnSdfSD;
                        stopTargSaccSdf{i} = (stopTargSaccSdf{i} - goTargOnSdfMean) / goTargOnSdfSD;
                        goTargFastCheckerSdf{i} = (goTargFastCheckerSdf{i} - goTargOnSdfMean) / goTargOnSdfSD;
                        goTargFastSaccSdf{i} = (goTargFastSaccSdf{i} - goTargOnSdfMean) / goTargOnSdfSD;
                    end
                    
                    
                    stopTargCheckerFn{i}      = nanmean(stopTargCheckerSdf{i}, 1);
                    stopTargSaccFn{i}      = nanmean(stopTargSaccSdf{i}, 1);
                    goTargFastCheckerFn{i}      = nanmean(goTargFastCheckerSdf{i}, 1);
                    goTargFastSaccFn{i}         = nanmean(goTargFastSaccSdf{i}, 1);
                    
                    stopTargCheckerFnSem{i}      = std(stopTargCheckerSdf{i}) / sqrt(size(stopTargCheckerSdf{i}, 1));
                    stopTargSaccFnSem{i}      = std(stopTargSaccSdf{i}) / sqrt(size(stopTargSaccSdf{i}, 1));
                    goTargFastCheckerFnSem{i}      = std(goTargFastCheckerSdf{i}) / sqrt(size(goTargFastCheckerSdf{i}, 1));
                    goTargFastSaccFnSem{i}         = std(goTargFastSaccSdf{i}) / sqrt(size(goTargFastSaccSdf{i}, 1));
                    
                    
                    
                    % Hanes et al 1998 t-test of spike rates 40 ms surrounding estimated ssrt
                    stopTargSpike{i}     = sum(iStopTargChecker.signal(iStopTargTrial, spikeWindow + iStopTargChecker.align + stopTargSsd(i) + iSsrt), 2);
                    goTargFastSpike{i}   = sum(iGoTargChecker.signal(iGoFastTrial, spikeWindow + iGoTargChecker.align + stopTargSsd(i) + iSsrt), 2);
                    [h,p,ci,sts]                        = ttest2(stopTargSpike{i}, goTargFastSpike{i});
                    
                    pValue40msStopTarg(i)    = p;
                    stats{i}     = sts;
                    
                    fprintf('ssd: %d  \tcolor: %.2f  \tgo v. stop: %.2f  %.2f sp, p = %.2f\n',...
                        stopTargSsd(i), stopTargCoh(i), mean(goTargFastSpike{i}), mean(stopTargSpike{i}), p);
                    
                    
                    
                    if plotFlag && i <= nPlot
                        figure(stopTargFig)
                        
                        % Data aligned on checkerboard onset
                        ax(i, colChkr) = axes('units', 'centimeters', 'position', [xAxesPosition(i, colChkr) yAxesPosition(i, colChkr) axisWidth axisHeight]);
                        set(ax(i, colChkr), 'ylim', [yLimMin yLimMax], 'xlim', [epochRangeChecker(1) epochRangeChecker(end)])
                        cla
                        hold(ax(i, colChkr), 'on')
                        plot(ax(i, colChkr), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
                        ttl = sprintf('SSD: %d  pMag: %.2f  nStop: %d', stopTargSsd(i), stopTargCoh(i), size(iStopTargChecker.signal, 1));
                        title(ttl)
                        
                        % Data aligned on response onset
                        ax(i, colSacc) = axes('units', 'centimeters', 'position', [xAxesPosition(i, colSacc) yAxesPosition(i, colSacc) axisWidth/2 axisHeight]);
                        set(ax(i, colSacc), 'ylim', [yLimMin yLimMax], 'xlim', [epochRangeSacc(1) epochRangeSacc(end)])
                        set(ax(i, colSacc), 'yticklabel', [], 'ycolor', [1 1 1])
                        cla
                        hold(ax(i, colSacc), 'on')
                        plot(ax(i, colSacc), [1 1], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
                        
                        %     if i > 1
                        %         set(ax(i, colChkr), 'yticklabel', [])
                        %         set(ax(i, colChkr), 'ycolor', [1 1 1])
                        %     end
                        
                        
                        plot(ax(i, colChkr), [stopTargSsd(i), stopTargSsd(i)], [yLimMin yLimMax], 'color', [.2 .2 .2], 'linewidth', 1)
                        plot(ax(i, colChkr), [stopTargSsd(i) + iSsrt, stopTargSsd(i) + iSsrt], [yLimMin yLimMax], '--', 'color', [0 0 0], 'linewidth', 1)
                        
                        iGoTargFastCheckerFn = goTargFastCheckerFn{i}(goTargFastCheckerAlign(i) + epochRangeChecker);
                        iStopTargCheckerFn = stopTargCheckerFn{i}(stopTargCheckerAlign(i) + epochRangeChecker);
                        iGoTargFastSaccFn = goTargFastSaccFn{i}(goTargFastSaccAlign(i) + epochRangeSacc);
                        iStopTargSaccFn = stopTargSaccFn{i}(stopTargSaccAlign(i) + epochRangeSacc);
                        
                        if options.plotSem
                            % GO Checker aligned
                            iGoTargFastCheckerFnSem = goTargFastCheckerFnSem{i}(goTargFastCheckerAlign(i) + epochRangeChecker);
                            fillX = [epochRangeChecker(1) : epochRangeChecker(end) epochRangeChecker(end) : -1 : epochRangeChecker(1)];
                            fillY = [iGoTargFastCheckerFn+iGoTargFastCheckerFnSem fliplr(iGoTargFastCheckerFn-iGoTargFastCheckerFnSem)];
                            fill(ax(i, colChkr), fillX, fillY, cMap(pSignalArray == stopTargCoh(i),:), 'facealpha', faceAlpha, 'linestyle' , 'none')
                            
                            % STOP TARG Checker aligned
                            iStopTargCheckerFnSem = stopTargCheckerFnSem{i}(stopTargCheckerAlign(i) + epochRangeChecker);
                            fillX = [epochRangeChecker(1) : epochRangeChecker(end) epochRangeChecker(end) : -1 : epochRangeChecker(1)];
                            fillY = [iStopTargCheckerFn+iStopTargCheckerFnSem fliplr(iStopTargCheckerFn-iStopTargCheckerFnSem)];
                            fill(ax(i, colChkr), fillX, fillY, stopColorMap(pSignalArray == stopTargCoh(i),:), 'facealpha', faceAlpha, 'linestyle' , 'none')
                            
                            % GO Sacc aligned
                            iGoTargFastSaccFnSem = goTargFastSaccFnSem{i}(goTargFastSaccAlign(i) + epochRangeSacc);
                            fillX = [epochRangeSacc(1) : epochRangeSacc(end) epochRangeSacc(end) : -1 : epochRangeSacc(1)];
                            fillY = [iGoTargFastSaccFn+iGoTargFastSaccFnSem fliplr(iGoTargFastSaccFn-iGoTargFastSaccFnSem)];
                            fill(ax(i, colSacc), fillX, fillY, cMap(pSignalArray == stopTargCoh(i),:), 'facealpha', faceAlpha, 'linestyle' , 'none')
                            
                            % STOP TARG Sacc aligned
                            iStopTargSaccFnSem = stopTargSaccFnSem{i}(stopTargSaccAlign(i) + epochRangeSacc);
                            fillX = [epochRangeSacc(1) : epochRangeSacc(end) epochRangeSacc(end) : -1 : epochRangeSacc(1)];
                            fillY = [iStopTargSaccFn+iStopTargSaccFnSem fliplr(iStopTargSaccFn-iStopTargSaccFnSem)];
                            fill(ax(i, colSacc), fillX, fillY, stopColorMap(pSignalArray == stopTargCoh(i),:), 'facealpha', faceAlpha, 'linestyle' , 'none')
                            
                            
                            plot(ax(i, colChkr), epochRangeChecker, iGoTargFastCheckerFn, 'color', cMap(pSignalArray == stopTargCoh(i),:), 'linewidth', goLineW)
                            iGoTargRTMean = round(nanmean(goTargFastCheckerEventLat{i}));
                            if iGoTargRTMean <= length(iGoTargFastCheckerFn)
                                plot(ax(i, colChkr), iGoTargRTMean, yLimMin, '.k','markersize', markSize)
                            end
                            
                            
                            plot(ax(i, colChkr), epochRangeChecker, iStopTargCheckerFn, '-', 'color', stopColorMap(pSignalArray == stopTargCoh(i),:), 'linewidth', stopLineW)
                            iStopTargRTMean = round(nanmean(stopTargCheckerEventLat{i}));
                            if iStopTargRTMean <= length(iStopTargCheckerFn)
                                plot(ax(i, colChkr), iStopTargRTMean, yLimMin, '.k','markersize', markSize)
                            end
                            
                            plot(ax(i, colSacc), epochRangeSacc, iGoTargFastSaccFn, 'color', cMap(pSignalArray == stopTargCoh(i),:), 'linewidth', goLineW)
                            
                            plot(ax(i, colSacc), epochRangeSacc, iStopTargSaccFn, 'color', stopColorMap(pSignalArray == stopTargCoh(i),:), 'linewidth', stopLineW)
                            
                            
                        end
                        
                    end
                end
                
            end
        end
        
        
        
        
        %             $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        %             $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        %             $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        %             $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        
        
        
        % Collect the data for later analyses
        Data(kUnitIndex).rf      = rf;
        
        % SSRT values
        Data(kUnitIndex).stopStopSsrt      = stopStopSsrt;
        Data(kUnitIndex).stopTargSsrt      = stopTargSsrt;
        
        %               Nonanceled Stop vs latency-matched (Fast) Go:
        Data(kUnitIndex).stopStopCoh      = stopStopCoh;
        Data(kUnitIndex).stopStopSsd      = stopStopSsd;
        Data(kUnitIndex).nStopStop      = nStopStop(usableStopStop);
        Data(kUnitIndex).stopTargCoh      = stopTargCoh;
        Data(kUnitIndex).stopTargSsd      = stopTargSsd;
        Data(kUnitIndex).nStopTarg      = nStopTarg(usableStopTarg);
        
        Data(kUnitIndex).stopTargSpike      = stopTargSpike;
        Data(kUnitIndex).stopTargCheckerData        = stopTargCheckerData;
        Data(kUnitIndex).stopTargCheckerAlign        = stopTargCheckerAlign;
        Data(kUnitIndex).stopTargCheckerEventLat        = stopTargCheckerEventLat;
        Data(kUnitIndex).stopTargSaccData  	= stopTargSaccData;
        Data(kUnitIndex).stopTargSaccAlign        = stopTargSaccAlign;
        Data(kUnitIndex).stopTargSaccEventLat        = stopTargSaccEventLat;
        
        Data(kUnitIndex).goTargFastSpike    = goTargFastSpike;
        Data(kUnitIndex).goTargFastCheckerData      = goTargFastCheckerData;
        Data(kUnitIndex).goTargFastCheckerAlign      = goTargFastCheckerAlign;
        Data(kUnitIndex).goTargFastCheckerEventLat     = goTargFastCheckerEventLat;
        Data(kUnitIndex).goTargFastSaccData 	= goTargFastSaccData;
        Data(kUnitIndex).goTargFastSaccAlign      = goTargFastSaccAlign;
        Data(kUnitIndex).goTargFastSaccEventLat     = goTargFastSaccEventLat;
        
        Data(kUnitIndex).pValue40msStopTarg     = pValue40msStopTarg;
        
        % Jeff's new analysis
        Data(kUnitIndex).pValuePostSsrtRt     = pValuePostSsrtRt;
        Data(kUnitIndex).pValuePostSsrt100     = pValuePostSsrt100;
        Data(kUnitIndex).targToStopDiffSD     = targToStopDiffSD;
        Data(kUnitIndex).preSaccDiffSD     = preSaccDiffSD;
        
        
        %               Canceled Stop vs latency-matched (Slow) Go:
        Data(kUnitIndex).stopStopCheckerData        = stopStopCheckerData;
        Data(kUnitIndex).stopStopCheckerAlign        = stopStopCheckerAlign;
        Data(kUnitIndex).stopStopCheckerEventLat        = stopStopCheckerEventLat;
        Data(kUnitIndex).stopStopCheckerSdf        = stopStopCheckerFn;
        
        Data(kUnitIndex).goTargSlowCheckerData      = goTargSlowCheckerData;
        Data(kUnitIndex).goTargSlowCheckerAlign      = goTargSlowCheckerAlign;
        Data(kUnitIndex).goTargSlowCheckerEventLat     = goTargSlowCheckerEventLat;
        Data(kUnitIndex).goTargSlowCheckerSdf     = goTargSlowCheckerFn;
        
        Data(kUnitIndex).goTargSlowSaccData 	= goTargSlowSaccData;
        Data(kUnitIndex).goTargSlowSaccAlign      = goTargSlowSaccAlign;
        Data(kUnitIndex).goTargSlowSaccEventLat     = goTargSlowSaccEventLat;
        Data(kUnitIndex).goTargSlowSaccSdf     = goTargSlowSaccFn;
        
        Data(kUnitIndex).stopStopSpike      = stopStopSpike;
        Data(kUnitIndex).goTargSlowSpike    = goTargSlowSpike;
        
        
        Data(kUnitIndex).stopTargCoh    = stopTargCoh;
        Data(kUnitIndex).stopTargSsd    = stopTargSsd;
        Data(kUnitIndex).stopTargCond    = stopTargCond;
        Data(kUnitIndex).stopStopCoh    = stopStopCoh;
        Data(kUnitIndex).stopStopSsd    = stopStopSsd;
        Data(kUnitIndex).stopStopCond    = stopStopCond;
        
        
        Data(kUnitIndex).cancelTimeSdf     = cancelTimeSdf;
        Data(kUnitIndex).cancelTimeDist     = cancelTimeDist;
        Data(kUnitIndex).pValue40msStopStop     = pValue40msStopStop;
        Data(kUnitIndex).cancelTime2Std    = cancelTime2Std;
        Data(kUnitIndex).cancelTime4Std    = cancelTime4Std;
        Data(kUnitIndex).cancelTime6Std    = cancelTime6Std;
        
        Data(kUnitIndex).inhibition    = dataInh;
        Data(kUnitIndex).pSignalArray    = pSignalArray;
        
        
        
        if printPlot
            if ANALYZE_CANCELED
                if ~isdir(fullfile(local_figure_path, subjectID, 'go_vs_canceled', options.ssrt, options.category))
                    mkdir(fullfile(local_figure_path, subjectID, 'go_vs_canceled', options.ssrt, options.category))
                end
                print(stopStopFig,fullfile(local_figure_path, subjectID, 'go_vs_canceled', options.ssrt, options.category, [sessionID, '_ccm_go_vs_canceled_',Unit(kUnitIndex).name, '.pdf']),'-dpdf', '-r300')
            end
            
            if ANALYZE_NONCANCELED
                if ~isdir(fullfile(local_figure_path, subjectID, 'go_vs_noncanceled', options.ssrt, options.category))
                    mkdir(fullfile(local_figure_path, subjectID, 'go_vs_noncanceled', options.ssrt, options.category))
                end
                print(stopTargFig,fullfile(local_figure_path, subjectID, 'go_vs_noncanceled', options.ssrt, options.category, [sessionID, '_ccm_go_vs_noncanceled_',Unit(kUnitIndex).name, '.pdf']),'-dpdf', '-r300')
            end
        end
        
        
    end % for kUnitIndex = 1 : nUnit
    
end % function



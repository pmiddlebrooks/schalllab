function Unit = maskbet_single_neuron(subjectID, sessionID, options)

%
% function Unit = maskbet_single_neuron(subjectID, sessionID, options)
%
% Single neuron analyses for reverse masking metacog task. Only plots the
% sdfs. To see rasters, use maskbet_single_neuron_rasters, which displays all
% conditions in a given epoch
%
% If called without any arguments, returns a default options structure.
% If options are input but one is not specified, it assumes default.
%
% input:
%   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%
% Possible options are (default listed first):
%     options.decOutcome       `    = 'collapse','each','correct','incorrect';
%     options.betOutcome            = 'collapse','each','high','low';
%     options.soa                 = 'collapse','each',<a vector of soa values>','mem';
%     options.decCollapseDir         = 'none','leftRight','upDown','all';
%     options.betCollapseDir         = 'none','leftRight','upDown','all';
%     options.directionStage      = 'each','decision','bet';
%     options.matchErrorDir       = 'false','true';
%     options.targetOrSaccadeDir  = 'target','saccade'


% Set default options or return a default options structure
if nargin < 3
    % Conditions to plot/analyze
    options.decOutcome        	= 'each';
    options.betOutcome        	= 'each';
    options.soa                 = 'collapse';
    options.decCollapseDir         = 'none';
    options.betCollapseDir         = 'none';
    options.directionStage      = 'each';
    options.matchErrorDir       = 'false';
    options.targetOrSaccadeDir  = 'target';
    
    options.plotFlag            = 'true';
    options.printFlot           = 'false';
    options.figureHandle      	= 23;
    
    % Return just the default options struct if no input
    if nargin == 0
        Unit           = options;
        return
    end
end


% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);

% For now, use the first reward click as the alignment for reward
trialData.rewardOn = cellfun(@(x) x(1), trialData.rewardOn);

maskAngle = SessionData.maskAngle;
betAngle = SessionData.betAngle;
nDecAngle = length(maskAngle); % Usually will be 4, sometimes 2
nBetAngle = length(betAngle); % For pitt data, always 2 -- but vandy tempo task allows 4 locations (2 each trial that can randomly vary across trials)


leftTargInd = (maskAngle < -89) & (maskAngle > -270) | ...
    (maskAngle > 90) & (maskAngle < 269);
decAngleLeftUp = max(maskAngle(leftTargInd));
decAngleLeftDown = min(maskAngle(leftTargInd));
decAngleRightUp = max(maskAngle(~leftTargInd));
decAngleRightDown = min(maskAngle(~leftTargInd));


clear Unit
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   constants
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MIN_RT = 80;
MAX_RT = 1200;
STD_MULTIPLE = 10; % This is generous, but that's ok because it's not a choice RT task (so there won't be high variation in RTs)

Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay  = 20;
% Kernel.method = 'gaussian';
% Kernel.sigma = 10;


% Which SOAs to analyze
if strcmp(options.soa, 'collapse')
    nSOA = 1;
    soaArray = 'collapse';
elseif strcmp(options.soa, 'each')
    nSOA = 4;
    soaArray = ExtraVar.soaArray;
elseif strcmp(options.soa, 'mem')
    nSOA = 1;
    soaArray = ExtraVar.soaArray(end);
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
end
switch options.betCollapseDir
    case {'leftRight','upDown'}
        nBetAngle = 2;
    case {'all'}
        nBetAngle = 1;
end



epochArray = {'decFixWindowEntered', 'decTargOn', 'decResponseOnset', 'betFixWindowEntered', 'betTargOn', 'betResponseOnset', 'rewardOn'};
nEpoch = length(epochArray);




% How many units were recorded?
% nUnit = size(unitArray, 2);           For nebby and shuffles, always just
% 1 unit recorded. Adjust this for multiple units later



% Get rid of trials with outlying RTs
[allRT, decRtOutlierTrial] = truncate_rt(trialData.decRT, MIN_RT, MAX_RT, STD_MULTIPLE);
trialData.decRT(decRtOutlierTrial) = nan;
[allRT, betRtOutlierTrial] = truncate_rt(trialData.betRT, MIN_RT, MAX_RT, STD_MULTIPLE);
trialData.betRT(betRtOutlierTrial) = nan;












% [minTrialDuration,i] = min(trialData.trialDuration(~strcmp(trialData.trialOutcome, 'fixationAbort')));
% [maxTrialDuration,i] = max(trialData.trialDuration(~strcmp(trialData.trialOutcome, 'fixationAbort')));
% nCondition = 5;
% sdfMax = zeros(nUnit, nEpoch, nSignal, nCondition);  % Keep track of maximum sdf values, for setting y-axis limits in plots

% for kUnit = 1 : nUnit
%     [a, kUnit] = ismember(unitArray{kUnit}, SessionData.spikeUnitArray);
kUnit = 1;

Unit(kUnit).subjectID = subjectID;
Unit(kUnit).sessionID = sessionID;
%     Unit(kUnit).name = unitArray{kUnit};
Unit(kUnit).soaArray = soaArray;







for iSOA = 1 : nSOA
    
        if strcmp(options.soa, 'collapse')
            selectOpt.soa = soaArray;
        else
            selectOpt.soa = soaArray(iSOA);
    end
    
    
    for jEpoch = 1 : length(epochArray)
        jEpochName = epochArray{jEpoch};
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         COLLECT THE RELEVANT TRIALS (based on options input)
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        % If (1) we're only concerned about where decisions were or (2) we're
        % separating decision and bet directions and it's a decision stage
        % epoch, use decision directions for the data and collapse across
        % bet directions
        if strcmp(options.directionStage, 'decision') || ...
                (strcmp(options.directionStage, 'each') && jEpoch <= 4)
            
            selectOpt.betSaccDir = 'collapse';
            selectOpt.betHighDir = 'collapse';
            
            % **************************************************
            % Decision Stage data
            % **************************************************
            for m = 1 : nDecAngle
                
                
                % Establish which decision stage target or saccade angle(s) to select:
                % If collapsing left/right, make dAngle include all angles in
                % hemifield.
                % Else dAngle is a single angle each iteration
                switch options.decCollapseDir
                    case 'none'
                        dAngle = maskAngle(m);
                    case 'leftRight'
                        if m == 1
                            dAngle = [decAngleLeftUp, decAngleLeftDown];
                        else
                            dAngle = [decAngleRightUp, decAngleRightDown];
                        end
                    case 'upDown'
                        if m == 1
                            dAngle = [decAngleLeftUp, decAngleRightUp];
                        else
                            dAngle = [decAngleLeftDown, decAngleRightDown];
                        end
                    case 'all'
                        dAngle = 'collapse';
                    otherwise
                        error('Not a valid options.decCollapseDir input')
                end
                
                switch options.targetOrSaccadeDir
                    case 'target'
                        selectOpt.decTargDir = dAngle;
                        selectOpt.decSaccDir = 'collapse';
                    case 'saccade'
                        selectOpt.decTargDir = 'collapse';
                        selectOpt.decSaccDir = dAngle;
                end
                
                
                
                
                % Collect each valid outcome and combine them below if
                % necessary
                selectOpt.decOutcome = 'targ';
                selectOpt.betOutcome = 'high';
                selectOpt.decSaccDir = dAngle; % This ensures a valid saccade was made to target (valid as determined by convert_translated_data.m)
                targHighTrial = maskbet_trial_selection(trialData, selectOpt);
                
                selectOpt.betOutcome = 'low';
                targLowTrial = maskbet_trial_selection(trialData, selectOpt);
                
                selectOpt.decOutcome = 'distractor';
                selectOpt.betOutcome = 'high';
                selectOpt.decSaccDir = 'collapse'; % For now accept saccades to any location
                distHighTrial = maskbet_trial_selection(trialData, selectOpt);
                
                selectOpt.betOutcome = 'low';
                distLowTrial = maskbet_trial_selection(trialData, selectOpt);
                
                
                concat_outcome_data
                % Possible categories
                % 	Unit(kUnit).soa(iSOA).(jEpochName).angle(m).all.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targ.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).dist.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).high.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).low.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targHigh.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targLow.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distHigh.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distLow.
                
            end % for dAngle = 1 : nAngle
            
            
            
            
            % If (1) we're only concerned about where bets were or (2) we're
            % separating decision and bet directions and it's a bet stage
            % epoch, use bet directions for the data and collapse across
            % decision directions
        elseif strcmp(options.directionStage, 'bet') || ...
                (strcmp(options.directionStage, 'each') && jEpoch > 4)
            
            selectOpt.decSaccDir = 'collapse';
            selectOpt.decTargDir = 'collapse';
            
            
            % **************************************************
            % Bet Stage data
            % **************************************************
            for m = 1 : nBetAngle
                % Establish which bet stage target or saccade angle(s) to select:
                % If collapsing left/right, make dAngle include all angles in
                % hemifield.
                % Else bAngle is a single angle each iteration
                switch options.betCollapseDir
                    case 'none'
                        bAngle = betAngle(m);
                    case 'leftRight'
                        if m == 1
                            bAngle = [betAngleLeftUp, betAngleLeftDown];
                        else
                            bAngle = [betAngleRightUp, betAngleRightDown];
                        end
                    case 'upDown'
                        if m == 1
                            bAngle = [betAngleLeftUp, betAngleRightUp];
                        else
                            bAngle = [betAngleLeftDown, betAngleRightDown];
                        end
                    case 'all'
                        bAngle = 'collapse';
                    otherwise
                        error('Not a valid options.betCollapseDir input')
                end
                
                % We want all bets made to the bAngle, so select trials with
                % bet saccades there (don't need to specify which target was high or low bet):
                selectOpt.betSaccDir = bAngle;
                
                
                % Debug and step through each of these to make sure the         ****************************************
                % selectOpt structure is right
                
                % Collect each valid outcome and combine them below if
                % necessary
                selectOpt.decOutcome = 'target';
                selectOpt.betOutcome = 'high';
                targHighTrial = maskbet_trial_selection(trialData, selectOpt);
                
                selectOpt.betOutcome = 'low';
                targLowTrial = maskbet_trial_selection(trialData, selectOpt);
                
                selectOpt.decOutcome = 'distractor';
                selectOpt.betOutcome = 'high';
                distHighTrial = maskbet_trial_selection(trialData, selectOpt);
                
                selectOpt.betOutcome = 'low';
                distLowTrial = maskbet_trial_selection(trialData, selectOpt);
                
                
                
                % Combine the data from various outcomes if needed
                if strcmp(options.decOutcome, 'collapse') && strcmp(options.betOutcome, 'collapse')
                    jTrial = [targHighTrial; targLowTrial; distHighTrial; distLowTrial];
                elseif strcmp(options.decOutcome, 'each') && strcmp(options.betOutcome, 'collapse')
                    targTrial = [targHighTrial; targLowTrial];
                    distTrial = [distHighTrial; distLowTrial];
                elseif strcmp(options.decOutcome, 'collapse') && strcmp(options.betOutcome, 'each')
                    highTrial = [targHighTrial; distHighTrial];
                    lowTrial = [targLowTrial; distLowTrial];
                elseif strcmp(options.decOutcome, 'each') && strcmp(options.betOutcome, 'each')
                    % do nothing- data is already divided into ch, cl, ih,
                    % and il trials
                    
                end
                
                concat_outcome_data
                % Possible categories
                %  	Unit(kUnit).soa(iSOA).(jEpochName).angle(m).all.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targ.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).dist.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targHigh.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targLow.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distHigh.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distLow.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).high.
                %     Unit(kUnit).soa(iSOA).(jEpochName).angle(m).low.
            end % b = 1 : nBetAngle
            
            
            
            
            
        end % options.directionStage
        
        
        
        
        
        
    end % for jEpoch = 1 : length(epochArray
end % for iSOA = 1 : nSOA









% Nested function to concatenate the relevant data whether it's in decision
% stage or bet stage
    function concat_outcome_data
        
        % Concatenate the data from various outcomes if needed
        % ****************************************************************
        
        
        % Collapsing across all outcomes
        if strcmp(options.decOutcome, 'collapse') && strcmp(options.betOutcome, 'collapse')
            jTrial = [targHighTrial; targLowTrial; distHighTrial; distLowTrial];
            
            % All
            alignList               = trialData.(jEpochName)(jTrial);
            [alignRas, alignInd]    = spike_to_raster(trialData.spikeData(jTrial, kUnit), alignList);
            
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).all.alignTime   = alignInd;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).all.raster      = alignRas;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).all.sdf         = spike_density_function(alignRas, Kernel);
            %             Unit(kUnit).soa(iSOA).(jEpochName).angle(m).all.sdfMean 	= nanmean(sdf, 1);
            
            
            
            % Collapsing across bets (target vs distrator decisions)
        elseif strcmp(options.decOutcome, 'each') && strcmp(options.betOutcome, 'collapse')
            targTrial = [targHighTrial; targLowTrial];
            distTrial = [distHighTrial; distLowTrial];
            
            % Target
            alignList               = trialData.(jEpochName)(targTrial);
            [alignRas, alignInd]    = spike_to_raster(trialData.spikeData(targTrial, kUnit), alignList);
            
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targ.alignTime   = alignInd;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targ.raster      = alignRas;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targ.sdf         = spike_density_function(alignRas, Kernel);
            %             Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targ.sdfMean 	= nanmean(sdf, 1);
            
            % Distractor
            alignList               = trialData.(jEpochName)(distTrial);
            [alignRas, alignInd]    = spike_to_raster(trialData.spikeData(distTrial, kUnit), alignList);
            
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).dist.alignTime   = alignInd;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).dist.raster      = alignRas;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).dist.sdf         = spike_density_function(alignRas, Kernel);
            %             Unit(kUnit).soa(iSOA).(jEpochName).angle(m).dist.sdfMean 	= nanmean(sdf, 1);
            
            
            % Collapsing across decisions (high vs. low bets)
        elseif strcmp(options.decOutcome, 'collapse') && strcmp(options.betOutcome, 'each')
            highTrial = [targHighTrial; distHighTrial];
            lowTrial = [targLowTrial; distLowTrial];
            
            % High
            alignList               = trialData.(jEpochName)(highTrial);
            [alignRas, alignInd]    = spike_to_raster(trialData.spikeData(highTrial, kUnit), alignList);
            
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).high.alignTime   = alignInd;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).high.raster      = alignRas;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).high.sdf         = spike_density_function(alignRas, Kernel);
            %             Unit(kUnit).soa(iSOA).(jEpochName).angle(m).high.sdfMean 	= nanmean(sdf, 1);
            
            % Low
            alignList               = trialData.(jEpochName)(lowTrial);
            [alignRas, alignInd]    = spike_to_raster(trialData.spikeData(lowTrial, kUnit), alignList);
            
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).low.alignTime   = alignInd;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).low.raster      = alignRas;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).low.sdf         = spike_density_function(alignRas, Kernel);
            %             Unit(kUnit).soa(iSOA).(jEpochName).angle(m).low.sdfMean 	= nanmean(sdf, 1);
            
            
            % Comparing all outcomes
        elseif strcmp(options.decOutcome, 'each') && strcmp(options.betOutcome, 'each')
            % Target High
            alignList               = trialData.(jEpochName)(targHighTrial);
            [alignRas, alignInd]    = spike_to_raster(trialData.spikeData(targHighTrial, kUnit), alignList);
            
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targHigh.alignTime   = alignInd;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targHigh.raster      = alignRas;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targHigh.sdf         = spike_density_function(alignRas, Kernel);
            %             Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targHigh.sdfMean 	= nanmean(sdf, 1);
            
            % Target Low
            alignList               = trialData.(jEpochName)(targLowTrial);
            [alignRas, alignInd]    = spike_to_raster(trialData.spikeData(targLowTrial, kUnit), alignList);
            
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targLow.alignTime   = alignInd;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targLow.raster      = alignRas;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targLow.sdf         = spike_density_function(alignRas, Kernel);
            %             Unit(kUnit).soa(iSOA).(jEpochName).angle(m).targLow.sdfMean 	= nanmean(sdf, 1);
            
            % Distractor High
            alignList               = trialData.(jEpochName)(distHighTrial);
            [alignRas, alignInd]    = spike_to_raster(trialData.spikeData(distHighTrial, kUnit), alignList);
            
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distHigh.alignTime   = alignInd;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distHigh.raster      = alignRas;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distHigh.sdf         = spike_density_function(alignRas, Kernel);
            % %             Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distHigh.sdfMean 	= nanmean(sdf, 1);
            
            % Distractor Low
            alignList               = trialData.(jEpochName)(distLowTrial);
            [alignRas, alignInd]    = spike_to_raster(trialData.spikeData(distLowTrial, kUnit), alignList);
            
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distLow.alignTime   = alignInd;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distLow.raster      = alignRas;
            Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distLow.sdf         = spike_density_function(alignRas, Kernel);
            %             Unit(kUnit).soa(iSOA).(jEpochName).angle(m).distLow.sdfMean 	= nanmean(sdf, 1);
            
        end
    end
end


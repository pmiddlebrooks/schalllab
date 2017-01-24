function Data = cmd_concat_neural_conditions(Unit, epochName, eventMarkName, conditionArray, angleArray, ssdArray)

% If singal strength or ssd index vectors are not input, assume user wants to
% collapse across all of them
if nargin < 4
    conditionArray = {'goTarg', 'stopTarg', 'stopStop'};
end
if nargin < 6
    ssdArray = Unit(1).ssdArray;
end
nUnit = length(Unit);
Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;

[c, ssdIndexArray] = ismember(ssdArray, Unit(1).ssdArray);
[c, angleIndexArray] = ismember(angleArray, Unit(1).targAngleArray);
% If multiple units are in the Unit struct, loop through each one.
for k = 1 : nUnit
    % Create a concatenated cell array of rasters according to the desired
    % input condiitons
    rasterCell = {};
    alignCell = {};
%     alignList = [];
    eventLatency = [];
    for c = 1 : length(conditionArray)
        condition = conditionArray{c};
        for i = 1 : length(angleIndexArray)
            angleIndex = angleIndexArray(i);
            
            
            if (strcmp(condition, 'goTarg') && ...
                    ~strcmp(epochName, 'stopSignalOn') && ...
                    ~strcmp(eventMarkName, 'stopSignalOn'))
%                 alignList = [alignList; Unit(k).angle(angleIndex).(condition).(epochName).alignTimeList];
                rasterCell = [rasterCell; Unit(k).angle(angleIndex).(condition).(epochName).raster];
                alignCell = [alignCell; Unit(k).angle(angleIndex).(condition).(epochName).alignTime];
                
                iEventLatency = get_event_latency(Unit(k).angle(angleIndex).(condition), epochName, eventMarkName);
                eventLatency = [eventLatency; iEventLatency];
            elseif strcmp(condition, 'stopTarg') || strcmp(condition, 'stopDist') || strcmp(condition, 'stopStop')
                for j = 1 : length(ssdIndexArray)
                    ssdIndex = ssdIndexArray(j);
                    
                    % Don't include stopStop trials if aligning on
                    % response onset
                    if ~(strcmp(epochName, 'responseOnset') && strcmp(condition, 'stopStop'))
%                         alignList = [alignList; Unit(k).angle(angleIndex).(condition).ssd(ssdIndex).(epochName).alignTimeList];
                        rasterCell = [rasterCell; Unit(k).angle(angleIndex).(condition).ssd(ssdIndex).(epochName).raster];
                        alignCell = [alignCell; Unit(k).angle(angleIndex).(condition).ssd(ssdIndex).(epochName).alignTime];
                        
                        % Don't get response onset event marks for stop correct
                        if ~(strcmp(eventMarkName, 'responseOnset') && strcmp(condition, 'stopStop'))
                            jEventLatency = get_event_latency(Unit(k).angle(angleIndex).(condition).ssd(ssdIndex), epochName, eventMarkName);
                            eventLatency = [eventLatency; jEventLatency];
                        end
                    end
                end % j = 1 : length(ssdIndex)
            end % if strcmp(condition)
        end % i = 1 : length(angleIndex)
    end % for c = 1 : length(conditionArray)
    
    % Get rasters for each set of trials
    [Data(k).raster, Data(k).align]     = align_raster_sets(rasterCell, alignCell);
    [Data(k).raster, eventLatency] = sort_rasters_events(Data(k).raster, eventLatency);
    Data(k).sdf   = nanmean(spike_density_function(Data(k).raster, Kernel), 1);
    Data(k).eventLatency = eventLatency;
    %     Data(k).alignTimeList = alignList;
end

end  % main function






function [eventLatency] = get_event_latency(UnitCondition, epochName, eventMarkName)

eventLatency = [];
alignTime = UnitCondition.(epochName).alignTime;
unitAlignList = UnitCondition.(epochName).alignTimeList;

if ~isempty(eventMarkName)
    eventAlignList = UnitCondition.(eventMarkName).alignTimeList;
else
    switch epochName
        case 'fixWindowEntered'
            % target onset
            eventAlignList = UnitCondition.targOn.alignTimeList;
        case 'targOn'
            % target onset
            eventAlignList = UnitCondition.checkerOn.alignTimeList;
        case 'stopSignalOn'
            % target onset
            eventAlignList = UnitCondition.checkerOn.alignTimeList;
        case 'responseOnset'
            % target onset
            eventAlignList = UnitCondition.checkerOn.alignTimeList;
        case 'rewardOn'
            % target onset
            eventAlignList = UnitCondition.responseOnset.alignTimeList;
    end
end
if ~isempty(alignTime)
    eventLatency = eventAlignList - unitAlignList;
end
end



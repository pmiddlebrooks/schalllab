function single_neuron_reward_outcomes(monkey, brainArea, sessionNumber, spikeUnit, target, plotFlag)

% function nostop_1back_neural(monkey, sessionNumber, brainArea, back1trial, current_trial, alignEvent, epochSpan)
%
%
%
% spikeUnit: PDP, DSP01a, DSP01b, ...
%
% target: 'each' or 'all'
%

rewardContingency = {'earned', 'unexpected', 'omitted'};
rewardColor = {'k', 'b', 'r'};

Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;


% Get the file name and path, joined in "matpath"
sessionNumberString = num2str(sessionNumber, '%03.0f');
dataPath = ['data/', monkey, '_', brainArea, '_', sessionNumberString];
load(dataPath)

% Obtain the spike data from the relevant trials
spikeUnitIndex = find(strcmp(spikeUnit, sessionData.collectedData.neurophysiology.singleUnit.names));
spikeData = trialData.singleUnit(:, spikeUnitIndex);


if plotFlag
    rewAx = axes('units', 'centimeters', 'Position',[12 11, 8, 5]);
    hold all;
end

% for iTarget = 1 : nTarget

alignEvent = 'rewardOnset';
displaySpan = get_display_span(alignEvent);


for iReward = 1 : length(rewardContingency);
    [trialList, alignTimeList] = reward_trial_selection(monkey, brainArea, sessionNumber, 'all', rewardContingency{iReward});
    
    % Get the aligned rasters
    [iRewardRaster, alignmentIndex] = spike_to_raster(spikeData(trialList, :), alignTimeList);
    
    [iRewardSDF, kernel] = spike_density_function(iRewardRaster, Kernel);
    iMeanSDF = nanmean(iRewardSDF, 1);
    if plotFlag
        plot(rewAx,iMeanSDF(alignmentIndex+displaySpan(1) : alignmentIndex+displaySpan(2)), 'color', rewardColor{iReward}, 'linewidth', 2)
        
    end % plotFlag
       
end % iReward




if plotFlag
    set(rewAx, 'xlim', [1 , displaySpan(2) - displaySpan(1)])
    set(rewAx, 'xtick', [1 : 200 : displaySpan(2)-displaySpan(1)], 'xtickLabel', [displaySpan(1)-1 : 200 : displaySpan(2)])
    yMax = ylim;
    
    plot(rewAx, [-displaySpan(1) -displaySpan(1)], [0 yMax(2)], 'k', 'linewidth', 2)
    legend(rewardContingency, 'location', 'northeast')
end


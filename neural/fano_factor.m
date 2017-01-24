function [fanoFactor] = fano_factor(nSpikeArray)

% [fanoFactor] = fano_factor(method, alignedRastersOrSDFs)
%
% Calculates the fano factor for an input matrix of spike rasters
% (alignedRasters: each row is a trial, each column is a millisecond in the
% trial). alginedRasters should already be parsed into the time window
% (i.e. bin) you want to analyze
%
% For now, only coding the fano factor with this script. Will want to
% include the ability to calculate the "normalized variance" a la
% Churchland et al. 2006 NatNeurosci (p. 3699), which is essentially a fano factor
% calculated on the spike rates (spike density funcitons) instead of the
% spike counts.
%
% ... something like the following:
%
% [fanoFactor] = fano_factor(method, alignedRastersOrSDFs)
%
% switch method
%     case 'spike count'
%         alignedRasters = alignedRastersOrSDFs;
%         spikesPerTrial = sum(alignedRasters);
%         spikesMean = mean(spikesPerTrial);
%         spikesVariance = var(spikesPerTrial);
%         
%         fanoFactor = spikesVariance / spikesMean;
%         
%         
%     case 'spike rate'
%         sdfs = alignedRastersOrSDFs;
%         % figure out the constant "c" here (based on how the sdfs were
%         % generated)
%         c = something;
%         sdfMean = mean(sdfs);
%         sdfVariance = var(sdfs);
%         normalizedVariance = c * sdfVariance / sdfMean;
%         fanoFactor = normalizedVariance;
% end


% spikesPerTrial = nansum(alignedRasters);
spikesMean = nanmean(nSpikeArray);
spikesVariance = nanvar(nSpikeArray);

fanoFactor = spikesVariance / spikesMean;





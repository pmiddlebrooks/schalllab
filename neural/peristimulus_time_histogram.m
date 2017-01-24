function [PSTH, binCenters] = peristimulus_time_histogram(alignedRaster, binWidth, plotFlag, spikeNumberOrRate)

% PSTH = peristimulus_time_hisogram(alignedRasters, binWidth)
%
% Creates a vector PSTH from an input matrix (or vector) of single trial
% rasters (or one trial raster)
% (alignedRaster). Size of spikeDensityFunction = size of alignedRasters.
%
% alignedRaster: the matrix of rasters or vector of raster.
%
% binWidth: can be a number denoting the desired width in ms, or the string
% 'optimal', which uses an algorithm (from Shimazaki & Shinomoto 2007) to
% optimize the binwidth for the particular alignedRaster input.
%
% plotFlag: 1 for plot, 0 for no plot
%
% spikeNumberOrRate (optional): input 'number' to obtain a PSTH with the spike counts
% in each bin, or 'rate' to obatain a PSTH with spike rates in each bin.
% Default is 'number.
%
% PSTH can be of the following form:
% --A single trial vector, if alignedRaster is a single trial raster.
% --An average PSTH from multiple trials, if alignedRaster is a matrix of trial rasters
% --A matrix of mulitple single PSTHs, if alignedRaster is a matrix of trial rasters
%
% binCenters: A vector of the centers (x-values) of each bin- can be used
% for plotting the data etc.
%

[nTrial rasterDuration] = size(alignedRaster);


plotCostFunction = 0;
if strcmp(binWidth, 'optimal')
    [binWidth, costFunction] = optimal_binwidth(alignedRaster);
    plotCostFunction = 1;
end

nBin            = floor(rasterDuration / binWidth);
binCenters      = (binWidth/2 : binWidth : nBin*binWidth);
PSTH            = zeros(1, nBin);
% Loop through the bins and add up the spikes in each bin
for iBin = 1 : nBin
    PSTH(iBin) = sum(sum(alignedRaster(:, binWidth*(iBin-1)+1 : binWidth*iBin)));
end
if nargin > 3
    if strcmp(spikeNumberOrRate, 'rate')
        PSTH = PSTH ./ nTrial ./ (binWidth * .001);
    end
end

binCenters = (binWidth/2 : binWidth : nBin*binWidth);

if plotFlag
    bar(binCenters, PSTH, 'hist');
    set(gcf, 'color', [1 1 1])
    xlabel('Time (sec)')
    ylabel('Firing Rate (sp/s)')
end

if plotCostFunction
    figure(2)
    plot(costFunction)
end







    function [binWidth, costFunction] = optimal_binwidth(alignedRaster)
        
        % This attempts to find an optimal binWidth to used. It's
        % algorithm 1 from Shimazaki & Shinomoto 2007 Neural Computation,
        % p.1508
        
        % Initialize the binWidth to the signal duration divided into 10 equal
        % bins
        rasterLength        = size(alignedRaster, 2);
        minimumBins         = 20;
        maximumBinWidth     = rasterLength / minimumBins;
        binWidthRange       = find(mod(maximumBinWidth, 1 : 1:maximumBinWidth) == 0);
        nBinRange           = rasterLength ./ binWidthRange;
        nSolution   = length(binWidthRange);
        costFunction        = zeros(nSolution, 1);
        for iSolutionIndex = 1 : nSolution
            iPSTH = zeros(1, nBinRange(iSolutionIndex));
            % Loop through the bins and add up the spikes in each bin
            for jBin = 1 : nBinRange(iSolutionIndex)
                iPSTH(jBin) = sum(sum(alignedRaster(:, binWidthRange(iSolutionIndex)*(jBin-1)+1 :  binWidthRange(iSolutionIndex)*jBin)));
            end
            %     iPSTH
            %     pause
            iMeanSpike = mean(iPSTH);
            iVarianceSpike = var(iPSTH);
            costFunction(iSolutionIndex) = (2*iMeanSpike - iVarianceSpike) / ((nTrial * binWidthRange(iSolutionIndex)) ^ 2);
        end
        [minCostFunction, minIndex] = min(costFunction);
        binWidth = binWidthRange(minIndex);
    end
end
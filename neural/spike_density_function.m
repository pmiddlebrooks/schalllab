function [spikeDensityFunction, kernelShape] = spike_density_function(alignedRasters, Kernel)

% spikeDensityFunction = spike_density_function(alignedRasters, Kernel)
%
% Creates a matrix of single trial spike density functions
% (spikeDensityFunction) from an input matrix of single trial rasters
% (alignedRasters). Size of spikeDensityFunction = size of alignedRasters.
%
% Kernel: A structure that determines how to convolve spikes in the input rasters:
%
% Kernel.method: can be 'gaussian' (default) or 'postsynaptic potential'
%
% If 'gaussian', Kernel needs the field sigma (gaussian width) or enter
% nothing to accept default sigma (10) value.
%
% Kernel.method = 'gaussian:
% Kernel.sigma = 10 (default).
%
% If 'postsynaptic potential', Kernel needs the growth and decay
% time constant values, or nothing to accept default growth (1) and decay
% (20) values.
%
% Kernel.method = 'postsynaptic potential':
% Kernel.growth = 1 (default).
% Kernel.decay = 20 (default).
%
%
% Gaussian kernel: Richmond et al 1987 JNeurophysiol, p.136
% Postsynaptic potential kernel: Thompson et al 1996 JNeurophysiol, p.4042
%
if isempty(alignedRasters)
    spikeDensityFunction = [];
    kernelShape = [];
    %     disp('spike_density_function.m:   No rasters input to compute a spike density funtion')
    return
end



% Set default values
if nargin < 2
    Kernel.method = 'gaussian';
    Kernel.sigma = 20;
%     Kernel.method = 'postsynaptic potential';
%             Kernel.growth = 1;
%             Kernel.decay = 20;
end

switch Kernel.method
    case 'gaussian'
        if ~isfield(Kernel, 'sigma')
            Kernel.sigma = 10;
        end
    case 'postsynaptic potential'
        if ~isfield(Kernel, 'growth')
            Kernel.growth = 1;
        end
        if ~isfield(Kernel, 'decay')
            Kernel.decay = 20;
        end
    otherwise
        error('Need to enter Kernel.method = ''guassian'' or ''postsynaptic potential''')
        return
end




[nTrial sdfDuration] = size(alignedRasters);
% Initialize spikeDensityFunction as zeros the same size as alignedRasters
spikeDensityFunction = zeros(nTrial, sdfDuration);
% Get the half length of the rasters duration- will use this later to cut
% off the sdf, which will have grown beyond the end (postsynaptic potential
% method) or beyond the beginning and the end (gaussian method), due to the
% convolution process
halfLength = ceil(sdfDuration / 2);
switch Kernel.method
    case 'gaussian'
        k = -halfLength : halfLength;
        kernelShape = normpdf(k, 0, Kernel.sigma);
        kernelCenter = ceil(length(kernelShape)/2);
        kernelShape = kernelShape(kernelCenter-50 : kernelCenter+50);
        
    case 'postsynaptic potential'
        kernelHalfLength = round(Kernel.decay * 8);
        kernelHalfTimes = 0 : kernelHalfLength;
        kernelHalf1 = zeros(1, kernelHalfLength);
        kernelHalf2 = (1 - (exp(-(kernelHalfTimes ./ Kernel.growth)))) .* (exp(-(kernelHalfTimes ./ Kernel.decay)));
        kernelHalf2 = kernelHalf2 ./ sum(kernelHalf2);
        kernelShape = [kernelHalf1, kernelHalf2];
end
% clf
% plot(kernel)
% pause
%         size(kernel)
for iTrial = 1 : nTrial
    convolvedSDF = conv(alignedRasters(iTrial, :), kernelShape);
    
    center = ceil(length(convolvedSDF) / 2 );
    iSDF = convolvedSDF(center - halfLength : (center + halfLength - 1));
    %     plot(iSDF)
    %     size(iSDF)
    % pause
    if length(iSDF) > sdfDuration
        iSDF = iSDF(1 : sdfDuration);
    end;
    spikeDensityFunction(iTrial, :) = iSDF;
end;

% Convert the sdf from spikes/ms to spikes/s
spikeDensityFunction = (spikeDensityFunction .* 1000);

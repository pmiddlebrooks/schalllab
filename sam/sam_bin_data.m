function [rtQ,cumProb,cumProbDefective,probMass,probMassDefective] = sam_bin_data(rt,prop,cumProb,minSize,dt)
% SAM_BIN_DATA Groups RT into bins
%  
% DESCRIPTION 
% RT bin edges are defined based on cumulative probabilities
% 
% SYNTAX
% [rtQ,cumProbDefective,probMass,probMassDefective] = SAM_BIN_DATA(rt,prop,cumProb,minSize)
% 
% rt        - reaction times (Nx1 double)
% prop      - proportion of trials in this category (1x1 double), with
%             category being either go trials or stop trials, e.g.
%             * correct (choice) go trials in the category of go trials 
%             * correct (choice) signal-respond trials in the category of
%             stop trials
% cumProb   - cumulative probabilities for which to compute RT quantiles (1xP double)
% minSize   - minimum number of trials for binning (1x1 or 1x2 double)
%             * if scalar, RT data with fewer than minSize trials are 
%             * if 1x2 double, RT data with fewer trials than 
%               + min(minSize) are grouped into one bin
%               + max(minSize) are grouped into two bins (median-split)
% dt        - time step of the model (1x1 double)
%
% EXAMPLES
% rt        = 400 + 100 * randn(500,1) + exprnd(150,500,1);
% prop      = 0.95;
% cumProb   = .1:.2:.9;
% minSize   = 40;
% dt        = 1;
% [rtQ,cumProbDefect,probMass,probMassDefective] = ...
% SAM_BIN_DATA(rt,prop,cumProb,minSize,1);
%
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Wed 28 Aug 2013 14:40:48 CDT by bram 
% $Modified: Wed 28 Aug 2013 15:28:13 CDT by bram
 
% CONTENTS 
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS AND SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Remove NaNs from RT vector, if any
% =========================================================================

if ~isempty(find(isnan(rt), 1))
    rt              = rt(~isnan(rt));
end

% Downscale to the temporal resolution of the model
rt                  = round(rt./dt).*dt;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. COMPUTE QUANTILES, PROBABILITIES, AND PROBABILITY MASSES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% RT quantiles
% =========================================================================
if isempty(rt)
    rt              = NaN;
    cumProb         = NaN;
    rtQ             = [];
elseif numel(rt) < min(minSize)
    cumProb         = NaN;
    rtQ             = [];
elseif numel(rt) < max(minSize)
    cumProb         = 0.5;
    rtQ             = quantile(rt,cumProb);
else
    rtQ             = quantile(rt,cumProb);
end

% Defective cumulative probabilities
% =========================================================================
cumProbDefective    = cumProb.*prop;

% Probability masses and defective probability masses
% =========================================================================
histCount           = histc(rt,[-Inf,rtQ,Inf]);
histCount           = histCount(1:end-1);
if max(histCount) == 0
    probMass        = histCount;
else
    probMass        = histCount./sum(histCount);
end
probMass            = probMass(:);
probMassDefective   = prop.*probMass;
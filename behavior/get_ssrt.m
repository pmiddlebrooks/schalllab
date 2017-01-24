function ssrt = get_ssrt(ssd,pRespond,nStop,goRT,WeibParams)
% SSRT
%   Find SSRT for a given data set.
%
%
% DESCRIPTION
%   Using the method of integration, the mean method, using mean SSD
%   see Hanes, Patterson, and Schall. JNeurophysiol. 1998.
%
% SYNTAX
%   method = ssrt(ssd,pRespond,nStop,goRT,WeibParams)
%
%
%   INPUT:
%       ssd       = vector containing all stop signal delays
%       pRespond   = p(noncanceled|SSD) for all stop signal delays
%       nStop         = the number of trials at each SSD
%       goRT        = reaction times for no stop trials
%
%   OPTIONAL INPUT:
%       WeibParams = four value vector containing optimum coeffecients that
%                    describe a Weibull curves best fit to the inhibition
%                    function such that:
%                            alpha  = WeibParams(1);
%                            beta   = WeibParams(2);
%                            gamma  = WeibParams(3);
%                            delta  = WeibParams(4);
%                    if no coefficients are provided or if the empty set is
%                    entered ([]) coefficients will be estimated using
%                    Weibull.m
%
%   OUTPUT:
%
%       ssrt: A structure containing:
%       
%           integration:        integration method SSRTs at each SSD 
%           integrationWeighted: weighted (by number of trials) integration
%                               method SSRT
%           integrationSimple:  Bisset & Logan simplified integration method
%           mean:               mean method
%           grand:              All of the above methods averaged together
%
%
%
%
% CONTENTS
% 1.Prepare data
% 2.Integration methods
%   2.1.Classic integration (one ssrt for each ssd)
%   2.2.Weighted integration
%   2.3.Simple integration
% 3. Simple method
% 3. Grand average






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1. Preliminary values setup
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pRespond    = reshape(pRespond, length(pRespond), 1);
nStop       = reshape(nStop, length(nStop), 1);
ssd         = reshape(ssd, length(nStop), 1);


% Sort RTs and get rid of nans
goRT = sort(goRT(~isnan(goRT)));

% Get the mean (and weighted mean) of ssds
ssdMean         = mean(ssd);
ssdMeanWeight   = sum(ssd .* nStop ./ sum(nStop));

% Mean goRT
meanRT = mean(goRT);










% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2. Integration method: Estimate individual ssrt at each SSD, and average
%  them to obtain final ssrt
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.1. Non-weighted average of ssrt (integration method)
% =========================================================================
p_goRT      = 1/length(goRT):1/length(goRT):1; %y-axis of a cumulative prob dist
ssrtEach    = nan(length(ssd), 1);

% Estimate SSRT at each SSD
for i = 1 : length(ssd)
  iSSD    = round(ssd(i));
  iPResp 	= pRespond(i);     
  
  % Don't calculate SSRT for a given SSD if:
  % 1. pNoncanceled doesn't exists (isnan)
  % 2. pNoncanceled is zero and the next SSD's pNoncanceled is also zero
  % 3. pnoncanceled is one and the previous SSD's pNoncanceled is also 1
  if isnan(iPResp) || ...
            iPResp == 0 || ...
      (i < length(ssd) && iPResp == 0 && pRespond(i+1) == 0) || ...
      (i > 1 && iPResp == 1 && pRespond(i-1) == 11 && pRespond(i-2) == 1)
    % ssrtEach(i) = nan
  else
    indexRT     = find(p_goRT >= iPResp,1);   %match estimated p(noncan|SSD) to p(RT)
    ssrtEach(i) = goRT(indexRT) - iSSD; %this is the SSRT by the method of integration
  end
  
end
ssrt.integration = ssrtEach;


% 2.2. Weighted average of ssrt (integration method)
% =========================================================================
ssrtIntWeight                   = nansum(ssrt.integration .* nStop ./ sum(nStop));
ssrt.integrationWeighted = ssrtIntWeight;





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.3 Bisset and Logan simplified integration method:
%   (goRT at the weighted mean probability of responding) - (weighted mean SSD)
% A single SSRT estimate for the entire distribution
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pNoncancelMeanWeight = nansum(pRespond .* nStop ./ sum(nStop));
indexRT     = find(p_goRT >= pNoncancelMeanWeight,1);   %match estimated p(noncan|SSD) to p(RT)

ssrt.integrationSimple = goRT(indexRT) - ssdMeanWeight;






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Mean method:
% A single SSRT for the entire distribution using the mean method, using a
% weibull curve genereated with a weighted SSD distribution
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if Weibull coefficients have not been provided find them
if nargin < 5, WeibParams = Weibull(ssd,pRespond,nStop); end
% if nargin < 5, [WeibParams best_fit] = Weibull_fast(ssd,pRespond,nStop); end

%get params
alpha = WeibParams(1);
beta  = WeibParams(2);
gamma = WeibParams(3);
delta = WeibParams(4);

%generate curve
time = 0:10000; %this needs be be further out than our furthest estimate of SSD using the Weibull fit so that we get the WHOLE distribution
WeibCurve = gamma - ((exp(-((time./alpha).^beta))).*(gamma-delta)); %reconstruct the Weibull curve with the best fit params

%rescale to 0 and 1
WeibCurve = WeibCurve - delta;
WeibCurve = WeibCurve ./ (gamma - delta);

SSDdensity = [nan diff(WeibCurve)] .* time; %convert CDF into PDF
meanSSD    = nansum(SSDdensity);            %integrate under the curve to find mean SSD

% method.mean.ssrt = meanRT - meanSSD;
ssrt.mean = meanRT - ssdMeanWeight;






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Grand average of each method:
% The average of each of the estimates of ssrts
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ssrt.grand = nanmean([...
  nanmean(ssrt.integration);...
  ssrt.integrationWeighted;...
  ssrt.integrationSimple;...
  ssrt.mean]);








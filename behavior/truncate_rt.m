function [rt, outlierTrial] = truncate_rt(rt, minRT, maxRT, nSTD)

% [rt, outlierTrial] = truncate_rt(rt, minRT, maxRT, nSTD)
%
% Truncate a vector of input RTs, such that:
%
% minRT < rt < maxRT
% (rt - nSTD*std(rt)) < rt < (rt + nSTD*std(rt))

if nargin < 4, nSTD = 3; end
if nargin < 3, maxRT = 1200; end
if nargin < 2, minRT = 120; end

beforeMinRT = rt < minRT;
afterMaxRT = rt > maxRT;

rtMean = nanmean(rt);
rtStd = nanstd(rt);

extraRangeRT = rt < (rtMean - nSTD*rtStd) | rt > (rtMean + nSTD*rtStd);

% outlierTrial = find(beforeMinRT | afterMaxRT | extraRangeRT);
outlierTrial = beforeMinRT | afterMaxRT | extraRangeRT;
rt(outlierTrial) = nan;


function ssdList = ssd_session_adjust(ssdList)

% ssdList = ssd_session_adjust(ssdList)
%
% This function bins SSDs that are temporally near enough each other to be
% considered "the same" SSD. It does this by taking a weighted sum of the
% relative SSDs and returning that "adjusted" value as the "new" SSD.
%
% ssdList is any list of SSDs (from a dataset, e.g.)



% SSDs must be separated by this much time to not be adjusted to match other SSDs near this time
minTime = 17;

ssdArray = unique(ssdList(~isnan(ssdList)));
nSsd = arrayfun( @(x)(length(find(ssdList==x))), ssdArray);

belowCriterion = (diff(ssdArray) <= minTime)'; % Indices below minTime
firstIndexOfRun = find(diff([0,belowCriterion,0]) == 1); % First index of each run of consequtive SSDs below minTime
lastIndexOfRun = find(diff([0,belowCriterion,0]) == -1); % First index of each run of consequtive SSDs below minTime

ssdWeighted = nan(length(firstIndexOfRun), 1);
for iRun = 1 : length(firstIndexOfRun)
    iSsd = ssdArray(firstIndexOfRun(iRun):lastIndexOfRun(iRun));  % Which run of SSDs is this?
    iSsdN = nSsd(firstIndexOfRun(iRun):lastIndexOfRun(iRun));  % How many of each of the SSDs were there
    ssdWeighted(iRun) = round(sum(iSsd .* iSsdN / sum(iSsdN))); % Weighted avg of the SSDs
    ssdList(ismember(ssdList, iSsd)) = ssdWeighted(iRun);  % Replace old SSDs with new weighted avg SSDs
end














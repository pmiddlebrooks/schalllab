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



% Iterate until a weighted sum of SSDs are produced. Require the time
% between SSDs to be at least "minTime".
criterionMet = false;

belowCriterion = find(diff(ssdArray) <= minTime);
if isempty(belowCriterion)
    criterionMet = true;
end

while ~criterionMet
    
    nSsd = arrayfun( @(x)(length(find(ssdList==x))), ssdArray);
    
    
    % Are there any that have runs of more than 2 SSDs that are less than
    % criteria? If so, we need to give up and keep one.
    remove = 1+find(diff(belowCriterion) < 2);
    belowCriterion(remove) = [];
    
    ssdIndAltered = [belowCriterion; belowCriterion+1];
    ssdKeep = setxor(ssdIndAltered, 1:length(ssdArray));
    
    ssdWeighted = nan(length(belowCriterion), 1);
    for i = 1 : length(belowCriterion)
        
        ssdInd = [belowCriterion(i) belowCriterion(i)+1];
        ssdWeighted(i) = round(sum(ssdArray(ssdInd) .* nSsd(ssdInd) / sum(nSsd(ssdInd))));
        ssdList(ssdList == ssdArray(ssdInd(1)) | ssdList == ssdArray(ssdInd(2))) = ssdWeighted(i);
    end
    
    % The new ssdArray
    ssdArray = sort([ssdArray(ssdKeep); ssdWeighted]);
    
    % Check to see if there are any SSDs between "minTime"
    belowCriterion = find(diff(ssdArray) <= minTime);
    if isempty(belowCriterion)
        criterionMet = true;
    end
    
end
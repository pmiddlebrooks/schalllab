function rf = ccm_find_saccade_rf(Data)
%
% function rf = find_saccade_rf(Data);
%

% analyses windows relative to their aligned event
preTargWindow       = -299 : -1;
postCheckerWindow   = 50 : 150;
preSaccWindow       = -49 : 0;



rf = 'none'; % Default is no RF

leftFixAlign        = Data.targOn.colorCoh(1).goTarg.alignTime;
rightFixAlign       = Data.targOn.colorCoh(end).goTarg.alignTime;
leftCheckerAlign    = Data.checkerOn.colorCoh(1).goTarg.alignTime;
rightCheckerAlign   = Data.checkerOn.colorCoh(end).goTarg.alignTime;
leftSaccAlign       = Data.responseOnset.colorCoh(1).goTarg.alignTime;
rightSaccAlign      = Data.responseOnset.colorCoh(end).goTarg.alignTime;

leftFixRate         = nansum(Data.targOn.colorCoh(1).goTarg.raster(:,leftFixAlign + preTargWindow), 2)  .* 1000 ./ length(preTargWindow);
rightFixRate        = nansum(Data.targOn.colorCoh(end).goTarg.raster(:,rightFixAlign + preTargWindow), 2)  .* 1000 ./ length(preTargWindow);
leftCheckerRate     = nansum(Data.checkerOn.colorCoh(1).goTarg.raster(:,leftCheckerAlign + postCheckerWindow), 2)  .* 1000 ./ length(postCheckerWindow);
rightCheckerRate 	= nansum(Data.checkerOn.colorCoh(end).goTarg.raster(:,rightCheckerAlign + postCheckerWindow), 2)  .* 1000 ./ length(postCheckerWindow);
leftSaccRate        = nansum(Data.responseOnset.colorCoh(1).goTarg.raster(:,leftSaccAlign + preSaccWindow), 2)  .* 1000 ./ length(preSaccWindow);
rightSaccRate       = nansum(Data.responseOnset.colorCoh(end).goTarg.raster(:,rightSaccAlign + preSaccWindow), 2)  .* 1000 ./ length(preSaccWindow);


% Is there leftward presaccadic activity?
leftRF = false; % default is no
if sum([leftFixRate; leftSaccRate])
% 1) Is presacc > fixation activity?
[h , p] = ttest2(leftFixRate , leftSaccRate , .05);
if h && mean(leftSaccRate) > mean(leftFixRate)
    % 1) If so, is presacc > checker activity?
%     [h , p] = ttest2(leftCheckerRate , leftSaccRate , .05);
    if mean(leftSaccRate) > mean(leftCheckerRate)
        leftRF = true;
    end
end
end

% Is there rightward presaccadic activity?
rightRF = false; % default is no
if sum([rightFixRate; rightSaccRate])
% 1) Is presacc > fixation activity?
[h , p] = ttest2(rightFixRate , rightSaccRate , .05);
if h && mean(rightSaccRate) > mean(rightFixRate)
    % 1) If so, is presacc > checker activity?
%     [h , p] = ttest2(rightCheckerRate , rightSaccRate , .05);
    if mean(rightSaccRate) > mean(rightCheckerRate)
        rightRF = true;
    end
end
end

if leftRF && ~rightRF
    rf = 'left';
elseif ~leftRF && rightRF
    rf = 'right';
elseif leftRF && rightRF  % if both sides have saccade activity, figure out which has more and use it as the RF
    if mean(leftSaccRate) >= mean(rightSaccRate)
        rf = 'left';
    else
        rf = 'right';
    end
end



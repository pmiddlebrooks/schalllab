function [type, rf] = mem_classify_cell(Data)

rf = nan; % For now, return one receptive field. May need to alter in future to include vis and mov RFs separately
vis = false;
mov = false;
alpha = .05;
baseEpoch = -99 : 0;
visEpoch = mem_epoch_range('targOn', 'analyze');
movEpoch = mem_epoch_range('responseOnset', 'analyze');


baseRasR = Data.rightTarg.targOn.raster(:, baseEpoch + Data.rightTarg.targOn.alignTime);
baseRasL = Data.leftTarg.targOn.raster(:, baseEpoch + Data.leftTarg.targOn.alignTime);
rightBaseSpikeRate = nansum(baseRasR, 2) ./ length(baseEpoch);
leftBaseSpikeRate = nansum(baseRasL, 2) ./ length(baseEpoch);




% Visual cell?
rightRas = Data.rightTarg.targOn.raster(:, visEpoch + Data.rightTarg.targOn.alignTime);
rightVisSpikeRate = nansum(rightRas, 2) ./ length(visEpoch);
[h, p, ci, stats] = ttest2(rightBaseSpikeRate, rightVisSpikeRate);
if p < alpha && nanmean(rightBaseSpikeRate) < nanmean(rightVisSpikeRate)
    vis = true;
end

leftRas = Data.leftTarg.targOn.raster(:, visEpoch + Data.leftTarg.targOn.alignTime);
leftVisSpikeRate = nansum(leftRas, 2) ./ length(visEpoch);
[h, p, ci, stats] = ttest2(leftBaseSpikeRate, leftVisSpikeRate);
if p < alpha && nanmean(leftBaseSpikeRate) < nanmean(leftVisSpikeRate)
    vis = true;
end

% Visual receptive field
if vis
    if nanmean(leftVisSpikeRate) > nanmean(rightVisSpikeRate)
        visRF = 'left';
    else
        visRF = 'right';
    end
end

    




% Movement cell?
rightRas = Data.rightTarg.responseOnset.raster(:, movEpoch + Data.rightTarg.responseOnset.alignTime);
rightMovSpikeRate = nansum(rightRas, 2) ./ length(movEpoch);
[h, p, ci, stats] = ttest2(rightBaseSpikeRate, rightMovSpikeRate);
if p < alpha && nanmean(rightBaseSpikeRate) < nanmean(rightMovSpikeRate)
    mov = true;
end

leftRas = Data.leftTarg.responseOnset.raster(:, movEpoch + Data.leftTarg.responseOnset.alignTime);
leftMovSpikeRate = nansum(leftRas, 2) ./ length(movEpoch);
[h, p, ci, stats] = ttest2(leftBaseSpikeRate, leftMovSpikeRate);
if p < alpha && nanmean(leftBaseSpikeRate) < nanmean(leftMovSpikeRate)
    mov = true;
end

% Movement receptive field
if mov
    if nanmean(leftMovSpikeRate) > nanmean(rightMovSpikeRate)
        movRF = 'left';
    else
        movRF = 'right';
    end
end



if ~vis && ~mov
    type = nan;
elseif vis && ~mov
    type = 'visual';
    rf = visRF;
elseif ~vis && mov
    type = 'movement';
    rf = movRF;
elseif vis && mov
    type = 'visuomovement';
    % For vis/mov cell, figure out biggest response and use that as
    % receptive field
    rfArray = {visRF visRF movRF movRF};
    [~, ind] = max([nanmean(leftVisSpikeRate), nanmean(rightVisSpikeRate), nanmean(leftMovSpikeRate), nanmean(rightMovSpikeRate)]);
    rf = rfArray{ind};
end
%     end % kUnitIndex
end


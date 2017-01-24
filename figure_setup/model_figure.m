function [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = model_figure(nRow, nColumn, figureHandle, screenOrSave)

switch screenOrSave
    case 'screen'
        % For screen viewing, use:
        cmWidth = 40;
        cmHeight = 15;
    case 'save'
        % For printing/pdf, use:
        cmWidth = 22;
        cmHeight = 10;
    otherwise
        fprintf('You need to specify ''screen'' or ''save'' for variable ''screenOrSave''/n')
        return;
end


interSpaceX = 1.5;
interSpaceY = 1.5;
leftRightMargin = 1.5;
topBottomMargin = 1.5;


axisWidth = (cmWidth - 2*leftRightMargin - interSpaceX * (nColumn-1)) / nColumn;
axisHeight = (cmHeight - 2*topBottomMargin - interSpaceY * (nRow - 1)) / nRow;


set(0, 'units', 'centimeters')
scrsz = get(0,'ScreenSize');

screenH = scrsz(4);
screenW = scrsz(3);
w = screenW;
h = screenH;


sumFig = figure(figureHandle);
clf


xAxesPosition = nan(nRow, nColumn);
yAxesPosition = nan(nRow, nColumn);
for iRow = 1 : nRow
    for iCol = 1 : nColumn
        xAxesPosition(iRow, iCol) = leftRightMargin + axisWidth*(iCol-1) + interSpaceX*(iCol - 1);
        yAxesPosition(iRow, iCol) = cmHeight - (topBottomMargin + axisHeight*iRow + interSpaceY*(iRow - 1));
    end
end

set(sumFig, 'units', 'centimeters','position', [w/3 h/3 cmWidth cmHeight], 'paperunits', 'centimeters', 'paperposition', [0 0 cmWidth cmHeight], 'color', [1 1 1])



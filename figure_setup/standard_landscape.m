function [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle)

% Figure background color
bkgColor = [.7 .7 .7];
% bkgColor = [1 1 1];

% standard letter paper size
cmHeight = 21.6;
cmWidth = 27.9;

sumFig = figure(figureHandle);


interSpaceX = .8;
% interSpaceX = 1.3;
interSpaceY = 1.5;
interSpaceY = .8;
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



xAxesPosition = nan(nRow, nColumn);
yAxesPosition = nan(nRow, nColumn);
for iRow = 1 : nRow
    for iCol = 1 : nColumn
        xAxesPosition(iRow, iCol) = leftRightMargin + axisWidth*(iCol-1) + interSpaceX*(iCol - 1);
        yAxesPosition(iRow, iCol) = cmHeight - (topBottomMargin + axisHeight*iRow + interSpaceY*(iRow - 1));
    end
end

set(sumFig,'PaperOrientation','landscape');
set(sumFig, 'units', 'centimeters','position', [w/3 h/3 cmWidth cmHeight], 'paperunits', 'centimeters', 'paperposition', [0 0 cmWidth cmHeight], 'color', bkgColor)

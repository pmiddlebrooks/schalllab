function [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle, squareAxes)


if nargin < 4, squareAxes = 0; end


set(0, 'units', 'centimeters')
scrsz = get(0,'ScreenSize');

scaleFactor = .8;
cmHeight = scrsz(4) * scaleFactor;
if squareAxes
    cmWidth = cmHeight;
else
    cmWidth = scrsz(3) * scaleFactor;
end


interSpaceX = 1;
interSpaceY = interSpaceX;
leftRightMargin = 2;
topMargin = 2;
if squareAxes
    bottomMargin = 2;
else
    bottomMargin = 2;
end


axisWidth = (cmWidth - 2*leftRightMargin - interSpaceX * (nColumn-1)) / nColumn;
axisHeight = (cmHeight - topMargin - bottomMargin - interSpaceY * (nRow - 1)) / nRow;




sumFig = figure(figureHandle);
% clf


xAxesPosition = nan(nRow, nColumn);
yAxesPosition = nan(nRow, nColumn);
for iRow = 1 : nRow
    for iCol = 1 : nColumn
        xAxesPosition(iRow, iCol) = leftRightMargin + axisWidth*(iCol-1) + interSpaceX*(iCol - 1);
        yAxesPosition(iRow, iCol) = cmHeight - (topMargin + axisHeight*iRow + interSpaceY*(iRow - 1));
    end
end

set(sumFig, 'units', 'centimeters','position', [cmWidth*(1-scaleFactor)/2 cmHeight*(1-scaleFactor)/2 cmWidth cmHeight], 'paperunits', 'centimeters', 'paperposition', [0 0 cmWidth cmHeight], 'color', [1 1 1])



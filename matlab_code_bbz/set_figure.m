function hFig = set_figure(figSpecs,paperSpecs,varargin)
%
%
% INPUTS
% figSpecs     - cell array with elements specifying figure properties
%                * width
%                * height
%                * units: 'centimeters','inches','points','pixels':
% paperSpecs   - cell array with elements specifying paper properties
%                * type: 'USLetter','A4','A5',etc.
%                * orientation: 'landscape','portrait'
% otherSpecs   - cell array with elements specifying other aspects of the figure
%                (optional)
%                * fontName
%                * fontSize
% Current usage:
% Make new figure with specified height and width
% hFig = set_figure(height,width);
%
% EXAMPLES
%
% hFig = set_figure({18,10,'centimeters'},paperSpecs)
%

% Set up figure
hFig = figure;

% Use figure specs for proportional positioning and scaling of figure onscreen
[fWidth,fHeight] = deal(figSpecs{1:2});
screenSize = get(0,'ScreenSize');

switch figSpecs{3}
   case {'centimeters','inches','points'}
      fWidthPix = fWidth.*unitsratio('inches',figSpecs{3}).*get(0,'ScreenPixelsPerInch');
      fHeightPix = fHeight.*unitsratio('inches',figSpecs{3}).*get(0,'ScreenPixelsPerInch');
   case 'pixels'
      fWidthPix = fWidth;
      fHeightPix = fHeight;
end

if fWidthPix > screenSize(3) | fHeightPix > screenSize(4)
   screenWH = screenSize(3:4);
   pos = get(hFig,'Position');
   pos(3:4) = min(screenWH./[fWidth,fHeight]).*[fWidth,fHeight];
   set(hFig,'Position',[pos(3)/2,pos(4)/2,pos(3),pos(4)])
else
   leftPos = screenSize(3)/2-fWidthPix/2;
   bottomPos = screenSize(4)/2-fHeightPix/2;
   set(hFig,'Position',[leftPos,bottomPos,fWidthPix,fHeightPix])
end

% Center figure on paper
switch figSpecs{3}
   case {'centimeters','inches','points'}
      set(hFig,'PaperType',paperSpecs{1});
      set(hFig,'PaperUnits',figSpecs{3});
end

switch lower(paperSpecs{2})
   case 'landscape'
      orient landscape
   case 'portrait'
      orient portrait
   otherwise
      orient portrait
end

% Paper settings
pDim = num2cell(get(gcf,'PaperSize'));
[pWidth,pHeight] = deal(pDim{:});
xLeft = (pWidth-fWidth)/2; 
yTop = (pHeight-fHeight)/2;
set(hFig,'PaperPosition',[xLeft yTop fWidth fHeight])

% Other settings
if nargin > 2
   otherSpecs = varargin{1};
   set(hFig, 'DefaultAxesFontName', otherSpecs{1});
   set(hFig, 'DefaultTextFontName', otherSpecs{1});
   set(hFig, 'DefaultTextFontSize', otherSpecs{2}); % [pt]
   set(hFig, 'DefaultAxesFontSize', otherSpecs{2}); % [pt]
else
   set(hFig, 'DefaultAxesFontName', 'Helvetica');
   set(hFig, 'DefaultTextFontName', 'Helvetica');
   set(hFig, 'DefaultTextFontSize', 10); % [pt]
   set(hFig, 'DefaultAxesFontSize', 10); % [pt]
end

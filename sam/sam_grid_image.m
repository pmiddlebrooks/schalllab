function sam_grid_image(X,varargin)
% SAM_GRID_IMAGE Displays matrix as an image with grid lines
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% SAM_GRID_IMAGE(X);
% SAM_GRID_IMAGE(X,addValue);
% 
%
% addValue      - logical, whether or not to enter number in matrix
%
%
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Mon 05 Aug 2013 09:48:23 CDT by bram 
% $Modified: Tue 06 Aug 2013 11:10:30 CDT by bram

 
% CONTENTS 
% 1. FIRST LEVEL HEADER 
%    1.1 Second level header 
%        1.1.1 Third level header 


if islogical(X)
  
  cMap = flipud(colormap('gray'));
  if islogical(X)
    % Display matrix X, ones are black
    hI = imagesc(X,[0 1]); colormap(cMap);
  elseif isnumeric(X)
    % Display matrix X, max value in matrix is black
    hI = imagesc(X,[0 max(X(:))]); colormap(cMap);
  end
  axis image off

elseif isnumeric(X)
  
  
  heatmaptext(X,'FontColor','k','Colorbar',false,'precision',3);colormap(ones(64,3));
  axis image
  axis ij
end

% Add raster
  xLim = get(gca,'XLim');
  yLim = get(gca,'YLim');
  x = xLim(1):1:xLim(2);
  y = yLim(1):1:yLim(2);
  gridxy(x,y,'Color','k','LineWidth',1,'LineStyle','-');

% % Add value
% [xx,yy] = meshgrid(1:size(X,1),1:size(X,2));
% 
% 
% cData = get(hI,'CData');
% iW = find(cData > 0.5);
% iK = find(cData <= 0.5);
% 
% if islogical(X)
%   arrayfun(@(a,b,c) text(a,b,sprintf('%d',c),'Color','w'),xx(iW),yy(iW),X(iW),'Uni',0);
%   arrayfun(@(a,b,c) text(a,b,sprintf('%d',c),'Color','k'),xx(iK),yy(iK),X(iK),'Uni',0);
% elseif isnumeric(X)
%   arrayfun(@(a,b,c) text(a,b,sprintf('%.3f',c),'Color','w'),xx(iW),yy(iW),X(iW),'Uni',0);
%   arrayfun(@(a,b,c) text(a,b,sprintf('%.3f',c),'Color','k'),xx(iK),yy(iK),X(iK),'Uni',0);
% end
% 
% 
set(gca,'XLim',xLim + [-1 1],'YLim',yLim + [-1 1]);
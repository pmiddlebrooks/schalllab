function [xO,yO] = sam_quantile_average(xI,yI,varargin)
% SAM_QUANTILE_AVERAGE <Synopsis of what this function does> 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% SAM_QUANTILE_AVERAGE; 
% xI            - independent data (NxM double or Nx1 cell)
% yI            - dependent data (NxM double or Nx1 cell)
% p             - cumulative probabilities (Nx1 double)
%
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Mon 16 Sep 2013 12:16:32 CDT by bram 
% $Modified: Mon 16 Sep 2013 12:55:08 CDT by bram

 
% CONTENTS 

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% #.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 
% #.#. Input checking
% ========================================================================= 


nSample = size(xI,1);

if isnumeric(xI)
  
  if nargin == 3
    p = varargin{1};
    xO = quantile(xI,p);
    yO = quantile(yI,p);
  elseif nargin == 4
    i = varargin{1};
    p = varargin{2};
    xO = quantile(xI(i),p);
    yO = quantile(yI(i),p);
  end
elseif iscell(xI)
  
  if nargin == 3
    p = varargin{1};
    xO = cellfun(@(a,b) quantile(a,b),xI,p,'Uni',0);
    yO = cellfun(@(a,b) quantile(a,b),yI,p,'Uni',0);
  elseif nargin == 4
    i = varargin{1};
    p = varargin{2};
    p = repmat({p},nSample,1);
    xO = cellfun(@(a,b,c) quantile(a(c),b),xI,p,i,'Uni',0);
    yO = cellfun(@(a,b,c) quantile(a(c),b),yI,p,i,'Uni',0);
    
  elseif nargin == 5
    p = varargin{1};
    et = varargin{2};
    tWindow = varargin{3};
    
    % Align time series
    dt    = unique(diff(xI));
    t1    = xI(1);
    
  elseif nargin == 6
    i = varargin{1};
    p = varargin{2};
    p = repmat({p},nSample,1);
        
    et = varargin{3};
    tWindow = varargin{4};
    tWindow = repmat({tWindow},nSample,1);
    
    dt    = cellfun(@(a) unique(diff(a)),xI,'Uni',0);
    t1    = cellfun(@(a,b) a(b(1)),xI,i,'Uni',0);
    
    
    % Align time series
    [xA,yA] = cellfun(@(a,b,c,d,e,f) sam_align_to_event(a(b),c,d,e,f), ...
                                      yI,i,et,dt,t1,tWindow,'Uni',0);
  
    % Compute quantile time points
    xO = cellfun(@(a,b,c) quantile(a(~isnan(b)),c),xA,yA,p,'Uni',0);
    
    % Find corresponding y-values through linear interpolation
    yO = cellfun(@(a,b,c) interp1(a,b,c),xA,yA,xO,'Uni',0);
    
  end
  
  
  if iscell(xO)
    xO = cell2mat(xO);
  end
  
  if iscell(yO)
    yO = cell2mat(yO);
  end
  
  % #.#.#. Quantile averaging across samples
  % -----------------------------------------------------------------------
  
  xO = nanmean(xO,1);
  yO = nanmean(yO,1);
  
  
end
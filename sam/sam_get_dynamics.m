function varargout = sam_get_dynamics(dynType,t,z,iTr,iUnit,et,p,tWin)
% SAM_GET_DYNAMICS <Synopsis of what this function does> 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
% 
% SYNTAX 
% Let there be M inputs, N units, P time points, Q trials, then
% dynType   - type of dynamics
% t         - time points (1xP double)
% z         - dynamics (NxQXP double)
% iTr       - trial array (1xQ logical)
% iUnit     - (1x1 double or 1xQ double)
% et        - event times (1xQ double)
% p         - cumulative probabilities (column vector, double)
% tWin      - time window (1x2 double)
%
% SAM_GET_DYNAMICS; 
%  
% EXAMPLES 
% [tOut,zOut] = sam_get_dynamics('qaverage',t,z,iTr,iN,et,p,tWin);
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Fri 20 Sep 2013 09:40:06 CDT by bram 
% $Modified: Fri 20 Sep 2013 09:40:06 CDT by bram 

 
% CONTENTS 
% 1. FIRST LEVEL HEADER 
%    1.1 Second level header 
%        1.1.1 Third level header 

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS & SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1.#. Process inputs
% ========================================================================= 

if islogical(iTr)
  iTr = find(iTr);
end

if islogical(iUnit)
  iUnit = find(iUnit);
end

% If dynType equals qav-sample, find out how many trials to sample
if strncmp(dynType,'qav-sample',10)
  % Read the number of trials to sample
  nSample = strread(dynType,'qav-sample%d');
  
  % Correct nSample if it is lower than the actual number of trials
  if nSample > numel(iTr)
    nSample = numel(iTr);
  end
  % Remove number of samples from dynType string
  dynType = 'qav-sample';
end


% 1.#.#. Specify dynamic variables
% ------------------------------------------------------------------------- 

% Numbers
nTr         = numel(iTr);         % Trials
nTP         = numel(t);           % Time points
nUnit       = size(z,1);          % Units

% Convert variables to cell arrays
tCell       = repmat({t},nTr,1);

if isscalar(iUnit)
  zCell       = reshape(z(iUnit,iTr,:),nTr,nTP);
  zCell       = mat2cell(zCell,ones(nTr,1),nTP);
elseif isvector(iUnit)
  zCell       = mat2cell(z(:,iTr,:),size(z,1),ones(nTr,1),nTP);
  zCell       = cellfun(@(a) reshape(a,nUnit,nTP),zCell(:),'Uni',0);
  zCell       = cellfun(@(a,b) a(b,:),zCell,num2cell(iUnit(:)),'Uni',0);
end

etCell      = num2cell(et(:));

% #.#. Quantile averaging
% ========================================================================= 

switch lower(dynType)
  case {'qav','qav-sample'}
    
    % #.#.#. Identify time points on which dynamics are recorded
    % ---------------------------------------------------------------------
    iRecorded = cellfun(@(a) ~isnan(a),zCell,'Uni',0);
    
    % #.#.#. Compute quantile average
    % ---------------------------------------------------------------------
    [tQ,zQ] = sam_quantile_average(tCell,zCell,iRecorded,p,etCell,tWin);
end

switch lower(dynType)
  case 'qav-sample'
    
    dtCell = repmat({unique(diff(t))},nSample,1);
    t1Cell = repmat({t(1)},nSample,1);
    tWinCell = repmat({tWin},nSample,1);
    
    [tS,zS] = cellfun(@(a,b,c,d,e) sam_align_to_event(a,b,c,d,e), ...
                       zCell(1:nSample),etCell(1:nSample),dtCell, ...
                       t1Cell,tWinCell,'Uni',0);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% #. OUTPUTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

switch lower(dynType)
  case 'qav'
    varargout{1} = tQ;
    varargout{2} = zQ;
  case 'qav-sample'
    varargout{1} = tQ;
    varargout{2} = zQ;
    varargout{3} = tS;
    varargout{4} = zS;
end
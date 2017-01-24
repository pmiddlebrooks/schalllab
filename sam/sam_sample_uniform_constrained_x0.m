function X0 = sam_sample_uniform_constrained_x0(N,LB,UB,varargin)
% Sample uniformly distributed starting points with bounds, linear, and/or nonlinear constraints 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% SAM_SAMPLE_UNIFORM_CONSTRAINED_X0(N,LB,UB);
% SAM_SAMPLE_UNIFORM_CONSTRAINED_X0(N,LB,UB,A,b);
% SAM_SAMPLE_UNIFORM_CONSTRAINED_X0(N,LB,UB,A,b,nonLinCon,solverType);
%  
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Tue 01 Oct 2013 13:47:56 CDT by bram 
% $Modified: Tue 01 Oct 2013 13:47:56 CDT by bram 

% CONTENTS 
% 1. FIRST LEVEL HEADER 
%    1.1 Second level header 
%        1.1.1 Third level header 

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. FIRST LEVEL HEADER 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

MAX_N_REP     = 10000;

% 1.1. Second level header
% ========================================================================= 

 
% 1.1.1. Third level header
% ------------------------------------------------------------------------- 


if nargin == 3
  linConA   = [];
  linConB   = [];
  nonLinCon = [];
elseif nargin == 5
  linConA   = varargin{1};
  linConB   = varargin{2};
  nonLinCon = [];
elseif nargin == 7
  linConA   = varargin{1};
  linConB   = varargin{2};
  nonLinCon = varargin{3};
  solverType = varargin{4};
end

iRep = 1;
nX = length(LB);

X0 = [];

while iRep < MAX_N_REP + 1
 
  % Sample N points between bounds
  X = repmat(LB,N,1) + repmat(UB-LB,N,1).*rand(N,nX);

  % Convert X to a cell with column vectors
  XCell = mat2cell(X,ones(N,1),nX);
  XCell = cellfun(@(x) x(:),XCell,'Uni',0);

  % Check linear constraints, if any
  if ~isempty(linConA)
    iLinCon = cell2mat(cellfun(@(x) all(linConA*x < linConB),XCell,'Uni',0));
  end

  % Check nonlinear constraints, if any
  if ~isempty(nonLinCon)
    switch lower(solverType)
      case 'fminsearchcon'
        % Inequality constraints
        c = cellfun(@(x) nonLinCon(x),XCell,'Uni',0);

      case {'fmincon','ga'}
        % Inequality and equality constraints
        [c,~] = cellfun(@(x) nonLinCon(x),XCell,'Uni',0);
    end
    iNonLinCon = cell2mat(cellfun(@(x) all(x < 0),c,'Uni',0));
  end

  % Identify suitable starting parameters
  if isempty(linConA) & isempty(nonLinCon)
    iSuitable = 1:N;
  elseif isempty(linConA) & ~isempty(nonLinCon)
    iSuitable = find(iNonLinCon);
  elseif ~isempty(linConA) & isempty(nonLinCon)
    iSuitable = find(iLinCon)
  elseif ~isempty(linConA) & ~isempty(nonLinCon)
    iSuitable = find(iLinCon & iNonLinCon);
  end
  
  % Check if number of data points 
  if numel(iSuitable) == N
    X0 = X(iSuitable,:);
    break
  elseif numel(iSuitable) == 0
    X0 = X0;
%     error('Check constraints; none of the samples met all linear and nonlinear constraints');
  else
    X0 = [X0;X(iSuitable,:)];
  end
    
    if size(X0,1) >= N
      X0 = X0(1:N,:);
      break
    else
      iRep = iRep + 1;
    end
  end

  if iRep == MAX_N_REP
    warning('%d times sampling %d starting points did not yield sufficient starting points',MAX_N_REP,N);
    X0 = X0(1:N,:);
  end
end
function [tO,zO] = sam_align_to_event(zI,et,dt,t1,tWinO)
% SAM_ALIGN_TO_EVENT <Synopsis of what this function does> 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
%
% Each row in the matrix is aligned to the value in the corresponding row
% in |et|.
%
%
% Let there be N time series of M time points with time steps equal to dt.
%
% SYNTAX 
% SAM_ALIGN_TO_EVENT; 
%  
%
% zI        - input time series (NxM)
% et        - event times to align the time series to (Nx1 double)
% t1        - time of the first time point (scalar)
% dt        - time series time step size (scalar,integer)
% tWin      - time window (1x2 double), relative to aligned event
%
% OUTPUT
% tO        - time vector of aligned data
% zO        - output time series
%
%
% EXAMPLES 
% myfun = @(a,b,c,d,e,t) a.*(t-d).^2 + b.*(t-d) + c + e*randn(size(t));
% tI = 0:1:1000;
% zI = [myfun(-0.004,0,1000,500,0,tI);myfun(-0.006,0,1000,400,0,tI)];
% et = [500;400]; 
% dt = 1;
% t1 = x(1);
% tWin = [-100,100];
% [tO,zO] = sam_align_to_event(zI,et,dt,t1,tWin);
% figure;
% subplot(1,2,1);
% plot(tI,zI);
% xlabel('time (ms)');
% subplot(1,2,2);
% plot(tO,zO);
% xlabel('time from peak (ms)');
%
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Tue 27 Aug 2013 12:04:20 CDT by bram 
% $Modified: Tue 27 Aug 2013 12:04:20 CDT by bram 

 
% CONTENTS 
% 1. FIRST LEVEL HEADER 
%    1.1 Second level header 
%        1.1.1 Third level header 

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS AND SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1.1. Sanity checks on inputs
% ========================================================================= 

% 1.1.1. Verify input size and class
% -------------------------------------------------------------------------
if ~isequal(size(zI,1),size(et,1));
  error('zI and et should have the same number of rows.');
end

if ~isscalar(dt) || ~isscalar(t1)
  error('dt and t1 should be scalars.');
end

if ~isequal(numel(tWinO),2)
  error('tWin should have 2 elements.');
end
 
% 1.1.2. Round input times to have correct temporal resolution
% -------------------------------------------------------------------------
et        = round(et(:)/dt)*dt;                   % Event time(s)
tWinO     = round(tWinO(:)'/dt)*dt;               % Output time window
t1        = round(t1/dt)*dt;                      % Time at data point 1

% 1.2. Specify variables
% ========================================================================= 

% Numbers
nTS       = size(zI,1);                           % Time series
nTPI      = size(zI,2);                           % Input time points
nTPO      = numel(tWinO(1):dt:tWinO(2));          % Output time points

% Time vectors
tI        = t1:dt:(t1 + (nTPI - 1)*dt);           % Input
tO        = tWinO(1):dt:tWinO(2);                 % Output
tA        = repmat(tI,nTS,1)-repmat(et,1,nTPI);   % Aligned data

% 1.3. Pre-allocate empty arrays
% ========================================================================= 
zO        = nan(nTS,nTPO);                        % Output time series

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. ALIGN THE DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for i = 1:nTS                                     % Loop over time series
  
  zO(i,:) = [nan(1,numel(find(tO < tA(i,1)))), ...% Front NaN-padding
             zI(i,tA(i,:) >= tO(1) & tA(i,:) <= tO(end)), ... % Data
             nan(1,numel(find(tO > tA(i,end))))]; % Back NaN-padding
  
end
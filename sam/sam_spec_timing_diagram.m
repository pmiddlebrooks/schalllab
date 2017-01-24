function [T,u] = sam_spec_timing_diagram(ons,dur,v,eta,se,dt,tWindow) %#codegen
% SAM_SPEC_TIMING_DIAGRAM Specifies timing diagram
%  
% SAM_SPEC_TIMING_DIAGRAM constructs input matrices that specify the accumulators that
% are driven and the time when drive happens. The function has two output
% variables: t and u. The variable t specifies a vector of time points with
% temporal resolution dt. The variable u specifies onset and duration of
% the stimuli. 
%
%
% Let there be L events, M inputs, and Q time points
%
%
% Numbers:
% M events
% N inputs
% Q time points
%
% pad ons and dur with zeros to have equal number of elements per input
% 
%
% N number of events (with a given onset and durations)
% M number of inputs
% P number of time points (depends on dt and time window)
%
% SYNTAX
% [t,u] = SAM_SPEC_TIMING_DIAGRAM(ons,dur,v,se,dt,tWindow);
% ons           - stimulus onsets (LxM double)
% dur           - stimulus durations (LxM double)
% v             - mean accumulation rates (Mx1 double)
%                 if empty, defaults to ones(M,1)
% eta           - between-trial variability in accumulation rate (Mx1 double)
%                 if empty, defaults to zeros(M,1)
% se            - extrinsic noise magnitudes (MxM double)
%                 if empty, defaults to zeros(M,1)
% dt            - time step size, in ms (scalar)
%                 if empty, defaults to 1
% tWindow       - specifies time at start and end of u (1x2 double)
%                 if empty, defaults to [0,2000]
%
% t             - time vector (1xQ double)
% u             - timing diagram (MxQ double)
%
% EXAMPLES
%
% Example 1:
% ons           = [0 0 0];
% dur           = [2000 2000 0];
% v             = [3,3,5]';
% eta           = [0,0,0]';
% se            = diag([0.5 0.5 0.5]);
% dt            = 10;
% tWindow       = [-500 2500];
%
% [t,u] = SAM_SPEC_TIMING_DIAGRAM(ons,dur,v,eta,se,dt,tWindow);
% figure;stairs(t,u');
% xlabel('time (ms)')
% ylabel('signal strength (a.u.)');
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Wed 24 Jul 2013 12:53:20 CDT by bram 
% $Modified: Fri 23 Aug 2013 11:19:24 CDT by bram

% CONTENTS
% 1.HANDLING OF INPUTS AND DEFINING VARIABLES
%   1.1. Handle inputs
%   1.2. Define variables based on inputs
%   1.3. Pre-allocate matrices
% 2.DEFINE STIMULUS AND INPUT MATRICES
%   2.1. Code onset (1) and offset (-1) in stimulus and input matrices
%   2.2. Specify duration of stimuli and inputs
%   2.3. Block inputs, if requested
%   2.4. Make sparse

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.1. Process inputs
% =========================================================================

if isempty(v)
  v = ones(size(ons,2),1);
end

if isempty(eta)
  eta = ones(size(ons,2),1);
end

if isempty(se)
  se = zeros(size(ons,2),1);
end

if isempty(dt)
  dt = 1;
elseif numel(dt) ~= 1
  error('dt should be a scalar.');
end

if isempty(tWindow)
  tWindow = [0 2000];
elseif ~isequal(size(tWindow),[1 2])
  error('tWindow should be a 1x2 double');
end

% 1.2. Specify variables based on inputs
% =========================================================================

% Time vector
T       = round((tWindow(1):dt:tWindow(2))'./dt)*dt;

% Numbers
l       = size(ons,1);                  % Events
m       = size(ons,2);                  % Inputs
q       = length(T);                    % Time points

% Event onsets
tOn    = round(ons/dt)*dt;

% Event offsets
tOff   = round((ons + dur)/dt)*dt;

% 1.3. Pre-allocate matrices
% =========================================================================
u           = zeros(m,q);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. SPECIFY STIMULUS AND INPUT MATRICES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.1. Specify event onset (1) and offset (-1) in timing diagrams
% =========================================================================

for h = 1:l % Event index
  for i = 1:m % Input index

    % Event onsets in stimulus timing diagram
    u(i,T == tOn(h,i)) = u(i,T == tOn(h,i)) + 1;

    % Event offsets in stimulus timing diagram
    u(i,T == tOff(h,i)) = u(i,T == tOff(h,i)) -1;

  end
end
  
% 2.2. Specify duration of stimuli and inputs
% =========================================================================
u   = cumsum(u,2);

% 2.3. Specify strength of inputs
% =========================================================================
u = diag((randn(m,1) .* eta) + v) * u;

% 2.4. Add extrinsic noise of weight se
% =========================================================================
u = se*randn(m,q) + u;

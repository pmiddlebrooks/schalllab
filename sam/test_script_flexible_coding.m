%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. SPECIFY TIMING MATRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4.1.1. Compute timing diagram
% -----------------------------------------------------------------
                                  % OUTPUT
[tStm, ...                        % - Time
 uStm] ...                        % - Strength of stimulus (t)
= sam_spec_timing_diagram ...     % FUNCTION
 ...                              % INPUT
(stimOns{iCnd,iTrType}(:)', ...   % - Stimulus onset time
 stimDur{iCnd,iTrType}(:)', ...   % - Stimulus duration
 [], ...                          % - Strength (default = 1);
 0, ...                           % - Magnitude of extrinsic noise
 dt, ...                          % - Time step
 timeWindow);                     % - Time window

% 4.2. Timing diagram of model inputs
% =====================================================================

% Onset and duration of the accumulation process
thisAccumOns = accumOns{iCnd,iTrType}(:)';
thisAccumDur = stimDur{iCnd,iTrType}(:)';

% 4.2.1. Adjust duration of the STOP process, if needed
% ---------------------------------------------------------------------
if iTrType > 1
switch durationSTOP
case 'trial'
  % STOP accumulation process lasts the entire trial
  thisAccumDur(iSTOP) = timeWindow(2) - thisAccumOns(iSTOP);
end
end
% 4.2.2. Compute timing diagram
% ---------------------------------------------------------------------

                                  % OUTPUT
[t, ...                               % - Time
u] ...                               % - Strength of model input
= sam_spec_timing_diagram ...        % FUNCTION
...                                  % INPUT
(thisAccumOns, ...                    % - Accumulation onset time
thisAccumDur, ...                    % - Stimulus duration
V{iCnd,iTrType}, ...                 % - Input strength
SE{iCnd,iTrType}, ...                % - Magnitude of extrinsic noise
dt, ...                              % - Time step
timeWindow);                         % - Time window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. SIMULATE TRIALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-defining the following variables reduces the overhead in the trial
% simulation MEX functions called below
n   = size(A,1);  % Number of units
m   = size(C,1);  % Number of inputs to units
p   = size(u,2);  % Number of time points
t1  = 1;          % First time point

% Pre-allocate response time and response arrays
rt  = inf(sum(N),nSim);

% Pre-allocate response and dynamics array
z = nan(sum(N),nSim,numel(t));

% Reduce overhea in parfor loops by getting rid of indexing (and
% using feval)
thisC  = C{iCnd};
thisSI = diag(SI{iCnd,iTrType});
thisZC = ZC{iCnd};
thisT  = t(:)';

% Reduce overhead in MEX functions by preallocating trial-dependent
% variables
thisRt    = inf(n,1);       % Response time
thisResp  = false(n,1);     % Response (i.e. whether a unit has reached zc)
thisZ     = nan(n,p);       % Activation

% Loop over trials
parfor iTr = 1:nSim
                                  % OUTPUT                          
  [rt(:,iTr), ...                 % - Response time
   ~, ...                         % - Responses
   z(:,iTr,:)] ...                % - Dynamics
   =  feval ...                   % FUNCTION
   ...                            % INPUT
  (trialSimFun, ...               % - Function handle
   u, ...                         % - Timing diagram of model inputs
   A, ...                         % - Endogenous connectivity matrix
   B, ...                         % - Extrinsic modulation matrix
   thisC, ...                     % - Exogenous connectivity matrix
   D, ...                         % - Intrinsic modulation matrix
   thisSI, ...                    % - Intrinsic noise matrix
   0 + (Z0-0).*rand(n,1), ...     % - Starting point matrix (uniformly distributed)
   thisZC, ...                    % - Threshold matrix
   ZLB, ...                       % - Activation lower bound matrix
   dt, ...                        % - Time step
   tau, ...                       % - Time scale
   thisT, ...                     % - Time points
   terminate, ...                 % - Termination matrix
   blockInput, ...                % - Blocked input matrix
   latInhib, ...                  % - Lateral inhibition matrix
   n, ...                         % - Number of units
   m, ...                         % - Number of inputs
   p, ...                         % - Number of time points
   t1, ...                        % - First time index
   thisRt, ...                    % - Array for current trial's RT
   thisResp, ...                  % - Array for current trial's response
   thisZ);                        % - Array for current trial's dynamics
end
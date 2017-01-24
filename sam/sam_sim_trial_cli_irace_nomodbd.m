function [rt,resp,z] = sam_sim_trial_cli_irace_nomodbd(u,A,~,C,~,SI, ...
                                                         Z0,ZC,ZLB,dt, ...
                                                         tau,T, ...
                                                         terminate,~,~, ...
                                                         n,~,p,t, ...
                                                         rt,resp, ...
                                                         z)
% Simulate trials: C as lateral inhibition, I as a race, no extr. and intr. modulation
% 
% DESCRIPTION 
% SAM trial simulation function, modeling choice as lateral inhibition, 
% inhibition as a race, and excluding extrinisic (B) or intrinisic (D) 
% modulation of connectivity.
% 
% SYNTAX 
% Let there be M inputs, N units, P time points, then
% u           - inputs to the accumulators (MxP double)
% A           - endogenous connectivity matrix (NxN double)
% B           - extrinsic modulation matrix (NxNxM double)
% C           - exogenous connectivity matrix (NxM double)
% D           - intrinsic modulation matrix (NxNxN double)
% SI          - intrinsic noise strength (MxM double)
% Z0          - starting value of activation (Nx1 double)
% ZC          - threshold on activation (Nx1 double)
% ZLB         - lower bound on activation (Nx1 double)
% dt          - time step (1x1 double)
% tau         - time scale (1x1 double)
% T           - time points (1xP double)
% terminate   - matrix indicating which units can terminate accumulation of
%               activation when they reach threshold (Nx1 logical)
% blockInput  - matrix indicating which units block which inputs when they
%               reach threshold (Nx1 logical)
% latInhib    - matrix indicating which elements in A remain 0 as long as
%               unit n (indexed by the columns of A) has not reached 
%               threshold (indexed by resp)
% n           - number of units (1x1 double)
% m           - number of inputs (1x1 double)
% p           - number of time points (1x1 double)
% t           - first time point (1x1 double)
% rt          - array for logging response time (Nx1 double)
% resp        - array for logging responses (Nx1 logical)
% z           - array for logging dynamics (NxP double)
%
% rt          - response times (Nx1 double)
% resp        - responses, inid (Nx1 logical)
% z           - dynamics (NxP double)
%
% [rt,resp,z] = SAM_SIM_TRIAL_CLI_IRACE_NOMODBD(u,A,~,C,~,SI,Z0,ZC, ...
%                                                 ZLB,dt,tau,T, ...
%                                                 terminate,~,~,n,~,p, ...
%                                                 t,rt,resp,z);
%
% EXAMPLES 
%
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Wed 24 Jul 2013 12:14:48 CDT by bram 
% $Modified: Wed 25 Sep 2013 11:00:45 CDT by bram

% Set starting values of z(t)
z(:,1)  = Z0;     

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. STOCHASTIC ACCUMULATION PROCESS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

while t < p - 1
  
  % Endogenous connectivity at time t (note that A is a function of t
  % because lateral inhibition kicks in once a unit has reached its
  % threshold)
  % -----------------------------------------------------------------------
%   At = A;
  
%   % Extrinsic modulation at time t
%   % -----------------------------------------------------------------------
%   Bt          = zeros(m,m);
%   for i = 1:m
%     Bt        = Bt + u(i,t)*B(:,:,i);
%   end

%   % Intrinsic modulation at time t
%   % -----------------------------------------------------------------------
%   Dt          = zeros(n,n);
%   for j = 1:n
%     Dt        = Dt + z(j,t)*D(:,:,j);
%   end

%   % Inhibition mechanism 1: block input(s), if any
%   % -----------------------------------------------------------------------
%   u(any(blockInput(:,resp),2),t) = 0;
  
%   % Inhibition mechanism 2: lateral inhibition
%   % -----------------------------------------------------------------------
% %   At(latInhib(:,~resp)) = 0; % This produces wrong indexing!
%   At(latInhib & logical(repmat(~resp',n,1))) = 0;
  
%   % Change in activation from time t to t + 1
%   % -----------------------------------------------------------------------
%   dzdt        = (At + Bt + Dt)  * z(:,t)      * dt/tau + ...  % 
%                 C               * u(:,t)      * dt/tau + ...  % Inputs
%                 SI             * randn(n,1)  * sqrt(dt/tau); % Noise (in)
              
  dzdt        = A   * z(:,t)      * dt/tau + ...   % Endogenous connectivity
                C   * u(:,t)      * dt/tau + ...   % Inputs
                SI  * randn(n,1)  * sqrt(dt/tau);  % Noise (in)
              
  % Log new activation level
  % -----------------------------------------------------------------------
  z(:,t+1)    = z(:,t) + dzdt;

  % Rectify activation if below zLB
  % -----------------------------------------------------------------------
  z(z(:,t+1) < ZLB,t+1) = ZLB(z(:,t+1) < ZLB);

  % Identify units that crossed threshold
  % -----------------------------------------------------------------------
  resp(z(:,t+1) > ZC) = true;
  
  % Determine time of crossing threshold
  % -----------------------------------------------------------------------
  rt(resp & isinf(rt)) = T(t + 1);
  
  % Break accumulation if termination criterion has been met
  % -----------------------------------------------------------------------
  if any(terminate(resp))
    break
  end
  
  % Update time
  % -----------------------------------------------------------------------
  t = t + 1;
  
end
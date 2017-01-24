function prd = sam_sim_expt(simGoal,X,SAM)
% Simulates response times and model dynamics for go and stop trials of the
% stop-signal task
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX
% prd = SAM_SIM_EXPT('optimize',X,SAM); 
% prd = SAM_SIM_EXPT('explore',X,SAM); 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Sat 21 Sep 2013 12:54:52 CDT by bram 
% $Modified: Mon 23 Sep 2013 20:54:22 CDT by bramzandbelt

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS AND SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1.1. Process inputs
% ========================================================================= 

nTrialCat     = size(SAM.optim.obs,1);

% 1.1.1. Experiment variables
% -------------------------------------------------------------------------------------------------------------------------

N             = SAM.expt.nRsp;
trialDur      = SAM.expt.trialDur;

nRsp          = SAM.expt.nRsp;
% nStm          = SAM.expt.nStm;

% 1.1.2. Model variables
% -------------------------------------------------------------------------------------------------------------------------

terminate     = SAM.model.mat.terminate;
blockInput    = SAM.model.mat.interClassBlockInp;
latInhib      = SAM.model.mat.interClassLatInhib;


dt            = SAM.model.accum.dt;
timeWindow    = SAM.model.accum.window;
zLB           = SAM.model.accum.zLB;
tau           = SAM.model.accum.tau;
accumTWindow  = SAM.model.accum.window;
durSTOP       = SAM.model.accum.durSTOP;

% We can control the random variability in starting point, non-decision time, and accumulation rate
if SAM.model.accum.randomZ0
  randomZ0Factor = 1;
else
  randomZ0Factor = realmin;
end

if SAM.model.accum.randomT0
  randomT0Factor = 1;
else
  randomT0Factor = realmin;
end

T       = (accumTWindow(1):dt:accumTWindow(2))';
p   = length(T);
t1  = 1;                    % First time point

% 1.1.3. Simulation variables
% -------------------------------------------------------------------------------------------------------------------------
nSim          = SAM.sim.n;
trialSimFun   = SAM.sim.fun.trial;
alignTWindow  = SAM.sim.tWindow;
rngSeedStage  = SAM.sim.rng.stage;
rngSeedId     = SAM.sim.rng.id;
qntls         = SAM.optim.cost.stat.cumProb;
minBinSize    = SAM.optim.cost.stat.minBinSize;
    
% 1.1.4. Optimization variables
% -------------------------------------------------------------------------------------------------------------------------
prd           = SAM.optim.prd;

% 1.2. Specify static variables
% =========================================================================

switch simGoal
  case 'explore'
    % Cumulative probabilities for quantile averaging of model dynamics
    CUM_PROB = 0:0.01:1;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. SEED THE RANDOM NUMBER GENERATOR (OPTIONAL)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

switch lower(rngSeedStage)
  case 'sam_sim_expt'

    % Note: MEX functions stay in memory until they are cleared.
    % Seeding of the random number generator should be accompanied by 
    % clearing MEX functions.

    clear(char(trialSimFun));
    rng(rngSeedId);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3. SIMULATE EXPERIMEMT
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Loop over trial categories
for iTrialCat = 1:nTrialCat
  
  trialCat    = SAM.optim.obs.trialCat{iTrialCat};
  stmOns      = SAM.optim.obs.onset{iTrialCat};
  stmDur      = SAM.optim.obs.duration{iTrialCat};
  
  iTargetGO       = find(SAM.optim.modelMat.iTarget{iTrialCat}{1});
  iNonTargetGO    = find(SAM.optim.modelMat.iNonTarget{iTrialCat}{1});
  iGO             = sort([iTargetGO(:);iNonTargetGO(:)]);
  iTargetSTOP     = find(SAM.optim.modelMat.iTarget{iTrialCat}{2});
  iNonTargetSTOP  = find(SAM.optim.modelMat.iNonTarget{iTrialCat}{2});
  iSTOP           = sort([iTargetSTOP(:);iNonTargetSTOP(:)]);
  
  % Pre-allocate response time and response arrays
  rt  = inf(sum(N),nSim);
  
  switch simGoal
    case 'explore'
      z = nan(sum(nRsp),nSim,p);
%       uLog = nan(sum(nStm),nSim,p);
  end
  
  % 3.1. Decode parameter vector
  % =====================================================================
  [endoConn, ...
   extrMod, ...
   exoConn, ...
   intrMod, ...
   V, ...
   ETA, ...
   SE, ...
   SI, ...
   Z0, ...
   ZC, ...
   T0] ...
   ...
   = sam_decode_x( ...
   ...
   SAM, ...
   X, ...
   iTrialCat);
  
  n   = size(endoConn,1);     % Number of units
  m   = size(exoConn,2);      % Number of inputs to units
  ZLB = zLB*ones(n,1);
  
  % 3.2. Simulate trials
  % =====================================================================
  
  % Make sure that these variables are row vectors
  stmOns  = stmOns(:)';
  stmDur  = stmDur(:)';
  T0      = T0(:)';
  
  if ~isempty(regexp(trialCat,'stopTrial.*', 'once'))
    switch lower(durSTOP)
      case 'trial'
        accumDurFactor = zeros(1,m);
        accumDurFactor(iSTOP) = 1;
      otherwise
        accumDurFactor = zeros(1,m);
    end
  else
    accumDurFactor = zeros(1,m);
  end

  for iTr = 1:nSim
    
    % 3.2.1. Timing diagram of model inputs
    % -----------------------------------------------------------------------------------------------------------------

    t0Var = randomT0Factor.*rand(1,m);
    accumOns = stmOns + T0 - T0.*t0Var;
    accumDur = stmDur - accumDurFactor.*(trialDur-stmDur-accumOns);
    
    [t, ...                                               % - Time
     u] ...                                               % - Strength of model input
    = sam_spec_timing_diagram_mex ...                     % FUNCTION
    ...                                                   % INPUT
    (accumOns, ...                                        % - Accumulation onset time
     accumDur, ...                                        % - Stimulus duration
     V, ...                                               % - Input strength
     ETA, ...                                             % - Trial-to-trial variability in input strength
     SE, ...                                              % - Magnitude of extrinsic noise
     dt, ...                                              % - Time step
     timeWindow);                                         % - Time window
    

    t = t(:)';

    % 3.2.2. Simulate trials
    % -----------------------------------------------------------------------------------------------------------------

    switch simGoal
      case 'optimize'
        rt(:,iTr) ...                  % - Response time
        =  feval ...                   % FUNCTION
        ...                            % INPUT
        (trialSimFun, ...               % - Function handle
        u, ...                         % - Timing diagram of model inputs
        endoConn, ...                  % - Endogenous connectivity matrix
        extrMod, ...                   % - Extrinsic modulation matrix
        exoConn, ...                   % - Exogenous connectivity matrix
        intrMod, ...                   % - Intrinsic modulation matrix
        SI, ...                        % - Intrinsic noise matrix
        Z0 - (Z0-zLB).*randomZ0Factor.*rand(n,1), ...                        % - Starting point matrix
        ZC, ...                        % - Threshold matrix
        ZLB, ...                       % - Activation lower bound matrix
        dt, ...                        % - Time step
        tau, ...                       % - Time scale
        t, ...                         % - Time points
        terminate, ...                 % - Termination matrix
        blockInput, ...                % - Blocked input matrix
        latInhib, ...                  % - Lateral inhibition matrix
        n, ...                         % - Number of units
        m, ...                         % - Number of inputs
        p, ...                         % - Number of time points
        t1, ...                        % - First time index
        inf(n,1), ...                  % - Array for current trial's RT
        false(n,1), ...                % - Array for current trial's response
        nan(n,p));                     % - Array for current trial's dynamics)
      
      case 'explore'
        
        [rt(:,iTr), ...                  % - Response time
         ~, ...                         % - Responses
         z(:,iTr,:)] ...                % - Dynamics  
        ...
        =  feval ...                   % FUNCTION
        ...                            % INPUT
        (trialSimFun, ...               % - Function handle
        u, ...                         % - Timing diagram of model inputs
        endoConn, ...                  % - Endogenous connectivity matrix
        extrMod, ...                   % - Extrinsic modulation matrix
        exoConn, ...                   % - Exogenous connectivity matrix
        intrMod, ...                   % - Intrinsic modulation matrix
        SI, ...                        % - Intrinsic noise matrix
        Z0 - (Z0-zLB).*randomZ0Factor.*rand(n,1), ...                        % - Starting point matrix
        ZC, ...                        % - Threshold matrix
        ZLB, ...                       % - Activation lower bound matrix
        dt, ...                        % - Time step
        tau, ...                       % - Time scale
        t, ...                         % - Time points
        terminate, ...                 % - Termination matrix
        blockInput, ...                % - Blocked input matrix
        latInhib, ...                  % - Lateral inhibition matrix
        n, ...                         % - Number of units
        m, ...                         % - Number of inputs
        p, ...                         % - Number of time points
        t1, ...                        % - First time index
        inf(n,1), ...                  % - Array for current trial's RT
        false(n,1), ...                % - Array for current trial's response
        nan(n,p));                     % - Array for current trial's dynamics)
      
%         uLog(:,iTr,:) = u;
    end
  end
  
  % 3.4. Classify trials
  % =====================================================================
  
  if ~isempty(regexp(trialCat,'goTrial.*', 'once'))
    
    % 4.4.1. Go correct trial: only one RT, produced by target GO unit
    % -------------------------------------------------------------------
    % Criteria (all should be met):
    % - Target GO unit produced an RT
    % - Target GO unit is the only unit having produced an RT
    iCorr = rt(iTargetGO,:) < Inf & sum(rt < Inf) == 1;
    
    % 4.4.2. Go commission error trial: any RT produced by a non-target GO unit
    % -------------------------------------------------------------------
    % Criteria (all should be met):
    % - At least one non-target GO unit has produced an RT
    iError = any(rt(iNonTargetGO,:) < Inf,1);

    % 4.4.3. Go omission error trial: no RT
    % -------------------------------------------------------------------
    % Criteria (all should be met):
    % - No unit has produced an RT
%     iOmit = sum(rt < Inf) == 0;

  elseif ~isempty(regexp(trialCat,'stopTrial.*', 'once'))
    
    if all(SAM.model.mat.endoConn.nonSelfOther == 0) % GO and STOP race
      % Stop success trial
      % -------------------------------------------------------------------
      % Criteria (*any* should be met)
      % - STOP unit produces a shorter RT than any other GO units
      % - None of the units produced an RT
      iCorr = rt(iTargetSTOP,:) < min(rt(iGO,:),[],1) | all(rt == Inf);
      
      % Stop failure trial
      % -------------------------------------------------------------------
      % Criteria
      % - The fastest GO unit has a shorter RT than any STOP unit
      iError = min(rt(iGO,:),[],1) < rt(iSTOP,:);
      
      % Note: no distinction made between correct and error choice in
      % Go task
      
%       iOmit = [];
    else
      % Stop success trial
      % -------------------------------------------------------------------
      % Criteria (all should be met):
      % - None of the GO units produced an RT
      iCorr = all(rt(iGO,:) == Inf);
      
      % Note: no distinction made between trials in which the STOP unit
      % did and did not produce and RT. We may be more stringent,
      % requiring that STOP finishes, or we could look at how
      % computation of SSRT is influenced by including trials in which
      % STOP did not finish versus when these trials are not included.
      
      % Stop failure trial
      % -------------------------------------------------------------------
      % Criteria:
      % - At least one GO unit has produced an RT
      iError = any(rt(iGO,:) < Inf);

      % Note: no distinction made between correct and error choice in
      % Go task

      % Note: no distinction made between trials in which the STOP unit
      % did and did not produce an RT (i.e. STOP units cannot terminate
      % the trial under 'blocked input' and 'lateral inhibition'
      % inhibition mechanisms)
      
%       iOmit = [];
    end
    
  end
  
  % 3.5. Log model predictions
  % =====================================================================

  % 3.5.1. Trial numbers
  % ---------------------------------------------------------------------
  
  % Compute
  nCorr = numel(find(iCorr));
  nError = numel(find(iError));
%   nOmit = numel(find(iOmit));
  
  % Log
  prd.nTotal(iTrialCat) = nCorr + nError;
  prd.nCorr(iTrialCat) = nCorr;
  prd.nError(iTrialCat) = nError;
%   prd.nOmit(iTrialCat) = nOmit;
  
  
  % 3.5.2. Trial probabilities
  % ---------------------------------------------------------------------
  
  % Compute
  pCorr = nCorr./nSim;
  pError = nError./nSim;
%   pOmit = nOmit./nSim;
  
  % Log
  prd.pTotal(iTrialCat) = pCorr + pError;
  prd.pCorr(iTrialCat) = pCorr;
  prd.pError(iTrialCat) = pError;
%   prd.pOmit(iTrialCat) = pOmit;
  
  % 3.5.3. Trial response times
  % ---------------------------------------------------------------------
  if ~isempty(regexp(trialCat,'goTrial.*', 'once'))
    if nCorr > 0
      rtCorr = sort(rt(iTargetGO,iCorr)) - stmOns(iTargetGO);
      
      prd.rtCorr{iTrialCat} = rtCorr;
    end
    
    if nError > 0
      rtError = sort(min(rt(iNonTargetGO,iError),[],1)) - stmOns(iTargetGO);
      
      % Note: stimOns for iTargetGO and iNonTargetGO are the same; I use iTargetGO instead
      % of iNonTargetGO because iTargetGO is always a scalar, iNonTargetGO not.
      
      prd.rtError{iTrialCat} = rtError;
    end
    
  elseif ~isempty(regexp(trialCat,'stopTrial.*', 'once'))
    if nCorr > 0
      rtCorr = sort(rt(iSTOP,iCorr)) - stmOns(iSTOP);
      
      prd.rtCorr{iTrialCat} = rtCorr;
    end
    
    if nError > 0
      rtError = sort(min(rt(iGO,iError),[],1)) - stmOns(iTargetGO);
      
      prd.rtError{iTrialCat} = rtError;
    end
    
  end
  
  % 3.5.4. Model matrices used for simulating this trial category
  % -----------------------------------------------------------------------------------------------------------------------
  if nCorr > 0
    
    [prd.rtQCorr{iTrialCat}, ...
     prd.pDefectiveCorr{iTrialCat}, ...
     prd.fCorr{iTrialCat}, ...
     prd.pMassCorr{iTrialCat}] = ...
     ...
     sam_bin_data(...
     prd.rtCorr{iTrialCat}, ...
     prd.pCorr(iTrialCat), ...
     prd.nCorr(iTrialCat), ...
     qntls, ...
     minBinSize);
    
  end
  
  if nError > 0
    
    [prd.rtQError{iTrialCat}, ...
     prd.pDefectiveError{iTrialCat}, ...
     prd.fError{iTrialCat}, ...
     prd.pMassError{iTrialCat}] = ...
     ...
     sam_bin_data(...
     prd.rtError{iTrialCat}, ...
     prd.pError(iTrialCat), ...
     prd.nError(iTrialCat), ...
     qntls, ...
     minBinSize);
    
  end
  
  % 3.5.5. Model matrices used for simulating this trial category
  % -----------------------------------------------------------------------------------------------------------------------
  
  prd.modelMat{iTrialCat}.endoConn  = endoConn;
  prd.modelMat{iTrialCat}.extrMod   = extrMod;
  prd.modelMat{iTrialCat}.exoConn   = exoConn;
  prd.modelMat{iTrialCat}.intrMod   = intrMod;
  prd.modelMat{iTrialCat}.Z0        = Z0;
  prd.modelMat{iTrialCat}.ZC        = ZC;
  prd.modelMat{iTrialCat}.V         = V;
  prd.modelMat{iTrialCat}.ETA       = ETA;
  prd.modelMat{iTrialCat}.T0        = T0;
  prd.modelMat{iTrialCat}.SE        = SE;
  prd.modelMat{iTrialCat}.SI        = SI;
  prd.modelMat{iTrialCat}.ZLB       = ZLB;
  
  
  % 3.5.6. Model inputts
  % -----------------------------------------------------------------------------------------------------------------------
  
  % Implement this
  
  
  % 3.5.7. Model dynamics
  % -----------------------------------------------------------------------------------------------------------------------
  switch simGoal
    case 'explore'
      if ~isempty(regexp(trialCat,'goTrial.*', 'once'))
        
        % Go correct quantile averaged dynamics
        if nCorr >= 1

          % Randomly sample one non-target GO unit
          thisINonTargetGO = randsample([iNonTargetGO(:);iNonTargetGO(:)],1);  % N.B. Repeating iNonTargetGO prevents that iNonTargetGO is a scalar, under which randsample samples from 1:iNonTargetGO

          % Get event times of go-signals and correct responses  
          etGo = repmat(stmOns(iTargetGO),nCorr,1);
          etResp = rt(iTargetGO,iCorr);

          % Get quantile averaged dynamics of targetGO aligned on go-signal
          [prd.dyn{iTrialCat}.corr.goStim.targetGO.qX, ...
           prd.dyn{iTrialCat}.corr.goStim.targetGO.qY, ...
           prd.dyn{iTrialCat}.corr.goStim.targetGO.sX, ...
           prd.dyn{iTrialCat}.corr.goStim.targetGO.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iCorr,iTargetGO,etGo,CUM_PROB,alignTWindow.go);

          % Get quantile averaged dynamics of nonTargetGO aligned on go-signal
          [prd.dyn{iTrialCat}.corr.goStim.nonTargetGO.qX, ...
           prd.dyn{iTrialCat}.corr.goStim.nonTargetGO.qY, ...
           prd.dyn{iTrialCat}.corr.goStim.nonTargetGO.sX, ...
           prd.dyn{iTrialCat}.corr.goStim.nonTargetGO.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iCorr,thisINonTargetGO,etGo,CUM_PROB,alignTWindow.go);

          % Get quantile averaged dynamics of targetGO aligned on response
          [prd.dyn{iTrialCat}.corr.resp.targetGO.qX, ...
           prd.dyn{iTrialCat}.corr.resp.targetGO.qY, ...
           prd.dyn{iTrialCat}.corr.resp.targetGO.sX, ...
           prd.dyn{iTrialCat}.corr.resp.targetGO.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iCorr,iTargetGO,etResp,CUM_PROB,alignTWindow.go);

          % Get quantile averaged dynamics of nonTargetGO aligned on response
          [prd.dyn{iTrialCat}.corr.resp.nonTargetGO.qX, ...
           prd.dyn{iTrialCat}.corr.resp.nonTargetGO.qY, ...
           prd.dyn{iTrialCat}.corr.resp.nonTargetGO.sX, ...
           prd.dyn{iTrialCat}.corr.resp.nonTargetGO.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iCorr,thisINonTargetGO,etResp,CUM_PROB,alignTWindow.go);
        end

        if nError >= 1

          % Get event times of go-signals and error responses  
          etGo = repmat(stmOns(iTargetGO),nError,1);
          [etResp,iNonTargetGOError]   = min(rt(:,iError),[],1);

          % Get quantile averaged dynamics of target GO aligned on go-signal
          [prd.dyn{iTrialCat}.error.goStim.targetGO.qX, ...
           prd.dyn{iTrialCat}.error.goStim.targetGO.qY, ...
           prd.dyn{iTrialCat}.error.goStim.targetGO.sX, ...
           prd.dyn{iTrialCat}.error.goStim.targetGO.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iError,iTargetGO,etGo,CUM_PROB,alignTWindow.go);

          % Get quantile averaged dynamics of error nontarget GO aligned on go-signal
          [prd.dyn{iTrialCat}.error.goStim.nonTargetGOError.qX, ...
           prd.dyn{iTrialCat}.error.goStim.nonTargetGOError.qY, ...
           prd.dyn{iTrialCat}.error.goStim.nonTargetGOError.sX, ...
           prd.dyn{iTrialCat}.error.goStim.nonTargetGOError.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iError,iNonTargetGOError,etGo,CUM_PROB,alignTWindow.go);

          % Get quantile averaged dynamics of target GO aligned on response
          [prd.dyn{iTrialCat}.error.resp.targetGO.qX, ...
           prd.dyn{iTrialCat}.error.resp.targetGO.qY, ...
           prd.dyn{iTrialCat}.error.resp.targetGO.sX, ...
           prd.dyn{iTrialCat}.error.resp.targetGO.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iError,iTargetGO,etResp,CUM_PROB,alignTWindow.resp);

          % Get quantile averaged dynamics of error nontarget GO aligned on go-signal
          [prd.dyn{iTrialCat}.error.resp.nonTargetGOError.qX, ...
           prd.dyn{iTrialCat}.error.resp.nonTargetGOError.qY, ...
           prd.dyn{iTrialCat}.error.resp.nonTargetGOError.sX, ...
           prd.dyn{iTrialCat}.error.resp.nonTargetGOError.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iError,iNonTargetGOError,etResp,CUM_PROB,alignTWindow.resp);
        end
        
      elseif ~isempty(regexp(trialCat,'stopTrial.*', 'once'))
        
        if nError >= 1 % Stop failure trial

          % Get event times of go-signals, stop-signals and error responses  
          etGo = repmat(stmOns(iTargetGO),nError,1);
          etStop = repmat(stmOns(iTargetSTOP),nError,1);
          [etResp,iRespGO]   = min(rt(iGO,iError),[],1);
          iRespGO = iGO(iRespGO);

          % Get quantile averaged dynamics of respGO aligned on go-signal
          [prd.dyn{iTrialCat}.error.goStim.respGO.qX, ...
           prd.dyn{iTrialCat}.error.goStim.respGO.qY, ...
           prd.dyn{iTrialCat}.error.goStim.respGO.sX, ...
           prd.dyn{iTrialCat}.error.goStim.respGO.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iError,iRespGO,etGo,CUM_PROB,alignTWindow.go);

          % Get quantile averaged dynamics of targetSTOP aligned on go-signal
          [prd.dyn{iTrialCat}.error.goStim.targetSTOP.qX, ...
           prd.dyn{iTrialCat}.error.goStim.targetSTOP.qY, ...
           prd.dyn{iTrialCat}.error.goStim.targetSTOP.sX, ...
           prd.dyn{iTrialCat}.error.goStim.targetSTOP.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iError,iTargetSTOP,etGo,CUM_PROB,alignTWindow.go);

          % Get quantile averaged dynamics of respGO aligned on stop-signal
          [prd.dyn{iTrialCat}.error.stopStim.respGO.qX, ...
           prd.dyn{iTrialCat}.error.stopStim.respGO.qY, ...
           prd.dyn{iTrialCat}.error.stopStim.respGO.sX, ...
           prd.dyn{iTrialCat}.error.stopStim.respGO.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iError,iRespGO,etStop,CUM_PROB,alignTWindow.stop);

          % Get quantile averaged dynamics of targetSTOP aligned on stop-signal
          [prd.dyn{iTrialCat}.error.stopStim.targetSTOP.qX, ...
           prd.dyn{iTrialCat}.error.stopStim.targetSTOP.qY, ...
           prd.dyn{iTrialCat}.error.stopStim.targetSTOP.sX, ...
           prd.dyn{iTrialCat}.error.stopStim.targetSTOP.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iError,iTargetSTOP,etStop,CUM_PROB,alignTWindow.stop);

          % Get quantile averaged dynamics of respGO aligned on response
          [prd.dyn{iTrialCat}.error.resp.respGO.qX, ...
           prd.dyn{iTrialCat}.error.resp.respGO.qY, ...
           prd.dyn{iTrialCat}.error.resp.respGO.sX, ...
           prd.dyn{iTrialCat}.error.resp.respGO.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iError,iRespGO,etResp,CUM_PROB,alignTWindow.resp);

          % Get quantile averaged dynamics of targetSTOP aligned on response
          [prd.dyn{iTrialCat}.error.resp.targetSTOP.qX, ...
           prd.dyn{iTrialCat}.error.resp.targetSTOP.qY, ...
           prd.dyn{iTrialCat}.error.resp.targetSTOP.sX, ...
           prd.dyn{iTrialCat}.error.resp.targetSTOP.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iError,iTargetSTOP,etResp,CUM_PROB,alignTWindow.resp);

        end

        if nCorr >= 1

          % Get event times of go-signals and stop-signals
          etGo = repmat(stmOns(iTargetGO),nCorr,1);
          etStop = repmat(stmOns(iTargetSTOP),nCorr,1);

          % Get quantile averaged dynamics of target GO aligned on go-signal
          [prd.dyn{iTrialCat}.corr.goStim.targetGO.qX, ...
           prd.dyn{iTrialCat}.corr.goStim.targetGO.qY, ...
           prd.dyn{iTrialCat}.corr.goStim.targetGO.sX, ...
           prd.dyn{iTrialCat}.corr.goStim.targetGO.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iCorr,iTargetGO,etGo,CUM_PROB,alignTWindow.go);

          % Get quantile averaged dynamics of target STOP aligned on go-signal
          [prd.dyn{iTrialCat}.corr.goStim.targetSTOP.qX, ...
           prd.dyn{iTrialCat}.corr.goStim.targetSTOP.qY, ...
           prd.dyn{iTrialCat}.corr.goStim.targetSTOP.sX, ...
           prd.dyn{iTrialCat}.corr.goStim.targetSTOP.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iCorr,iTargetSTOP,etGo,CUM_PROB,alignTWindow.go);

          % Get quantile averaged dynamics of target GO aligned on stop-signal
          [prd.dyn{iTrialCat}.corr.stopStim.targetGO.qX, ...
           prd.dyn{iTrialCat}.corr.stopStim.targetGO.qY, ...
           prd.dyn{iTrialCat}.corr.stopStim.targetGO.sX, ...
           prd.dyn{iTrialCat}.corr.stopStim.targetGO.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iCorr,iTargetGO,etStop,CUM_PROB,alignTWindow.stop);

          % Get quantile averaged dynamics of target STOP aligned on stop-signal
          [prd.dyn{iTrialCat}.corr.stopStim.targetSTOP.qX, ...
           prd.dyn{iTrialCat}.corr.stopStim.targetSTOP.qY, ...
           prd.dyn{iTrialCat}.corr.stopStim.targetSTOP.sX, ...
           prd.dyn{iTrialCat}.corr.stopStim.targetSTOP.sY] ...
          = sam_get_dynamics('qav-sample100',T,z,iCorr,iTargetSTOP,etStop,CUM_PROB,alignTWindow.stop);

        end
      end
  end
end
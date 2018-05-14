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
% -------------------------------------------------------------------------

trialDur      = SAM.expt.trialDur;

nRsp          = SAM.expt.nRsp;

% 1.1.2. Model variables
% -------------------------------------------------------------------------

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

T             = (accumTWindow(1):dt:accumTWindow(2))';
p             = length(T);

% 1.1.3. Simulation variables
% -------------------------------------------------------------------------------------------------------------------------
nSim          = SAM.sim.n;
trialSimFun   = SAM.sim.fun.trial;
alignTWindow  = SAM.sim.tWindow;
rngSeedStage  = SAM.sim.rng.stage;
rngSeedId     = SAM.sim.rng.id;

% 1.1.4. Optimization variables
% -------------------------------------------------------------------------------------------------------------------------
prd           = SAM.optim.prd;

% 1.2. Specify static variables
% =========================================================================

t1            = 1;                    % First time point

switch simGoal
    case 'explore'
        % Cumulative probabilities for quantile averaging of model dynamics
        CUM_PROB = 0:0.01:1;
end

% 1.3. Specify dynamic variables
% =========================================================================

% Maximum number of stimuli and response, across conditions
maxNRsp           = max(cell2mat(nRsp(:)),[],1);

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
    
    %     if iTrialCat == nTrialCat
    %         disp('Here we are');
    %     end
    
    trialCat    = SAM.optim.obs.trialCat{iTrialCat};
    stmOns      = SAM.optim.obs.onset{iTrialCat};
    stmDur      = SAM.optim.obs.duration{iTrialCat};
    
    trueTarget  = SAM.optim.modelMat.iTarget{iTrialCat};
    trueNonTarget = SAM.optim.modelMat.iNonTarget{iTrialCat};
    
    % Identify targets inside each cell
    iTarget     = cellfun(@(in1) find(in1),trueTarget,'Uni',0);
    
    % Identify non-targets inside each cell
    iNonTarget  = cellfun(@(in1) find(in1),trueNonTarget,'Uni',0);
    
    % Number of elements to add to get to corresponding index
    nOffset     = [{0},cellfun(@(in1) numel(in1),trueTarget(1:end-1),'Uni',0)];
    
    iTargetGO       = iTarget{1} + nOffset{1};
    iNonTargetGO    = iNonTarget{1} + nOffset{1};
    iGO             = sort([iTargetGO(:);iNonTargetGO(:)]);
    
    iTargetSTOP     = iTarget{2} + nOffset{2};
    iNonTargetSTOP  = iNonTarget{2} + nOffset{2};
    iSTOP           = sort([iTargetSTOP(:);iNonTargetSTOP(:)]);
    
    % Pre-allocate response time and response arrays
    rt  = inf(sum(maxNRsp),nSim);
    
    switch simGoal
        case 'explore'
            z = nan(sum(maxNRsp),nSim,p);
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
    
    parfor iTr = 1:nSim
        %   for iTr = 1:nSim
        
        switch lower(rngSeedStage)
            case 'sam_sim_expt'
                rng(rngSeedId.Seed + iTr, 'twister');
        end
        
        % 3.2.1. Timing diagram of model inputs
        % ---------------------------------------------------------------------
        
        t0Var = randomT0Factor.*rand(1,m);
        accumOns = stmOns + T0 - T0.*t0Var;
        accumDur = stmDur + accumDurFactor.*(trialDur-stmDur-accumOns);
        
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
        % ---------------------------------------------------------------------
        
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
        
        % 4.4.1. Go trial, correct choice: only one RT, produced by target GO unit
        % -------------------------------------------------------------------
        % Criteria (all should be met):
        % - Target GO unit produced an RT
        % - Target GO unit is the only unit having produced an RT
        iGoCCorr = rt(iTargetGO,:) < Inf & sum(rt < Inf) == 1;
        
        % 4.4.2. Go trial, choice error: any RT produced by a non-target GO unit
        % -------------------------------------------------------------------
        % Criteria (all should be met):
        % - At least one non-target GO unit has produced an RT
        iGoCError = any(rt(iNonTargetGO,:) < Inf,1);
        
    elseif ~isempty(regexp(trialCat,'stopTrial.*', 'once'))
        
        if all(SAM.model.mat.endoConn.nonSelfOther == 0) % GO and STOP race
            % Stop success trial
            % -------------------------------------------------------------------
            % Criteria (*any* should be met)
            % - STOP unit produces a shorter RT than any other GO units
            % - None of the units produced an RT
            iStopICorr = rt(iTargetSTOP,:) < min(rt(iGO,:),[],1) | all(rt == Inf);
            
            % Stop failure trial, correct choice
            % -------------------------------------------------------------------
            % Criteria
            % - The fastest GO unit has a shorter RT than any STOP unit
            % - Target GO unit produced an RT
            % - Target GO unit is the only unit having produced an RT
            iStopIErrorCCorr = rt(iTargetGO,:) <= rt(iSTOP,:) & all(isinf(rt(iNonTargetGO,:)),1);
            % Stop failure trial, choice error
            % -------------------------------------------------------------------
            iStopIErrorCError = min(rt(iNonTargetGO,:),[],1) <= rt(iSTOP,:) & any(rt(iNonTargetGO,:) < Inf,1);
            
        else
            % Stop success trial
            % -------------------------------------------------------------------
            % Criteria (all should be met):
            % - None of the GO units produced an RT
            iStopICorr = all(rt(iGO,:) == Inf,1);
            
            % Note: no distinction made between trials in which the STOP unit
            % did and did not produce and RT. We may be more stringent,
            % requiring that STOP finishes, or we could look at how
            % computation of SSRT is influenced by including trials in which
            % STOP did not finish versus when these trials are not included.
            
            % Stop failure trial, correct choice
            % Criteria:
            % - Target GO unit produced an RT
            % -------------------------------------------------------------------
            iStopIErrorCCorr = rt(iTargetGO,:) < Inf & all(isinf(rt(iNonTargetGO,:)),1);
            
            % Stop failure trial, choice error
            % -------------------------------------------------------------------
            % Criteria:
            % - At least one non-target GO unit has produced an RT
            iStopIErrorCError = any(rt(iNonTargetGO,:) < Inf,1);
            
            % Note: no distinction made between trials in which the STOP unit
            % did and did not produce an RT (i.e. STOP units cannot terminate
            % the trial under 'blocked input' and 'lateral inhibition'
            % inhibition mechanisms)
            
        end
        
    end
    
    % 3.5. Log model predictions
    % =====================================================================
    
    % 3.5.1. Trial numbers
    % ---------------------------------------------------------------------
    
    if ~isempty(regexp(trialCat,'goTrial.*', 'once'))
        % Compute
        nGoCCorr                          = numel(find(iGoCCorr));
        nGoCError                         = numel(find(iGoCError));
        
        % Log
        prd.nTotal(iTrialCat)             = nGoCCorr + nGoCError;
        prd.nGoCCorr(iTrialCat)           = nGoCCorr;
        prd.nGoCError(iTrialCat)          = nGoCError;
        
    elseif ~isempty(regexp(trialCat,'stopTrial.*', 'once'))
        % Compute
        nStopICorr                        = numel(find(iStopICorr));
        nStopIErrorCCorr                  = numel(find(iStopIErrorCCorr));
        nStopIErrorCError                 = numel(find(iStopIErrorCError));
        
        % Log
        prd.nTotal(iTrialCat)             = nStopICorr + nStopIErrorCCorr + nStopIErrorCError;
        prd.nStopICorr(iTrialCat)         = nStopICorr;
        prd.nStopIErrorCCorr(iTrialCat)   = nStopIErrorCCorr;
        prd.nStopIErrorCError(iTrialCat)  = nStopIErrorCError;
    end
    
    % 3.5.2. Trial probabilities
    % ---------------------------------------------------------------------
    
    if ~isempty(regexp(trialCat,'goTrial.*', 'once'))
        % Compute
        pGoCCorr                          = nGoCCorr/nSim;
        pGoCError                         = nGoCError/nSim;
        
        % Log
        prd.pTotal(iTrialCat)             = pGoCCorr + pGoCError;
        prd.pGoCCorr(iTrialCat)           = pGoCCorr;
        prd.pGoCError(iTrialCat)          = pGoCError;
        
    elseif ~isempty(regexp(trialCat,'stopTrial.*', 'once'))
        % Compute
        pStopICorr                        = nStopICorr/nSim;
        pStopIErrorCCorr                  = nStopIErrorCCorr/nSim;
        pStopIErrorCError                 = nStopIErrorCError/nSim;
        
        % Log
        prd.pTotal(iTrialCat)             = pStopICorr + pStopIErrorCCorr + pStopIErrorCError;
        prd.pStopICorr(iTrialCat)         = pStopICorr;
        prd.pStopIErrorCCorr(iTrialCat)   = pStopIErrorCCorr;
        prd.pStopIErrorCError(iTrialCat)  = pStopIErrorCError;
        
    end
    
    % 3.5.3. Trial response times
    % ---------------------------------------------------------------------
    if ~isempty(regexp(trialCat,'goTrial.*', 'once'))
        if nGoCCorr > 0
            % Compute
            rtGoCCorr = sort(rt(iTargetGO,iGoCCorr)) - stmOns(iTargetGO);
            
            % Log
            prd.rtGoCCorr{iTrialCat} = rtGoCCorr;
        end
        
        if nGoCError > 0
            % Compute
            rtGoCError = sort(min(rt(iNonTargetGO,iGoCError),[],1)) - stmOns(iTargetGO);
            
            % Note: stimOns for iTargetGO and iNonTargetGO are the same; I use
            % iTargetGO instead of iNonTargetGO because iTargetGO is always a
            % scalar, iNonTargetGO not.
            
            % Log
            prd.rtGoCError{iTrialCat} = rtGoCError;
        end
        
    elseif ~isempty(regexp(trialCat,'stopTrial.*', 'once'))
        if nStopICorr > 0
            % Compute
            rtStopICorr = sort(rt(iSTOP,iStopICorr)) - stmOns(iSTOP);
            
            % Log
            prd.rtStopICorr{iTrialCat} = rtStopICorr;
        end
        
        if nStopIErrorCCorr > 0
            % Compute
            rtStopIErrorCCorr = sort(min(rt(iGO,iStopIErrorCCorr),[],1)) - stmOns(iTargetGO);
            
            % Log
            prd.rtStopIErrorCCorr{iTrialCat} = rtStopIErrorCCorr;
        end
        
        if nStopIErrorCError > 0
            % Compute
            rtStopIErrorCError = sort(min(rt(iGO,iStopIErrorCError),[],1)) - stmOns(iTargetGO);
            
            % Log
            prd.rtStopIErrorCError{iTrialCat} = rtStopIErrorCError;
        end
        
    end
    
    % 3.5.4. Model matrices used for simulating this trial category
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
    
    
    % 3.5.5. Model inputts
    % -----------------------------------------------------------------------------------------------------------------------
    
    % Implement this
    
    
    % 3.5.6. Model dynamics
    % -----------------------------------------------------------------------------------------------------------------------
    switch simGoal
        case 'explore'
            
            % Randomly sample a non-target GO unit
            thisINonTargetGO = randsample([iNonTargetGO(:);iNonTargetGO(:)],1);  % N.B. Repeating iNonTargetGO prevents that iNonTargetGO is a scalar, under which randsample samples from 1:iNonTargetGO
            
            if ~isempty(regexp(trialCat,'goTrial.*', 'once'))
                
                % Go correct quantile averaged dynamics
                if nGoCCorr >= 1
                    
                    % Get event times of go-signals and correct responses
                    etGo = repmat(stmOns(iTargetGO),nGoCCorr,1);
                    etResp = rt(iTargetGO,iGoCCorr);
                    
                    % Get quantile averaged dynamics of targetGO aligned on go-signal
                    [prd.dyn{iTrialCat}.goCCorr.goStim.targetGO.qX, ...
                        prd.dyn{iTrialCat}.goCCorr.goStim.targetGO.qY, ...
                        prd.dyn{iTrialCat}.goCCorr.goStim.targetGO.sX, ...
                        prd.dyn{iTrialCat}.goCCorr.goStim.targetGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iGoCCorr,iTargetGO,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of nonTargetGO aligned on go-signal
                    [prd.dyn{iTrialCat}.goCCorr.goStim.nonTargetGO.qX, ...
                        prd.dyn{iTrialCat}.goCCorr.goStim.nonTargetGO.qY, ...
                        prd.dyn{iTrialCat}.goCCorr.goStim.nonTargetGO.sX, ...
                        prd.dyn{iTrialCat}.goCCorr.goStim.nonTargetGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iGoCCorr,thisINonTargetGO,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of targetGO aligned on response
                    [prd.dyn{iTrialCat}.goCCorr.resp.targetGO.qX, ...
                        prd.dyn{iTrialCat}.goCCorr.resp.targetGO.qY, ...
                        prd.dyn{iTrialCat}.goCCorr.resp.targetGO.sX, ...
                        prd.dyn{iTrialCat}.goCCorr.resp.targetGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iGoCCorr,iTargetGO,etResp,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of nonTargetGO aligned on response
                    [prd.dyn{iTrialCat}.goCCorr.resp.nonTargetGO.qX, ...
                        prd.dyn{iTrialCat}.goCCorr.resp.nonTargetGO.qY, ...
                        prd.dyn{iTrialCat}.goCCorr.resp.nonTargetGO.sX, ...
                        prd.dyn{iTrialCat}.goCCorr.resp.nonTargetGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iGoCCorr,thisINonTargetGO,etResp,CUM_PROB,alignTWindow.go);
                end
                
                if nGoCError >= 1
                    
                    % Get event times of go-signals and error responses
                    etGo = repmat(stmOns(iTargetGO),nGoCError,1);
                    [etResp,iNonTargetGOError]   = min(rt(:,iGoCError),[],1);
                    
                    % Get quantile averaged dynamics of target GO aligned on go-signal
                    [prd.dyn{iTrialCat}.goCError.goStim.targetGO.qX, ...
                        prd.dyn{iTrialCat}.goCError.goStim.targetGO.qY, ...
                        prd.dyn{iTrialCat}.goCError.goStim.targetGO.sX, ...
                        prd.dyn{iTrialCat}.goCError.goStim.targetGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iGoCError,iTargetGO,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of error nontarget GO aligned on go-signal
                    [prd.dyn{iTrialCat}.goCError.goStim.nonTargetGOError.qX, ...
                        prd.dyn{iTrialCat}.goCError.goStim.nonTargetGOError.qY, ...
                        prd.dyn{iTrialCat}.goCError.goStim.nonTargetGOError.sX, ...
                        prd.dyn{iTrialCat}.goCError.goStim.nonTargetGOError.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iGoCError,iNonTargetGOError,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of target GO aligned on response
                    [prd.dyn{iTrialCat}.goCError.resp.targetGO.qX, ...
                        prd.dyn{iTrialCat}.goCError.resp.targetGO.qY, ...
                        prd.dyn{iTrialCat}.goCError.resp.targetGO.sX, ...
                        prd.dyn{iTrialCat}.goCError.resp.targetGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iGoCError,iTargetGO,etResp,CUM_PROB,alignTWindow.resp);
                    
                    % Get quantile averaged dynamics of error nontarget GO aligned on go-signal
                    [prd.dyn{iTrialCat}.goCError.resp.nonTargetGOError.qX, ...
                        prd.dyn{iTrialCat}.goCError.resp.nonTargetGOError.qY, ...
                        prd.dyn{iTrialCat}.goCError.resp.nonTargetGOError.sX, ...
                        prd.dyn{iTrialCat}.goCError.resp.nonTargetGOError.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iGoCError,iNonTargetGOError,etResp,CUM_PROB,alignTWindow.resp);
                end
                
            elseif ~isempty(regexp(trialCat,'stopTrial.*', 'once'))
                
                if nStopIErrorCCorr >= 1 % Stop failure trial, correct choice
                    
                    % Get event times of go-signals, stop-signals and error responses
                    etGo = repmat(stmOns(iTargetGO),nStopIErrorCCorr,1);
                    etStop = repmat(stmOns(iTargetSTOP),nStopIErrorCCorr,1);
                    [etResp,iRespGO]   = min(rt(iGO,iStopIErrorCCorr),[],1);
                    iRespGO = iGO(iRespGO);
                    
                    % Get quantile averaged dynamics of respGO aligned on go-signal
                    [prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.respGO.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.respGO.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.respGO.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.respGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCCorr,iRespGO,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of nonTargetGO aligned on go-signal
                    [prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.nonTargetGO.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.nonTargetGO.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.nonTargetGO.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.nonTargetGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCCorr,thisINonTargetGO,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of targetSTOP aligned on go-signal
                    [prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.targetSTOP.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.targetSTOP.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.targetSTOP.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.goStim.targetSTOP.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCCorr,iTargetSTOP,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of respGO aligned on stop-signal
                    [prd.dyn{iTrialCat}.stopIErrorCCorr.stopStim.respGO.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.stopStim.respGO.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.stopStim.respGO.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.stopStim.respGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCCorr,iRespGO,etStop,CUM_PROB,alignTWindow.stop);
                    
                    % Get quantile averaged dynamics of targetSTOP aligned on stop-signal
                    [prd.dyn{iTrialCat}.stopIErrorCCorr.stopStim.targetSTOP.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.stopStim.targetSTOP.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.stopStim.targetSTOP.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.stopStim.targetSTOP.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCCorr,iTargetSTOP,etStop,CUM_PROB,alignTWindow.stop);
                    
                    % Get quantile averaged dynamics of respGO aligned on response
                    [prd.dyn{iTrialCat}.stopIErrorCCorr.resp.respGO.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.resp.respGO.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.resp.respGO.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.resp.respGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCCorr,iRespGO,etResp,CUM_PROB,alignTWindow.resp);
                    
                    % Get quantile averaged dynamics of targetSTOP aligned on response
                    [prd.dyn{iTrialCat}.stopIErrorCCorr.resp.targetSTOP.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.resp.targetSTOP.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.resp.targetSTOP.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCCorr.resp.targetSTOP.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCCorr,iTargetSTOP,etResp,CUM_PROB,alignTWindow.resp);
                    
                end
                
                if nStopIErrorCError >= 1 % Stop failure trial, correct choice
                    
                    % Get event times of go-signals, stop-signals and error responses
                    etGo = repmat(stmOns(iTargetGO),nStopIErrorCError,1);
                    etStop = repmat(stmOns(iTargetSTOP),nStopIErrorCError,1);
                    [etResp,iRespGO]   = min(rt(iGO,iStopIErrorCError),[],1);
                    iRespGO = iGO(iRespGO);
                    
                    % Get quantile averaged dynamics of respGO aligned on go-signal
                    [prd.dyn{iTrialCat}.stopIErrorCError.goStim.respGO.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.goStim.respGO.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.goStim.respGO.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.goStim.respGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCError,iRespGO,etGo,CUM_PROB,alignTWindow.go);
                    
                    %           % Get quantile averaged dynamics of nonTargetGO aligned on go-signal
                    %           [prd.dyn{iTrialCat}.stopIErrorCError.goStim.nonTargetGO.qX, ...
                    %            prd.dyn{iTrialCat}.stopIErrorCError.goStim.nonTargetGO.qY, ...
                    %            prd.dyn{iTrialCat}.stopIErrorCError.goStim.nonTargetGO.sX, ...
                    %            prd.dyn{iTrialCat}.stopIErrorCError.goStim.nonTargetGO.sY] ...
                    %           = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCError,thisINonTargetGO,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of targetSTOP aligned on go-signal
                    [prd.dyn{iTrialCat}.stopIErrorCError.goStim.targetSTOP.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.goStim.targetSTOP.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.goStim.targetSTOP.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.goStim.targetSTOP.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCError,iTargetSTOP,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of respGO aligned on stop-signal
                    [prd.dyn{iTrialCat}.stopIErrorCError.stopStim.respGO.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.stopStim.respGO.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.stopStim.respGO.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.stopStim.respGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCError,iRespGO,etStop,CUM_PROB,alignTWindow.stop);
                    
                    % Get quantile averaged dynamics of targetSTOP aligned on stop-signal
                    [prd.dyn{iTrialCat}.stopIErrorCError.stopStim.targetSTOP.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.stopStim.targetSTOP.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.stopStim.targetSTOP.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.stopStim.targetSTOP.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCError,iTargetSTOP,etStop,CUM_PROB,alignTWindow.stop);
                    
                    % Get quantile averaged dynamics of respGO aligned on response
                    [prd.dyn{iTrialCat}.stopIErrorCError.resp.respGO.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.resp.respGO.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.resp.respGO.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.resp.respGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCError,iRespGO,etResp,CUM_PROB,alignTWindow.resp);
                    
                    % Get quantile averaged dynamics of targetSTOP aligned on response
                    [prd.dyn{iTrialCat}.stopIErrorCError.resp.targetSTOP.qX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.resp.targetSTOP.qY, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.resp.targetSTOP.sX, ...
                        prd.dyn{iTrialCat}.stopIErrorCError.resp.targetSTOP.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopIErrorCError,iTargetSTOP,etResp,CUM_PROB,alignTWindow.resp);
                    
                end
                
                if nStopICorr >= 1
                    
                    % Get event times of go-signals and stop-signals
                    etGo = repmat(stmOns(iTargetGO),nStopICorr,1);
                    etStop = repmat(stmOns(iTargetSTOP),nStopICorr,1);
                    
                    % Get quantile averaged dynamics of target GO aligned on go-signal
                    [prd.dyn{iTrialCat}.stopICorr.goStim.targetGO.qX, ...
                        prd.dyn{iTrialCat}.stopICorr.goStim.targetGO.qY, ...
                        prd.dyn{iTrialCat}.stopICorr.goStim.targetGO.sX, ...
                        prd.dyn{iTrialCat}.stopICorr.goStim.targetGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopICorr,iTargetGO,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of nontarget GO aligned on go-signal
                    [prd.dyn{iTrialCat}.stopICorr.goStim.nonTargetGO.qX, ...
                        prd.dyn{iTrialCat}.stopICorr.goStim.nonTargetGO.qY, ...
                        prd.dyn{iTrialCat}.stopICorr.goStim.nonTargetGO.sX, ...
                        prd.dyn{iTrialCat}.stopICorr.goStim.nonTargetGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopICorr,thisINonTargetGO,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of target STOP aligned on go-signal
                    [prd.dyn{iTrialCat}.stopICorr.goStim.targetSTOP.qX, ...
                        prd.dyn{iTrialCat}.stopICorr.goStim.targetSTOP.qY, ...
                        prd.dyn{iTrialCat}.stopICorr.goStim.targetSTOP.sX, ...
                        prd.dyn{iTrialCat}.stopICorr.goStim.targetSTOP.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopICorr,iTargetSTOP,etGo,CUM_PROB,alignTWindow.go);
                    
                    % Get quantile averaged dynamics of target GO aligned on stop-signal
                    [prd.dyn{iTrialCat}.stopICorr.stopStim.targetGO.qX, ...
                        prd.dyn{iTrialCat}.stopICorr.stopStim.targetGO.qY, ...
                        prd.dyn{iTrialCat}.stopICorr.stopStim.targetGO.sX, ...
                        prd.dyn{iTrialCat}.stopICorr.stopStim.targetGO.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopICorr,iTargetGO,etStop,CUM_PROB,alignTWindow.stop);
                    
                    % Get quantile averaged dynamics of target STOP aligned on stop-signal
                    [prd.dyn{iTrialCat}.stopICorr.stopStim.targetSTOP.qX, ...
                        prd.dyn{iTrialCat}.stopICorr.stopStim.targetSTOP.qY, ...
                        prd.dyn{iTrialCat}.stopICorr.stopStim.targetSTOP.sX, ...
                        prd.dyn{iTrialCat}.stopICorr.stopStim.targetSTOP.sY] ...
                        = sam_get_dynamics('qav-sample100',T,z,iStopICorr,iTargetSTOP,etStop,CUM_PROB,alignTWindow.stop);
                    
                end
            end
    end
end
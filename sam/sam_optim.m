function SAM = sam_optim(SAM,iStartVal)
% Optimizes stochastic accumulator model to account for observations
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% SAM_OPTIM; 
%  
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Sat 21 Sep 2013 12:53:48 CDT by bram 
% $Modified: Sat 21 Sep 2013 19:40:18 CDT by bram

 
% CONTENTS 
% 1.PROCESS INPUTS AND SPECIFY VARIABLES 
%   1.1.Process inputs
%   1.2. Pre-allocate empty arrays
% 2.SPECIFY PRECURSOR AND PARAMETER-INDEPENDENT MODEL MATRICES
% 3.CHARACTERIZE OBSERVED DATA
%   3.1. Organize observations
%   3.2. Compute response time bin statistics
% 4.OPTIMIZE MODEL
%   4.1.Seed the random number generator
%   4.2.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS AND SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% 1.1. Process inputs
% ========================================================================= 

% Scope of the simulation
optimScope      = SAM.sim.scope;

% Stage at which the random number generator is seeded
rngSeedStage  = SAM.sim.rng.stage;

% Random number generator seed identifier
rngSeedId     = SAM.sim.rng.id;

% Solver type
solverType    = SAM.optim.solver.type;

% Solver options
solverOpts    = SAM.optim.solver.opts;

% Job ID
jobIDDigits   = SAM.compCluster.jobID;

% Starting values
% ---------------------------------------------------------------------
X0            = SAM.optim.x0(iStartVal,:);

% Constraints
% ---------------------------------------------------------------------
LB            = SAM.optim.constraint.bound.LB;
UB            = SAM.optim.constraint.bound.UB;
A             = SAM.optim.constraint.linear.A;
b             = SAM.optim.constraint.linear.b;
Aeq           = SAM.optim.constraint.linear.Aeq;
beq           = SAM.optim.constraint.linear.beq
nonlcon       = SAM.optim.constraint.nonlinear.nonLinCon;

% Cost function 
% ---------------------------------------------------------------------
costFun       = SAM.optim.cost.fun;

% Logging
% ---------------------------------------------------------------------

% Iteration log file
iterLogFile   = SAM.optim.log.iterLogFile;

% Iteration lof frequency
iterLogFreq   = SAM.optim.log.iterLogFreq;

% Final log file
finalLogFile  = SAM.optim.log.finalLogFile;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. SEED THE RANDOM NUMBER GENERATOR (OPTIONAL)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

switch lower(rngSeedStage)
  case 'sam_optim'

    % Note: MEX functions stay in memory until they are cleared.
    % Seeding of the random number generator should be accompanied by 
    % clearing MEX functions.

    clear(char(SAM.sim.fun.trial));
    rng(rngSeedId);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% #. START THE PARALLEL POOL
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% In unique directory to prevent collision of parallel jobs
% e.g. see: http://www.mathworks.com/matlabcentral/answers/97141-why-am-i-unable-to-start-a-local-matlabpool-from-multiple-matlab-sessions-that-use-a-shared-preferen
c = parcluster();
if isfield(SAM,'compCluster')
	c.NumWorkers = SAM.compCluster.nProcessors;
else
    c.NumWorkers = 1;
end
[~,homeDir] = system('echo $HOME');
homeDir = strtrim(homeDir);
release = version('-release')
tempDir = fullfile(homeDir,'.matlab','local_cluster_jobs',release);
if exist(tempDir) ~= 7
    mkdir(tempDir)
end
t = fullfile([tempname(tempDir),'_JobID_',jobIDDigits]);
mkdir(t);
c.JobStorageLocation=t;
tWait = 1+60*rand();
pause(tWait);
myPool = parpool(c);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. OPTIMIZE MODEL
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 2.1. Switch between solvers
% =========================================================================

switch lower(solverType)
  % 2.1.1. Simplex
  % -----------------------------------------------------------------------
  case 'fminsearchbnd'
    
    history = nan(solverOpts.MaxIter + 1,numel(X0) + 1);
    
    solverOpts.OutputFcn        = @myoutput;
    
    [X, ...
     fVal, ...
     exitFlag, ...
     solverOutput] ...
     ...
     = fminsearchbnd ...
     (@(X)  ...
     costFun ...
     (X, ...
      SAM), ...
      ...
      X0, ...
      LB, ...
      UB, ...
      solverOpts);
    
      % Also get alternative cost value
      [~,altCost] = sam_cost(X,SAM);
  
      % Save the final log file
      save(finalLogFile,'X','fVal','altCost','exitFlag','solverOutput','history');
      
      % Remove iteration log file (history is also saved in final log file)
      delete(iterLogFile);
%       
%       varargout{1} = X;
%       varargout{2} = fVal;
%       varargout{3} = exitFlag;
%       varargout{4} = solverOutput;
%       varargout{5} = history;
%   
  case 'fminsearchcon'
    
    history = nan(solverOpts.MaxIter + 1,numel(X0) + 1);
    
    solverOpts.OutputFcn        = @myoutput;
    
    tS = tic;
    
    [X, ...
     fVal, ...
     exitFlag, ...
     solverOutput] ...
     ...
     = fminsearchcon ...
     (@(X)  ...
     costFun ...
     (X, ...
      SAM), ...
      ...
      X0, ...
      LB, ...
      UB, ...
      A, ...
      b, ...
      nonlcon, ...
      solverOpts);
    
    tElapse = toc(tS);
    
    % Also get alternative cost value
    [~,altCost] = sam_cost(X,SAM);
  
    % Save the final log file
    save(finalLogFile,'X','fVal','altCost','exitFlag','solverOutput','history');
    
    % Remove unused lines (containing only NaNs)
    iFirstNanLine = find(all(isnan(history),2),1,'first');
    history = history(1:iFirstNanLine-1,:);

    % Save the final log file
    save(finalLogFile,'X','fVal','exitFlag','solverOutput','history','tElapse');

    % Remove iteration log file (history is also saved in final log file)
    delete(iterLogFile);
      
%       varargout{1} = X;
%       varargout{2} = fVal;
%       varargout{3} = exitFlag;
%       varargout{4} = solverOutput;
%       varargout{5} = history;
 
  % 2.1.2. Differential evolution
  % -----------------------------------------------------------------------
  case 'de'
  
  % 2.1.3. Genetic algorithm
  % -----------------------------------------------------------------------
  case 'ga'

%     nX = numel(LB);
%     
%     [X, ...
%      fVal, ...
%      exitFlag, ...
%      solverOutput, ...
%      pop, ...
%      cost] = ga(@(X)  costFun ...
%                      (X, ...
%                       SAM, ...
%                       ), ...
%                       ...
%                       nX, ...
%                       linConA, ...
%                       linConB, ...
%                       [], ...
%                       [], ...
%                       LB, ...
%                       UB, ...
%                       nonLinCon, ...
%                       solverOpts);

  % 2.1.4. Simulated annealing
  % -----------------------------------------------------------------------
  case 'sa'
% 
%     [X, ...
%      fVal, ...
%      exitFlag, ...
%      solverOutput] = simulannealbnd(@(X)  costFun(simGoal, ...
%                                                   X0, ...
%                                                   STATE, ...
%                                                   obs, ...
%                                                   prd, ...
%                                                   SAM, ...
%                                                   nSim, ...
%                                                   simFun, ...
%                                                   optimScope, ...
%                                                   VCor, ...
%                                                   VIncor, ...
%                                                   S, ...
%                                                   terminate, ...
%                                                   blockInput, ...
%                                                   latInhib, ...
%                                                   doPlot), ...
%                                                   ...
%                                                   X0, ...
%                                                   LB, ...
%                                                   UB, ...
%                                                   solver.options);
end

% Shut down the parallel pool
% =========================================================================
delete(myPool);

% Remove temporary directory
system(['rm -r ',t]);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3. OUTPUT
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

SAM.estim.X             = X;
SAM.estim.fVal          = fVal;
SAM.estim.fValAlt       = altCost;
SAM.estim.exitFlag      = exitFlag;
SAM.estim.solverOutput  = solverOutput;
SAM.estim.tElapse       = tElapse;

function stop = myoutput(x,optimvalues,state)
  stop = false;
  if strncmpi(state,'iter',4)
    history(optimvalues.iteration + 1,:) = [optimvalues.fval,x];
    
    % Save the iteration log file if
    if ismultiple(optimvalues.iteration,iterLogFreq)
      save(iterLogFile,'history');
    end
  end
end

  function out = ismultiple(iter,freq)
    out = freq*round(double(iter)/freq) == iter;
  end
end

function solverOpts = sam_get_solver_opts(solverType)
% SAM_GET_SOLVER_OPTS <Synopsis of what this function does> 
%  
% DESCRIPTION 
% Specifies the optimization structure |solverOpts| for each of the
% following solver types:
% |'de'|              - Differential evolution
% |'fminsearchbnd'|   - Simplex w/ bounds
% |'fminsearchcon'|   - Simplex w/ bounds and (non)linear constraints
% |'ga'|              - Genetic algorithm
% |'sa'|              - Simulated annealing
%  
% SYNTAX 
% solverOpts = SAM_GET_SOLVER_OPTS(solverType);
%  
% EXAMPLES 
% solverOpts = SAM_GET_SOLVER_OPTS('fminsearchcon');
% 
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Wed 12 Feb 2014 12:42:28 CST by bram 
% $Modified: Wed 12 Feb 2014 12:42:28 CST by bram 

% 1.1. Get default options
% ========================================================================= 
switch lower(solverType)
    
    % Differential evolution
    case 'de'
      error('Don''t know which option structure to use. Implement this.');
    
    % Simplex with bounds  
    case 'fminsearchbnd'
      solverOpts            = optimset(@fminsearch);
    
    % Simplex with bounds and linear and nonlinear constraints
    case 'fminsearchcon'
      solverOpts            = optimset(@fminsearch);
  
    % Genetic algorithm
    case 'ga'
      solverOpts            = gaoptimset(@ga);
    
    % Simulated annealing
    case 'sa'
      solverOpts            = saoptimset(@simulannealbnd);
end

% 1.2. Modify solver-generic options
% ========================================================================= 

% Set display to iteration and plot functions off
solverOpts.Display        = 'iter';
solverOpts.PlotFcns       = {[]};

% 1.3. Modify solver-generic options
% ========================================================================= 

switch lower(solverType)
      
  % Differential evolution
  case 'de'
    %     OPTIONS:
    %  These options may be specified using parameter, value pairs or by
    %  passing a structure. Defaults are shown in parentheses.
    %   popsize        - total number of individuals. (100)
    %   generations    - (10)
    %   strategy       - mutation strategy (see MUTATE for options)
    %   step_weight    - stepsize weight (between 0 and 2) to apply to
    %                    differentials when mutating parameters (0.85)
    %   crossover      - crossover probability constant (between 0 and
    %                    1).  Percentage of random new mutated
    %                    parameters to use in the new population (1)
    %   range_bound    - boolean indicating whether parameters are
    %                    strictly bound by the values in ranges (true)
    %   start_file     - path to a MAT-file containing two variables:
    %                     fitness
    %                     parameters
    %   collect_erfvec - if true, erf values will be saved. (false)

  % Simplex with bounds  
  case 'fminsearchbnd'

      solverOpts.MaxFunEvals    = 150000;
      solverOpts.MaxIter        = 50;
      solverOpts.TolFun         = 1e-5;
      solverOpts.TolX           = 1e-5;

  % Simplex with bounds and linear and nonlinear constraints    
  case 'fminsearchcon'

    solverOpts.MaxFunEvals      = '500*numberofvariables';
    solverOpts.MaxIter          = '500*numberofvariables';
    solverOpts.TolFun           = 1e-6;
    solverOpts.TolX             = 1e-6;

  % Genetic algorithm  
  case 'ga'

    popSize=30;
    nOfColony=1;
    popVec=ones(1,nOfColony)*popSize;

    solverOpts.PopInitRange     = [LB;UB];
    solverOpts.PopulationSize   = popVec;
    solverOpts.EliteCount       = floor(popSize*.2);
    solverOpts.Generations      = 2;
    solverOpts.CrossoverFcn     = {@crossoverscattered};
    solverOpts.MutationFcn      = {@mutationadaptfeasible};
    solverOpts.SelectionFcn     = {@selectionroulette};
    solverOpts.Vectorized       = 'off';
    solverOpts.PlotFcns         = {@gaplotbestf,@gaplotbestindiv};

  % Simulated annealing
  case 'sa'

    solverOpts.PlotFcns         = {@optimplotx, ...
                                   @optimplotfval, ...
                                   @optimplotfunccount};
    
end
function [VCor,VIncor,S,terminate,blockInput,latInhib] = sam_spec_general_mat(SAM)
% Specifies precursor and parameter-independent model matrices
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% [VCor,VIncor,S,terminate,blockInput,latInhib] = SAM_SPEC_GENERAL_MAT(SAM); 
%
% VCor            
% VIncor
% S
% terminate
% blockInput
% latInhib
%
%
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Sat 21 Sep 2013 12:09:45 CDT by bram 
% $Modified: Sat 21 Sep 2013 12:09:45 CDT by bram 

% CONTENTS 
% 1.PROCESS INPUTS AND SPECIFY VARIABLES
%   1.1.Process inputs
%   1.2.Specify dynamic variables
% 2.SPECIFY PRECURSOR MODEL MATRICES
%   2.1.Precursor matrix for accumulation rates to target units
%   2.2.Precursor matrix for accumulation rates to nontarget units
%   2.3.Precursor matrix for extrinsic and intrinsic noise levels
% 3.SPECIFY PARAMETER-INDEPENDENT MODEL MATRICES
%   3.1.Termination matrix
%   3.2.Blocked-input matrix
%   3.3.Lateral inhibition matrix

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS AND SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1.1. Process inputs
% =========================================================================
          
% Number of conditions
nCnd      = SAM.des.expt.nCnd;

% Number of stop-signal delays
nSsd      = SAM.des.expt.nSsd;

% Scope of simulations
optimScope  = SAM.sim.scope;

% Type of inhibition mechanism
inhibMechType = SAM.des.inhibMech.type;

% Number of GO and STOP units
nGo       = SAM.des.nGO;
nStop     = SAM.des.nSTOP;

% Indices of GO, GO target, and STOP units
iGO       = SAM.des.iGO;
iGOT      = SAM.des.iGOT;
iGONT     = SAM.des.iGONT;

% 1.2. Specify dynamic variables
% =========================================================================

switch lower(optimScope)
  case 'go'
    % Number of units
    N = nGo;
    
    % Number of model inputs (go and stop stimuli)
    M = nGo;
    
  case {'stop','all'}
    
    % Number of units
    N = [nGo nStop];
    
    % Number of model inputs (go and stop stimuli)
    M = [nGo nStop];
    
end

% Cell array of true arrays for GO and STOP units
trueN     = arrayfun(@(x) true(x,1),N,'Uni',0);

% Cell array of true arrays for Go and Stop stimuli
trueM     = arrayfun(@(x) true(x,1),M,'Uni',0);

% Cell array of false arrays for Go and Stop stimuli
falseM    = arrayfun(@(x) false(x,1),M,'Uni',0);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. SPECIFY PRECURSOR MODEL MATRICES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 2.1. Precursor matrix for accumulation rates to target units
% =========================================================================

switch lower(optimScope)
  case 'go'
    VCor    = cell(nCnd,1);
  case {'stop','all'}
    VCor    = cell(nCnd,1 + nSsd);
    
    for iCnd = 1:nCnd
      stopTrMotif = falseM;
      
      % GO target unit index
      stopTrMotif{1}(iGOT{iCnd}) = true;
      
      % STOP target unit index
      stopTrMotif{2}(1) = true;
      
      % Fill in precursor matrix for accumulation rates to target unit(s) on
      % Stop trials
      VCor(iCnd,2:end) = cellfun(@(a) blkdiag(stopTrMotif{:}),VCor(iCnd,2:end),'Uni',0);
    
    end
end

for iCnd = 1:nCnd
  goTrMotif = falseM;
  
  % GO target unit index
  goTrMotif{1}(iGOT{iCnd}) = true;
 
  % Fill in precursor matrix for accumulation rates to target unit(s) on
  % Go trials
  VCor(iCnd,1)     = cellfun(@(a) blkdiag(goTrMotif{:}),VCor(iCnd,1),'Uni',0);
end

% 2.2. Precursor matrix for accumulation rates to nontarget units
% =========================================================================

switch lower(optimScope)
  case 'go'
    VIncor  = cell(nCnd,1);
  case {'stop','all'}
    VIncor  = cell(nCnd,1 + nSsd);
end

for iCnd = 1:nCnd
  goTrMotif = falseM;
  
  % GO non-target unit indices
  goTrMotif{1}(iGONT{iCnd}) = true;
  
  % Fill in precursor matrix for accumulation rates to nontarget unit(s) on
  % Go trials
  VIncor(iCnd,:) = cellfun(@(a) blkdiag(goTrMotif{:}),VIncor(iCnd,1),'Uni',0);
end

% 2.3. Precursor matrix for extrinsic and intrinsic noise levels
% =========================================================================

switch lower(optimScope)
  case 'go'
    S = cell(nCnd,1);

  case {'stop','all'}
    S = cell(nCnd,1 + nSsd);

    for iCnd = 1:nCnd
      stopTrMotif = falseM;
      
      % GO unit indices
      stopTrMotif{1}(iGO{iCnd}) = true;
      
      % STOP unit index
      stopTrMotif{2}(1) = true;

      % Fill in precursor matrix for noise on Stop trials
      S(iCnd,2:end) = cellfun(@(a) blkdiag(stopTrMotif{:}),S(iCnd,2:end),'Uni',0);
      
    end
    
end

for iCnd = 1:nCnd
  goTrMotif = falseM;
  
  % GO unit indices
  goTrMotif{1}(iGO{iCnd}) = true;
  
  % Fill in precursor matrix for noise on Go trials
  S(iCnd,1)       = cellfun(@(a) blkdiag(goTrMotif{:}),S(iCnd,1),'Uni',0);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3. SPECIFY PARAMETER-INDEPENDENT MODEL MATRICES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 3.1. Termination matrix
% =========================================================================

switch lower(optimScope)
  case 'go'
    % GO units can terminate the accumulation process
    terminate = blkdiag(trueN{:})*true';
  case {'stop','all'}
    switch inhibMechType
      case 'race'
        % GO and STOP units can terminate the accumulation process
        terminate = blkdiag(trueN{:})*[true true]';
      case 'bi'
        % Only GO units can terminate the accumulation process
        terminate = blkdiag(trueN{:})*[true false]';
      case 'li'
        % Only GO units can terminate the accumulation process
        terminate = blkdiag(trueN{:})*[true false]';
    end
end

% Make sure that terminateis a logical array
terminate     = logical(terminate);

% 3.2. Blocked-input matrix
% =========================================================================  

switch lower(optimScope)
  case 'go'
    % GO units can terminate the accumulation process
     blockInput = (blkdiag(trueM{:})*false') * ...
                  (blkdiag(trueN{:})*false')';
  case {'stop','all'}
    switch inhibMechType
      case 'race'
         % Input is never blocked
         blockInput = (blkdiag(trueM{:})*[false false]') * ...
                      (blkdiag(trueN{:})*[false false]')';
      case 'bi'
        % Stop units block input to Go units when they reach threshold
        blockInput = (blkdiag(trueM{:})*[true false]') * ...
                     (blkdiag(trueN{:})*[false true]')';
      case 'li'
        % Input is never blocked
        blockInput = (blkdiag(trueM{:})*[false false]') * ...
                     (blkdiag(trueN{:})*[false false]')';
    end
end

% Make sure that blockInput is a logical array
blockInput    = logical(blockInput);

% 3.3. Lateral inhibition matrix
% =========================================================================    

switch lower(optimScope)
  case 'go'
    % No STOP units, so there are no connections that changes when eithe accumulator reaches threshold
    latInhib = (blkdiag(trueN{:})*false') *  ...
               (blkdiag(trueN{:})*false')'; ...
  case {'stop','all'}
    switch inhibMechType
      case 'race'
        % STOP=>GO and GO=>STOP connections do not change once GO and STOP
        % units reach threshold, respecitively
        latInhib = (blkdiag(trueN{:})*[false false]') *  ...   % To Go units
                   (blkdiag(trueN{:})*[false false]')' ...     % From Stop units
                   + ...                                      % AND
                   (blkdiag(trueN{:})*[false false]') *  ...   % From Go units
                   (blkdiag(trueN{:})*[false false]')';     % To Stop units
      case 'bi'
        % STOP=>GO and GO=>STOP connections do not change once GO and STOP
        % units reach threshold, respecitively
        latInhib = (blkdiag(trueN{:})*[false false]') *  ...   % To Go units
                   (blkdiag(trueN{:})*[false false]')' ...     % From Stop units
                   + ...                                      % AND
                   (blkdiag(trueN{:})*[false false]') *  ...   % From Go units
                   (blkdiag(trueN{:})*[false false]')';     % To Stop units
      case 'li'
        % Lateral inhibition from STOP=>GO is contingent on STOP reaching 
        % threshold
        % Lateral inhibition from GO=>STOP is present at all times
        
        latInhib = (blkdiag(trueN{:})*[true false]') *  ...   % To Go units
                   (blkdiag(trueN{:})*[false true]')' ...     % From Stop units
                   + ...                                      % AND
                   (blkdiag(trueN{:})*[false false]') *  ...   % From Go units
                   (blkdiag(trueN{:})*[false false]')';        % To Stop units
    end
end

% Make sure that latInhib is a logical array
latInhib      = logical(latInhib);
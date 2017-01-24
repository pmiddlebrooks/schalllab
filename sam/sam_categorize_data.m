function obs = sam_categorize_data(SAM)
% SAM_CATEGORIZE_DATA Categorizes behavioral data into different trial types
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% obs = SAM_CATEGORIZE_DATA(file); 
% 
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Wed 12 Mar 2014 14:16:31 CDT by bram 
% $Modified: Wed 12 Mar 2014 14:16:31 CDT by bram 

% CONTENTS
% 1. PROCESS INPUTS AND SPECIFY VARIABLES
% 2. LOAD DATA
% 3. DEFINE TRIAL CATEGORIES AND TAGS
% 4. CATEGORIZE DATA
% 5. COMPUTE DESCRIPTIVES


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESSING INPUTS AND SPECIFYING VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1.1. Process inputs
% =========================================================================

file              = SAM.io.behavFile;

nStm              = SAM.expt.nStm;                                            
nRsp              = SAM.expt.nRsp;
nCnd              = SAM.expt.nCnd;                                            
nSsd              = SAM.expt.nSsd;                                            
stmOns            = SAM.expt.stmOns;                                        
stmDur            = SAM.expt.stmDur;                                       

dt                = SAM.model.accum.dt;
modelToFit        = SAM.model.variants.toFit;

optimScope        = SAM.sim.scope;

cumProb           = SAM.optim.cost.stat.cumProb;
minBinSize        = SAM.optim.cost.stat.minBinSize;

% 1.2. Dynamic variables
% =========================================================================
    
% Miscellaneous
% -------------------------------------------------------------------------

% Maximum number of stimuli and response, across conditions
maxNStm           = max(cell2mat(nStm(:)),[],1);
maxNRsp           = max(cell2mat(nRsp(:)),[],1);

taskFactors       = [maxNStm;maxNRsp;nCnd,nCnd];

trueNStm          = arrayfun(@(x) true(x,1),maxNStm,'Uni',0);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. LOAD BEHAVIORAL DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the variable data from the behavior file
load(file,'data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. DEFINE TRIAL CATEGORIES AND TAGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Categories and tags depend on the scope of the optimization
% 'go'      - optimizing GO parameters involves no-signal trials only 
% 'stop'    - optimizing STOP parameters involves stop-signal trials only 
% 'all'     - optimizing GO and STOP parameters involves no-signal and stop-signal trials

% 3.1. Identify the signatures and the number of factors varying
% =========================================================================

% 3.1.1. Signatures
% -------------------------------------------------------------------------
signatureGO         = any(modelToFit.features(:,:,1),2);
signatureGO(2) = 1;  % constrain the model to fit both responses (so a version that varies drift rate only across conditions will also fit both response sides- otherwise it won't)
switch lower(optimScope)
  case {'all','stop'}
    signatureSTOP   = any(modelToFit.features(:,:,2),2);
end

% 3.1.2. Number of factors varying
% -------------------------------------------------------------------------
nFactGO             = numel(taskFactors(signatureGO,1));
nFactGOStmRsp       = numel(taskFactors(signatureGO(1:2),1));

switch lower(optimScope)
  case {'stop','all'}
    nFactSTOP       = numel(taskFactors(signatureSTOP,1));
    nFactSTOPStmRsp = numel(taskFactors(signatureSTOP(1:2),1));
end

% 3.2. Specify design matrix containing factor settings
% =========================================================================
% Go trials will at least be broken down by condition. If GO parameters 
% vary accross stimuli or responses, then trials are further broken down
% according to these factors. 
%
% Stop trials will will at least be broken down by stop-signal delay and
% condition. If GO or STOP parameters vary accross stimuli or responses,
% then trials are further broken down according to these factors.

switch lower(optimScope)
  case {'go','all'}
    goDesMat        = fullfact([taskFactors(signatureGO(1:2),1);nCnd])';
end

switch lower(optimScope)
  case {'stop','all'}
    if ~isequal(signatureGO(1:2),[0 0]') && ...
       ~isequal(signatureSTOP(1:2),[0 0]')
      % GO parameters vary across stimuli and/or responses
      % STOP parameters vary across stimuli and/or responses
      stopDesMat    = fullfact([nSsd; ...
                                taskFactors(signatureGO(1:2),1); ...
                                taskFactors(signatureSTOP(1:2),2); ...
                                nCnd])';
      
    elseif ~isequal(signatureGO(1:2),[0 0]') && ...
            isequal(signatureSTOP(1:2),[0 0]')  
      % GO parameters vary across stimuli and/or responses
      % STOP parameters do not vary across stimuli and/or responses
      stopDesMat    = fullfact([nSsd; ...
                                taskFactors(signatureGO(1:2),1); ...
                                nCnd])';
      
    elseif isequal(signatureGO(1:2),[0 0]') && ...
          ~isequal(signatureSTOP(1:2),[0 0]')  
      % GO parameters do not vary across stimuli and/or responses
      % STOP parameters vary across stimuli and/or responses
      stopDesMat    = fullfact([nSsd; ...
                                taskFactors(signatureSTOP(1:2),1); ...
                                nCnd])';
      
    elseif isequal(signatureGO(1:2),[0 0]') && ...
           isequal(signatureSTOP(1:2),[0 0]')  
      % GO parameters do not vary across stimuli and/or responses
      % STOP parameters do not vary across stimuli and/or responses
      stopDesMat    = fullfact([nSsd;nCnd])';
      
    end
end

% 3.3. Specify inputs for trial tag functions
% =========================================================================
% For the trial tags, the variation across stimuli and responses, if any,
% is omitted. Also, the design matrices are converted to cell arrays.

% Remove any variation across stimuli or responses from trial tag data

switch lower(optimScope)
  case {'go','all'}
    goTrialCat                      = goDesMat;
    if nFactGOStmRsp > 0
      goTrialCat(1:nFactGOStmRsp,:) = [];
    end
end

switch lower(optimScope)
  case {'stop','all'}
    stopTrialCat                    = stopDesMat;
    if nFactGOStmRsp > 0 || nFactSTOPStmRsp > 0
      iRemove                       = setdiff(1:size(stopTrialCat), ...
                                             [1,size(stopTrialCat)]);
      stopTrialCat(iRemove,:)       = [];
    end
end

% 3.3.2. Convert to cell arrays
% -------------------------------------------------------------------------

switch lower(optimScope)
  case {'go','all'}
    if numel(goTrialCat) == 1
      % If there is just one trial category (1 condition)
      goTrialCatCell    = {1};
    else
      goTrialCatCell    = mat2cell(goTrialCat, ...
                                   size(goTrialCat,1), ...
                                   ones(size(goTrialCat,2),1));
    end
end

switch lower(optimScope)
  case {'stop','all'}
    if numel(stopTrialCat) == 2
      % If there is just one trial category (1 ssd, 1 condition)
      stopTrialCatCell  = {1};
    else
      stopTrialCatCell  = mat2cell(stopTrialCat, ...
                                   size(stopTrialCat,1), ...
                                   ones(size(stopTrialCat,2),1));
    end
end

% 3.4. Specify inputs for parameter tag functions
% =========================================================================

% 3.4.1. GO parameters
% -------------------------------------------------------------------------

switch lower(optimScope)
  case {'go','all'}
    
    GOCat                           = goDesMat;
    
    % Correct GO parameters for:
    % - absence of variation across conditions
    if ~signatureGO(3)
      GOCat(end,:)                  = [];
    end
end

switch lower(optimScope)    
  case {'stop','all'}
    if ~isequal(signatureGO(1:2),[0 0]') && ...
       ~isequal(signatureSTOP(1:2),[0 0]')
      % GO parameters vary across stimuli and/or responses
      % STOP parameters vary across stimuli and/or responses
      
      % Correct GO parameters for:
      % - SSDs
      % - STOP parameter levels
      % - absence of variation across conditions
      
      GOCatStopTrial                              = stopDesMat;
      GOCatStopTrial([1,nFactGOStmRsp+2:end-1],:) = [];
      if ~signatureGO(3)
        GOCatStopTrial(end,:)                     = [];
      end
      
    elseif ~isequal(signatureGO(1:2),[0 0]') && ...
            isequal(signatureSTOP(1:2),[0 0]')  
      % GO parameters vary across stimuli and/or responses
      % STOP parameters do not vary across stimuli and/or responses
      
      % Correct GO parameters for:
      % - SSDs
      % - absence of variation across conditions
      
      GOCatStopTrial                              = stopDesMat;
      GOCatStopTrial(1,:)                         = [];
      if ~signatureGO(3)
        GOCatStopTrial(end,:)                     = [];
      end
      
    elseif isequal(signatureGO(1:2),[0 0]') && ...
          ~isequal(signatureSTOP(1:2),[0 0]')  
      % GO parameters do not vary across stimuli and/or responses
      % STOP parameters vary across stimuli and/or responses
      
      % Correct GO parameters for:
      % - SSDs
      % - absence of variation across conditions
      
      GOCatStopTrial                              = stopDesMat;
      GOCatStopTrial(1:nFactSTOPStmRsp+1,:)       = [];
      if ~signatureGO(3)
        GOCatStopTrial(end,:)                     = [];
      end
      
    elseif isequal(signatureGO(1:2),[0 0]') && ...
           isequal(signatureSTOP(1:2),[0 0]')  
      % GO parameters do not vary across stimuli and/or responses
      % STOP parameters do not vary across stimuli and/or responses
      
      % Correct GO parameters for:
      % - SSDs
      % - absence of variation across conditions
      
      GOCatStopTrial                              = stopDesMat;
      GOCatStopTrial(1,:)                         = [];
      if ~signatureGO(3)
        GOCatStopTrial(end,:)                     = [];
      end
    end
end

% 3.3.2. STOP parameters
% -------------------------------------------------------------------------

switch lower(optimScope)
  case {'stop','all'}
    if ~isequal(signatureGO(1:2),[0 0]') && ...
       ~isequal(signatureSTOP(1:2),[0 0]')
      % GO parameters vary across stimuli and/or responses
      % STOP parameters vary across stimuli and/or responses
      
      % Correct STOP parameters for:
      % - SSDs
      % - GO parameter levels
      % - absence of variation across conditions
      
      STOPCat                            = stopDesMat;
      STOPCat(1:nFactGOStmRsp+1,:)       = [];
      if ~signatureSTOP(3)
        STOPCat(end,:)                   = [];
      end
      
    elseif ~isequal(signatureGO(1:2),[0 0]') && ...
            isequal(signatureSTOP(1:2),[0 0]')  
      % GO parameters vary across stimuli and/or responses
      % STOP parameters do not vary across stimuli and/or responses
      
      % Correct STOP parameters for:
      % - SSDs
      % - GO parameter levels
      % - absence of variation across conditions
      
      STOPCat                                     = stopDesMat;
      STOPCat(1:nFactGOStmRsp+1,:)                = [];
      if ~signatureSTOP(3)
        STOPCat(end,:)                            = [];
      end 
      
    elseif isequal(signatureGO(1:2),[0 0]') && ...
          ~isequal(signatureSTOP(1:2),[0 0]')  
      % GO parameters do not vary across stimuli and/or responses
      % STOP parameters vary across stimuli and/or responses
      
      % Correct STOP parameters for:
      % - SSDs
      % - absence of variation across conditions
      
      STOPCat                                     = stopDesMat;
      STOPCat(1,:)                                = [];
      if ~signatureSTOP(3)
        STOPCat(end,:)                            = [];
      end
      
    elseif isequal(signatureGO(1:2),[0 0]') && ...
           isequal(signatureSTOP(1:2),[0 0]')  
      % GO parameters do not vary across stimuli and/or responses
      % STOP parameters do not vary across stimuli and/or responses
      
      % Correct STOP parameters for:
      % - SSDs
      % - absence of variation across conditions
      
      STOPCat                                     = stopDesMat;
      STOPCat(1,:)                                = [];
      if ~signatureSTOP(3)
        STOPCat(end,:)                            = [];
      end
    end
end

% 3.4.3. Convert to cell arrays
% -------------------------------------------------------------------------

switch lower(optimScope)
  case {'go','all'}
    if numel(GOCat) == 1
      % If there is just one trial category
      GOCatCell               = {1};
    else
      GOCatCell               = mat2cell(GOCat, ...
                                         size(GOCat,1), ...
                                         ones(size(GOCat,2),1));
    end
end

switch lower(optimScope)
  case {'stop','all'}
    
    if numel(STOPCat) == 1
      % If there is just one trial category
      GOCatStopTrialCell      = {1};
      STOPCatCell             = {1};
    else
      GOCatStopTrialCell      = mat2cell(GOCatStopTrial, ...
                                         size(GOCatStopTrial,1), ...
                                         ones(size(GOCatStopTrial,2),1));
      STOPCatCell             = mat2cell(STOPCat, ...
                                         size(STOPCat,1), ...
                                         ones(size(STOPCat,2),1));
    end
end

% 3.5. Specify functions function handles for producing tags
% =========================================================================
% Trial tag
%   Example: goTrial_c1_{GO:r1c1}
%   Example: stopTrial_ssd1_c1_{GO:c1}_{STOP}
%
% Parameter tag
%   Example: {GO:r1c1}
%   Example: {STOP}

% 3.5.1. Trial function handles
% -------------------------------------------------------------------------

switch lower(optimScope)
  case {'go','all'}
    funTrialGo    = @(in1) sprintf('goTrial_c%d',in1);
end

switch lower(optimScope)
  case {'stop','all'}
    funTrialStop  = @(in1) sprintf('stopTrial_ssd%d_c%d',in1);
end

% 3.5.2. Parameter function handles
% -------------------------------------------------------------------------
if isequal(signatureGO,[0 0 0]')
  funGO       = @(in1) sprintf('{GO%c}',in1);
elseif isequal(signatureGO,[1 0 0]')
  funGO       = @(in1) sprintf('{GO:s%d}',in1);
elseif isequal(signatureGO,[0 1 0]')
  funGO       = @(in1) sprintf('{GO:r%d}',in1);
elseif isequal(signatureGO,[0 0 1]')
  funGO       = @(in1) sprintf('{GO:c%d}',in1);
  funGO       = @(in1) sprintf('{GO:r%d,c%d}',in1);
elseif isequal(signatureGO,[1 1 0]')
  funGO       = @(in1) sprintf('{GO:s%d,r%d}',in1);
elseif isequal(signatureGO,[1 0 1]')
  funGO       = @(in1) sprintf('{GO:s%d,c%d}',in1);
elseif isequal(signatureGO,[0 1 1]')
  funGO       = @(in1) sprintf('{GO:r%d,c%d}',in1);
elseif isequal(signatureGO,[1 1 1]')
  funGO       = @(in1) sprintf('{GO:s%d,r%d,c%d}',in1);
end

switch lower(optimScope)
  case {'stop','all'}
    if isequal(signatureSTOP,[0 0 0]')
      funSTOP       = @(in1) sprintf('{STOP%c}',in1);
    elseif isequal(signatureSTOP,[1 0 0]')
      funSTOP       = @(in1) sprintf('{STOP:s%d}',in1);
    elseif isequal(signatureSTOP,[0 1 0]')
      funSTOP       = @(in1) sprintf('{STOP:r%d}',in1);
    elseif isequal(signatureSTOP,[0 0 1]')
      funSTOP       = @(in1) sprintf('{STOP:c%d}',in1);
    elseif isequal(signatureSTOP,[1 1 0]')
      funSTOP       = @(in1) sprintf('{STOP:s%d,r%d}',in1);
    elseif isequal(signatureSTOP,[1 0 1]')
      funSTOP       = @(in1) sprintf('{STOP:s%d,c%d}',in1);
    elseif isequal(signatureSTOP,[0 1 1]')
      funSTOP       = @(in1) sprintf('{STOP:r%d,c%d}',in1);
    elseif isequal(signatureSTOP,[1 1 1]')
      funSTOP       = @(in1) sprintf('{STOP:s%d,r%d,c%d}',in1);
    end
end

% 3.6. Produce trial tags
% =========================================================================
% Trial tags specify condition, and for stop trials, stop-signal delay.
% Parameter tags specify all factors that are varied.

switch lower(optimScope)
  case 'go'
    tag         = cellfun(@(in1,in2) [funTrialGo(in1),'_',funGO(in2)], ...
                          goTrialCatCell, ...
                          GOCatCell, ...
                          'Uni',0)';
    
    nTrialCatGo = numel(tag);
    nTrialCat   = numel(tag);
  case 'stop'
    tag         = cellfun(@(in1,in2,in3) [funTrialStop(in1),'_',funGO(in2),'_',funSTOP(in3)], ...
                         stopTrialCatCell,GOCatStopTrialCell,STOPCatCell,'Uni',0)';
    
    nTrialCatGo = 0;
    nTrialCat   = numel(tag);
  case 'all'
    tagGo       = cellfun(@(in1,in2) [funTrialGo(in1),'_',funGO(in2)], ...
                          goTrialCatCell,GOCatCell,'Uni',0)';
    tagStop     = cellfun(@(in1,in2,in3) [funTrialStop(in1),'_',funGO(in2),'_',funSTOP(in3)], ...
                          stopTrialCatCell,GOCatStopTrialCell,STOPCatCell,'Uni',0)';
    
    tag         = [tagGo;tagStop];
    
    nTrialCatGo = numel(tagGo);;
    nTrialCat   = numel(tag);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. CLASSIFY TRIALS, COMPUTE DESCRIPTIVES, AND SAVE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4.1. Pre-allocate arrays for logging
% =====================================================================

% Dataset array
obs             = dataset({cell(nTrialCat,1),'trialCat'}, ...
                          {cell(nTrialCat,1),'onset'}, ...
                          {cell(nTrialCat,1),'duration'}, ...
                          {nan(nTrialCat,1),'ssd'}, ...
                          {nan(nTrialCat,1),'nTotal'}, ...
                          {nan(nTrialCat,1),'nGoCCorr'}, ...
                          {nan(nTrialCat,1),'nGoCError'}, ...
                          {nan(nTrialCat,1),'nStopICorr'}, ...
                          {nan(nTrialCat,1),'nStopIErrorCCorr'}, ...
                          {nan(nTrialCat,1),'nStopIErrorCError'}, ...
                          {nan(nTrialCat,1),'pTotal'}, ...
                          {nan(nTrialCat,1),'pGoCCorr'}, ...
                          {nan(nTrialCat,1),'pGoCError'}, ...
                          {nan(nTrialCat,1),'pStopICorr'}, ...
                          {nan(nTrialCat,1),'pStopIErrorCCorr'}, ...
                          {nan(nTrialCat,1),'pStopIErrorCError'}, ...
                          {cell(nTrialCat,1),'rtGoCCorr'}, ...
                          {cell(nTrialCat,1),'rtGoCError'}, ...
                          {cell(nTrialCat,1),'rtStopICorr'}, ...
                          {cell(nTrialCat,1),'rtStopIErrorCCorr'}, ...
                          {cell(nTrialCat,1),'rtStopIErrorCError'}, ...
                          {cell(nTrialCat,1),'rtQGoCCorr'}, ...
                          {cell(nTrialCat,1),'rtQGoCError'}, ...
                          {cell(nTrialCat,1),'rtQStopICorr'}, ...
                          {cell(nTrialCat,1),'rtQStopIErrorCCorr'}, ...
                          {cell(nTrialCat,1),'rtQStopIErrorCError'}, ...
                          {cell(nTrialCat,1),'cumProbGoCCorr'}, ...
                          {cell(nTrialCat,1),'cumProbGoCError'}, ...
                          {cell(nTrialCat,1),'cumProbStopICorr'}, ...
                          {cell(nTrialCat,1),'cumProbStopIErrorCCorr'}, ...
                          {cell(nTrialCat,1),'cumProbStopIErrorCError'}, ...
                          {cell(nTrialCat,1),'cumProbDefectiveGoCCorr'}, ...
                          {cell(nTrialCat,1),'cumProbDefectiveGoCError'}, ...
                          {cell(nTrialCat,1),'cumProbDefectiveStopICorr'}, ...
                          {cell(nTrialCat,1),'cumProbDefectiveStopIErrorCCorr'}, ...
                          {cell(nTrialCat,1),'cumProbDefectiveStopIErrorCError'}, ...
                          {cell(nTrialCat,1),'probMassGoCCorr'}, ...
                          {cell(nTrialCat,1),'probMassGoCError'}, ...
                          {cell(nTrialCat,1),'probMassStopICorr'}, ...
                          {cell(nTrialCat,1),'probMassStopIErrorCCorr'}, ...
                          {cell(nTrialCat,1),'probMassStopIErrorCError'}, ...
                          {cell(nTrialCat,1),'probMassDefectiveGoCCorr'}, ...
                          {cell(nTrialCat,1),'probMassDefectiveGoCError'}, ...
                          {cell(nTrialCat,1),'probMassDefectiveStopICorr'}, ...
                          {cell(nTrialCat,1),'probMassDefectiveStopIErrorCCorr'}, ...
                          {cell(nTrialCat,1),'probMassDefectiveStopIErrorCError'});                    

for iTrialCat = 1:nTrialCat

  obs.trialCat{iTrialCat} = tag{iTrialCat};

  % 4.2. Classify on the basis of trial (go/stop) and experimental factors
  % (stimulus, response, condition)
  % =======================================================================
  % Note that trials are always categorized by condition, irrespective of
  % whether GO parameters vary between conditions.
  
%   % If this is a Go trial
%   if ~isempty(regexp(tag{iTrialCat},'goTrial.*', 'once'))
% 
%     if isequal(signatureGO,[0 0 0]')
%       iSelect = find(data.stm2    == 0 & ...
%                      data.cnd     == goTrialCat(end,iTrialCat));
%     elseif isequal(signatureGO,[1 0 0]')
%       iSelect = find(data.stm2    == 0 & ...
%                      data.stm1    == GOCat(1,iTrialCat) & ...
%                      data.cnd     == goTrialCat(end,iTrialCat));
%     elseif isequal(signatureGO,[0 1 0]')
%       iSelect = find(data.stm2    == 0 & ...
%                      data.resp    == GOCat(1,iTrialCat) & ...
%                      data.cnd     == goTrialCat(end,iTrialCat));
%     elseif isequal(signatureGO,[0 0 1]')
%       iSelect = find(data.stm2    == 0  & ...
%                      data.cnd     == goTrialCat(end,iTrialCat));
%     elseif isequal(signatureGO,[1 1 0]')
%       iSelect = find(data.stm2    == 0 & ...
%                      data.stm1    == GOCat(1,iTrialCat) & ...
%                      data.resp    == GOCat(2,iTrialCat) & ...
%                      data.cnd     == goTrialCat(end,iTrialCat));
%     elseif isequal(signatureGO,[1 0 1]')
%       iSelect = find(data.stm2    == 0 & ...
%                      data.stm1    == GOCat(1,iTrialCat)  & ...
%                      data.cnd     == goTrialCat(end,iTrialCat));
%     elseif isequal(signatureGO,[0 1 1]')
%       iSelect = find(data.stm2    == 0 & ...
%                      data.resp    == GOCat(1,iTrialCat) & ...
%                      data.cnd     == goTrialCat(end,iTrialCat));
%     elseif isequal(signatureGO,[1 1 1]')
%       iSelect = find(data.stm2    == 0 & ...
%                      data.stm1    == GOCat(1,iTrialCat) & ...
%                      data.resp    == GOCat(2,iTrialCat) & ...
%                      data.cnd     == goTrialCat(end,iTrialCat));
%     end
  % If this is a Go trial
  if ~isempty(regexp(tag{iTrialCat},'goTrial.*', 'once'))

    if isequal(signatureGO,[0 0 0]')
      iSelect = find(data.stm2    == 0 & ...
                     data.cnd     == goTrialCat(end,iTrialCat));
    elseif isequal(signatureGO,[1 0 0]')
      iSelect = find(data.stm2    == 0 & ...
                     data.stm1    == GOCat(1,iTrialCat) & ...
                     data.cnd     == goTrialCat(end,iTrialCat));
    elseif isequal(signatureGO,[0 1 0]')
      iSelect = find(data.stm2    == 0 & ...
                     data.rsp1    == GOCat(1,iTrialCat) & ...
                     data.cnd     == goTrialCat(end,iTrialCat));
    elseif isequal(signatureGO,[0 0 1]')
      iSelect = find(data.stm2    == 0  & ...
                     data.cnd     == goTrialCat(end,iTrialCat));
    elseif isequal(signatureGO,[1 1 0]')
      iSelect = find(data.stm2    == 0 & ...
                     data.stm1    == GOCat(1,iTrialCat) & ...
                     data.rsp1    == GOCat(2,iTrialCat) & ...
                     data.cnd     == goTrialCat(end,iTrialCat));
    elseif isequal(signatureGO,[1 0 1]')
      iSelect = find(data.stm2    == 0 & ...
                     data.stm1    == GOCat(1,iTrialCat)  & ...
                     data.cnd     == goTrialCat(end,iTrialCat));
    elseif isequal(signatureGO,[0 1 1]')
      iSelect = find(data.stm2    == 0 & ...
                     data.rsp1    == GOCat(1,iTrialCat) & ...
                     data.cnd     == goTrialCat(end,iTrialCat));
    elseif isequal(signatureGO,[1 1 1]')
      iSelect = find(data.stm2    == 0 & ...
                     data.stm1    == GOCat(1,iTrialCat) & ...
                     data.rsp1    == GOCat(2,iTrialCat) & ...
                     data.cnd     == goTrialCat(end,iTrialCat));
    end

  % If this is a Stop trial 
  elseif ~isempty(regexp(tag{iTrialCat},'stopTrial.*', 'once'))

    % Correspondence of indexing between iTrialCat and the columns of the
    % matrices GOCatStopTrial, stopTrialCat, and STOPCat depends on 
    % optimization scope.
    
    switch lower(optimScope)
      case 'stop'
        iCol = iTrialCat;
      case 'all'
        iCol = iTrialCat - nTrialCatGo;
    end
    
%     % Select trials based on GO criteria
%     if isequal(signatureGO,[0 0 0]')
%       iSelectGo = find(data.subj);
%     elseif isequal(signatureGO,[1 0 0]')
%       iSelectGo = find(data.stm1    == GOCatStopTrial(1,iCol));
%     elseif isequal(signatureGO,[0 1 0]')
%       iSelectGo = find(data.resp    == GOCatStopTrial(1,iCol));
%     elseif isequal(signatureGO,[0 0 1]')
%       iSelectGo = find(data.cnd     == GOCatStopTrial(1,iCol));
%     elseif isequal(signatureGO,[1 1 0]')
%       iSelectGo = find(data.stm1    == GOCatStopTrial(1,iCol) & ...
%                        data.resp    == GOCatStopTrial(2,iCol));
%     elseif isequal(signatureGO,[1 0 1]')
%       iSelectGo = find(data.stm1    == GOCatStopTrial(1,iCol) & ...
%                        data.cnd     == GOCatStopTrial(2,iCol));
%     elseif isequal(signatureGO,[0 1 1]')
%       iSelectGo = find(data.resp    == GOCatStopTrial(1,iCol) & ...
%                        data.cnd     == GOCatStopTrial(2,iCol));
%     elseif isequal(signatureGO,[1 1 1]')
%       iSelectGo = find(data.stm1    == GOCatStopTrial(1,iCol) & ...
%                        data.resp    == GOCatStopTrial(2,iCol) & ...
%                        data.cnd     == GOCatStopTrial(3,iCol));
%     end
    % Select trials based on GO criteria
    if isequal(signatureGO,[0 0 0]')
      iSelectGo = find(data.subj);
    elseif isequal(signatureGO,[1 0 0]')
      iSelectGo = find(data.stm1    == GOCatStopTrial(1,iCol));
    elseif isequal(signatureGO,[0 1 0]')
      iSelectGo = find(data.rsp1    == GOCatStopTrial(1,iCol));
    elseif isequal(signatureGO,[0 0 1]')
      iSelectGo = find(data.cnd     == GOCatStopTrial(1,iCol));
    elseif isequal(signatureGO,[1 1 0]')
      iSelectGo = find(data.stm1    == GOCatStopTrial(1,iCol) & ...
                       data.rsp1    == GOCatStopTrial(2,iCol));
    elseif isequal(signatureGO,[1 0 1]')
      iSelectGo = find(data.stm1    == GOCatStopTrial(1,iCol) & ...
                       data.cnd     == GOCatStopTrial(2,iCol));
    elseif isequal(signatureGO,[0 1 1]')
      iSelectGo = find(data.rsp1    == GOCatStopTrial(1,iCol) & ...
                       data.cnd     == GOCatStopTrial(2,iCol));
    elseif isequal(signatureGO,[1 1 1]')
      iSelectGo = find(data.stm1    == GOCatStopTrial(1,iCol) & ...
                       data.rsp1    == GOCatStopTrial(2,iCol) & ...
                       data.cnd     == GOCatStopTrial(3,iCol));
    end

    % Select trials based on STOP criteria
    % N.B. For Go trials, with overt responses, selection for different
    % responses is done based on resp. For Stop trials, this is not
    % possible, because there is no overt response of the Stop/Ignore
    % process. Therefore, selection happens based on rsp2, i.e. the
    % accumulator that SHOULD have reached threshold first. Obviously, this
    % needs to be changed in the future.
    
    if isequal(signatureSTOP,[0 0 0]')
      iSelectStop = find(data.stm2    > 0 & ...
                         data.iSSD    == stopTrialCat(1,iCol) & ...
                         data.cnd     == stopTrialCat(end,iCol));
    elseif isequal(signatureSTOP,[1 0 0]')
      iSelectStop = find(data.stm2    == 1 & ...
                         data.iSSD    == stopTrialCat(1,iCol) & ...
                         data.stm2    == STOPCat(1,iCol) & ...
                         data.cnd     == stopTrialCat(end,iCol));
    elseif isequal(signatureSTOP,[0 1 0]')
      iSelectStop = find(data.stm2    > 0 & ...
                         data.iSSD    == stopTrialCat(1,iCol) & ...
                         data.rsp2    == STOPCat(1,iCol) & ...
                         data.cnd     == stopTrialCat(end,iCol));
    elseif isequal(signatureSTOP,[0 0 1]')
      iSelectStop = find(data.stm2    > 0 & ...
                         data.iSSD    == stopTrialCat(1,iCol) & ...
                         data.cnd     == stopTrialCat(end,iCol));
    elseif isequal(signatureSTOP,[1 1 0]')
      iSelectStop = find(data.iSSD    == stopTrialCat(1,iCol) & ...
                         data.stm2    == STOPCat(1,iCol) & ...
                         data.rsp2    == STOPCat(2,iCol));
    elseif isequal(signatureSTOP,[1 0 1]')
      iSelectStop = find(data.iSSD    == stopTrialCat(1,iCol) & ...
                         data.stm2    == STOPCat(1,iCol) & ...
                         data.cnd     == stopTrialCat(end,iCol));
    elseif isequal(signatureSTOP,[0 1 1]')
      iSelectStop = find(data.stm2    > 0 & ...
                         data.iSSD    == stopTrialCat(1,iCol) & ...
                         data.rsp2    == STOPCat(1,iCol) & ...
                         data.cnd     == stopTrialCat(end,iCol));
    elseif isequal(signatureSTOP,[1 1 1]')
      iSelectStop = find(data.iSSD    == stopTrialCat(1,iCol) & ...
                         data.stm2    == STOPCat(1,iCol) & ...
                         data.rsp2    == STOPCat(2,iCol) & ...
                         data.cnd     == stopTrialCat(end,iCol));
    end

    % Only keep trials satisfying both criteria
    iSelect = intersect(iSelectGo,iSelectStop);

  end

  % 4.3. Narrow down the classification based on trial type
  % =======================================================================
  % For go trials, distinguish:
  % - correct choice
  % - error choice
  % For stop trials, distinguish:
  % - signal-inhibit
  % - signal-respond, correct choice
  % - signal-respond, error choice
  
  if ~isempty(regexp(tag{iTrialCat},'goTrial.*', 'once'))
      
      iGoCCorr              = intersect(iSelect,find(data.rsp1 == data.resp));
      iGoCError             = intersect(iSelect,find(data.rsp1 ~= data.resp));
      
  elseif ~isempty(regexp(tag{iTrialCat},'stopTrial.*', 'once'))
      
      iStopICorr            = intersect(iSelect,find(data.rt == 0));
      iStopIErrorCCorr      = intersect(iSelect,find(data.rt > 0 & data.rsp1 == data.resp));
      iStopIErrorCError     = intersect(iSelect,find(data.rt > 0 & data.rsp1 ~= data.resp));
  end
  
  % Number of trials
  % -----------------------------------------------------------------------
  if ~isempty(regexp(tag{iTrialCat},'goTrial.*', 'once'))
      % Compute
      nGoCCorr                          = numel(find(iGoCCorr));
      nGoCError                         = numel(find(iGoCError));
      nGoTotal                          = nGoCCorr + nGoCError;
      
      % Log
      obs.nTotal(iTrialCat)             = nGoTotal;
      obs.nGoCCorr(iTrialCat)           = nGoCCorr;
      obs.nGoCError(iTrialCat)          = nGoCError;
      
  elseif ~isempty(regexp(tag{iTrialCat},'stopTrial.*', 'once'))
      % Compute
      nStopICorr                        = numel(find(iStopICorr));
      nStopIErrorCCorr                  = numel(find(iStopIErrorCCorr));
      nStopIErrorCError                 = numel(find(iStopIErrorCError));
      nStopTotal                        = nStopICorr + nStopIErrorCCorr + nStopIErrorCError;
      
      % Log
      obs.nTotal(iTrialCat)             = nStopTotal;
      obs.nStopICorr(iTrialCat)         = nStopICorr;
      obs.nStopIErrorCCorr(iTrialCat)   = nStopIErrorCCorr;
      obs.nStopIErrorCError(iTrialCat)  = nStopIErrorCError;
  end
  
  % Trial probabilities
  % -----------------------------------------------------------------------
  if ~isempty(regexp(tag{iTrialCat},'goTrial.*', 'once'))
      % Compute
      if nGoTotal > 0
          pGoCCorr                          = nGoCCorr/nGoTotal;
          pGoCError                         = nGoCError/nGoTotal;
      else
          pGoCCorr                          = 0;
          pGoCError                         = 0;
      end
      
      % Log
      obs.pTotal(iTrialCat)             = pGoCCorr + pGoCError;
      obs.pGoCCorr(iTrialCat)           = pGoCCorr;
      obs.pGoCError(iTrialCat)          = pGoCError;
      
  elseif ~isempty(regexp(tag{iTrialCat},'stopTrial.*', 'once'))
      % Compute
      if nStopTotal > 0
          pStopICorr                        = nStopICorr/nStopTotal;
          pStopIErrorCCorr                  = nStopIErrorCCorr/nStopTotal;
          pStopIErrorCError                 = nStopIErrorCError/nStopTotal;
      else
          pStopICorr                        = 0;
          pStopIErrorCCorr                  = 0;
          pStopIErrorCError                 = 0;
      end
      
      % Log
      obs.pTotal(iTrialCat)             = pStopICorr + pStopIErrorCCorr + pStopIErrorCError;
      obs.pStopICorr(iTrialCat)         = pStopICorr;
      obs.pStopIErrorCCorr(iTrialCat)   = pStopIErrorCCorr;
      obs.pStopIErrorCError(iTrialCat)  = pStopIErrorCError;
      
  end
  
  % Response time
  % -----------------------------------------------------------------------
  
  if ~isempty(regexp(tag{iTrialCat},'goTrial.*', 'once'))
    if nGoCCorr > 0
        % Compute
        rtGoCCorr = sort(data.rt(iGoCCorr));
        
        % Log
        obs.rtGoCCorr{iTrialCat} = rtGoCCorr;
    end
    
    if nGoCError > 0
        % Compute
        rtGoCError = sort(data.rt(iGoCError));
        
        % Note: stimOns for iTargetGO and iNonTargetGO are the same; I use 
        % iTargetGO instead of iNonTargetGO because iTargetGO is always a 
        % scalar, iNonTargetGO not.
      
        % Log
        obs.rtGoCError{iTrialCat} = rtGoCError;
    end
    
  elseif ~isempty(regexp(tag{iTrialCat},'stopTrial.*', 'once'))
    if nStopICorr > 0
        % Compute
        rtStopICorr = sort(data.rt(iStopICorr));
        
        % Log
        obs.rtStopICorr{iTrialCat} = rtStopICorr;
    end
    
    if nStopIErrorCCorr > 0
        % Compute
        rtStopIErrorCCorr = sort(data.rt(iStopIErrorCCorr));
        
        % Log
        obs.rtStopIErrorCCorr{iTrialCat} = rtStopIErrorCCorr;
    end
    
    if nStopIErrorCError > 0
        % Compute
        rtStopIErrorCError = sort(data.rt(iStopIErrorCError));
      
        % Log
        obs.rtStopIErrorCError{iTrialCat} = rtStopIErrorCError;
    end
    
  end
  
  % Stop-signal delay
  % -----------------------------------------------------------------------
  if ~isempty(regexp(tag{iTrialCat},'stopTrial.*', 'once'))
    if numel(unique(nonnans(data.ssd(iSelect)))) > 1
      error('More than one SSD detected');
    end
    obs.ssd(iTrialCat)      = max([0,unique(nonnans(data.ssd(iSelect)))]);
  end

  % RT bin data
  % -----------------------------------------------------------------------
  if ~isempty(regexp(tag{iTrialCat},'goTrial.*', 'once'))
      
      if nGoCCorr > 0
          [obs.rtQGoCCorr{iTrialCat}, ...
           obs.cumProbGoCCorr{iTrialCat}, ...   
           obs.cumProbDefectiveGoCCorr{iTrialCat}, ...
           obs.probMassGoCCorr{iTrialCat}, ...
           obs.probMassDefectiveGoCCorr{iTrialCat}] = ...
           sam_bin_data(rtGoCCorr,pGoCCorr,cumProb,minBinSize,dt);
      end
      
      if nGoCError > 0
          [obs.rtQGoCError{iTrialCat}, ...
           obs.cumProbGoCError{iTrialCat}, ...
           obs.cumProbDefectiveGoCError{iTrialCat}, ...
           obs.probMassGoCError{iTrialCat}, ...
           obs.probMassDefectiveGoCError{iTrialCat}] = ...
           sam_bin_data(rtGoCError,pGoCError,cumProb,minBinSize,dt); 
      end
      
  elseif ~isempty(regexp(tag{iTrialCat},'stopTrial.*', 'once'))
      
      if nStopICorr > 0
          
          obs.probMassStopICorr{iTrialCat} = 1;
          obs.probMassDefectiveStopICorr{iTrialCat} = pStopICorr;
          
      end
      
      if nStopIErrorCCorr > 0
          [obs.rtQStopIErrorCCorr{iTrialCat}, ...
           obs.cumProbStopIErrorCCorr{iTrialCat}, ...
           obs.cumProbDefectiveStopIErrorCCorr{iTrialCat}, ...
           obs.probMassStopIErrorCCorr{iTrialCat}, ...
           obs.probMassDefectiveStopIErrorCCorr{iTrialCat}] = ...
           sam_bin_data(rtStopIErrorCCorr,pStopIErrorCCorr,cumProb,minBinSize,dt);
      end
      
      if nStopIErrorCError > 0
          [obs.rtQStopIErrorCError{iTrialCat}, ...
           obs.cumProbStopIErrorCError{iTrialCat}, ...
           obs.cumProbDefectiveStopIErrorCError{iTrialCat}, ...
           obs.probMassStopIErrorCError{iTrialCat}, ...
           obs.probMassDefectiveStopIErrorCError{iTrialCat}] = ...
           sam_bin_data(rtStopIErrorCError,pStopIErrorCError,cumProb,minBinSize,dt);
      end
      
  end

  % Stimulus onsets and durations
  % -----------------------------------------------------------------------
  if ~isempty(regexp(tag{iTrialCat},'goTrial.*', 'once'))
    obs.onset{iTrialCat}    = blkdiag(trueNStm{:})*[stmOns(1) 0]';
    obs.duration{iTrialCat} = blkdiag(trueNStm{:})*[stmDur(1) 0]';
  elseif ~isempty(regexp(tag{iTrialCat},'stopTrial.*', 'once'))
    obs.onset{iTrialCat}    = blkdiag(trueNStm{:})*[stmOns(1) stmOns(1) + obs.ssd(iTrialCat)]';
    obs.duration{iTrialCat} = blkdiag(trueNStm{:})*(stmDur)';
  end
  
end
function obs = sam_categorize_data(SAM)
% SAM_PROCESS_RAW_DATA <Synopsis of what this function does> 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
%
% Here is a description of the data
% 
% Subjects performed 12 sessions
% In each session there were three conditions: 2-choice, 4-choice, 6-choice
%
% Conditions were presented in three blocks. The order of blocks within
% each session was randomized
%
% Each choice condition started with a practice block of 36 trials without
% stop-signals. This no-stop-signal block was followed by  another practice
% block of 36 trials with stop-signals. After the two practice blocks,
% there were two experimental blocks of 120 trials.
%
%
%
%
%
%
%
%
% SYNTAX 
% SAM_PROCESS_RAW_DATA; 
%  
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Fri 23 Aug 2013 12:08:07 CDT by bram 
% $Modified: Fri 23 Aug 2013 12:08:07 CDT by bram 
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESSING INPUTS AND SPECIFYING VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1.1. Process inputs
% =========================================================================

rawDataDir        = SAM.io.rawDataDir;
workDir           = SAM.io.workDir;

nStm              = SAM.expt.nStm;                                            
nRsp              = SAM.expt.nRsp;
nCnd              = SAM.expt.nCnd;                                            
nSsd              = SAM.expt.nSsd;                                            
trialDur          = SAM.expt.trialDur;
stmOns            = SAM.expt.stmOns;                                        
stmDur            = SAM.expt.stmDur;                                       

modelToFit        = SAM.model.variants.toFit;

qntls             = SAM.optim.cost.stat.cumProb;
minBinSize        = SAM.optim.cost.stat.minBinSize;

% 1.2. Dynamic variables
% =========================================================================

% Maximum number of stimuli and response, across conditions
maxNStm           = max(cell2mat(nStm(:)),[],1);
maxNRsp           = max(cell2mat(nRsp(:)),[],1);

taskFactors       = [maxNStm;maxNRsp;nCnd,nCnd];

% Miscellaneous
% -------------------------------------------------------------------------
trueNStm          = arrayfun(@(x) true(x,1),maxNStm,'Uni',0);

nClass            = numel(SAM.model.general.classNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. CATEGORIZE DATA BASED ON MODEL FEATURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Go trials
% =========================================================================================================================
signatureGo   = any(modelToFit.features(:,:,1),2);

if ~isequal(signatureGo,[0 0 0]')
  combiGo       = fullfact(taskFactors(signatureGo,1))';
else
  combiGo       = 1;
end

if isequal(signatureGo,[0 0 0]')
  funGO = @(in1) sprintf('{GO}',in1);
elseif isequal(signatureGo,[1 0 0]')
  funGO = @(in1) sprintf('{GO:s%d}',in1);
elseif isequal(signatureGo,[0 1 0]')
  funGO = @(in1) sprintf('{GO:r%d}',in1);
elseif isequal(signatureGo,[0 0 1]')
  funGO = @(in1) sprintf('{GO:c%d}',in1);
elseif isequal(signatureGo,[1 1 0]')
  funGO = @(in1) sprintf('{GO:s%d,r%d}',in1);
elseif isequal(signatureGo,[1 0 1]')
  funGO = @(in1) sprintf('{GO:s%d,c%d}',in1);
elseif isequal(signatureGo,[0 1 1]')
  funGO = @(in1) sprintf('{GO:r%d,c%d}',in1);
elseif isequal(signatureGo,[1 1 1]')
  funGO = @(in1) sprintf('{GO:s%d,r%d,c%d}',in1);
end

if all(combiGo(:) == 1)
  combiCellGo = mat2cell(combiGo,1,ones(size(combiGo,2),1));
else
  combiCellGo = mat2cell(combiGo,size(combiGo,1),ones(size(combiGo,2),1));
end

tagGo = cellfun(@(in1) ['goTrial_',funGO(in1)],combiCellGo,'Uni',0);

% Stop trials
% =========================================================================================================================
signatureStop         = any(modelToFit.features(:,:,2),2);

if ~isequal(signatureGo,[0 0 0]') && ~isequal(signatureStop,[0 0 0]')
  combiStop             = fullfact([nSsd;taskFactors(signatureGo,1);taskFactors(signatureStop,2)])';
elseif ~isequal(signatureGo,[0 0 0]') && isequal(signatureStop,[0 0 0]')  
  combiStop             = fullfact([nSsd;taskFactors(signatureGo,1)])';
elseif isequal(signatureGo,[0 0 0]') && ~isequal(signatureStop,[0 0 0]')  
  combiStop             = fullfact([nSsd;taskFactors(signatureStop,1)])';
elseif isequal(signatureGo,[0 0 0]') && isequal(signatureStop,[0 0 0]')  
  combiStop             = fullfact([nSsd,1,1])';
end

if isequal(signatureStop,[0 0 0]')
  funSTOP = @(in1) sprintf('{STOP}',in1);
elseif isequal(signatureStop,[1 0 0]')
  funSTOP = @(in1) sprintf('{STOP:s%d}',in1);
elseif isequal(signatureStop,[0 1 0]')
  funSTOP = @(in1) sprintf('{STOP:r%d}',in1);
elseif isequal(signatureStop,[0 0 1]')
  funSTOP = @(in1) sprintf('{STOP:c%d}',in1);
elseif isequal(signatureStop,[1 1 0]')
  funSTOP = @(in1) sprintf('{STOP:s%d,r%d}',in1);
elseif isequal(signatureStop,[1 0 1]')
  funSTOP = @(in1) sprintf('{STOP:s%d,c%d}',in1);
elseif isequal(signatureStop,[0 1 1]')
  funSTOP = @(in1) sprintf('{STOP:r%d,c%d}',in1);
elseif isequal(signatureStop,[1 1 1]')
  funSTOP = @(in1) sprintf('{STOP:s%d,r%d,c%d}',in1);
end

nFactGo = numel(taskFactors(signatureGo,1));
nFactStop = numel(taskFactors(signatureStop,1));

if all(all(combiStop(2:end,:) == 1))
  combiCellStop = mat2cell(combiStop,[1;1;1],ones(size(combiStop,2),1));
else
  combiCellStop = mat2cell(combiStop,[1;nFactGo;nFactStop],ones(size(combiStop,2),1));
end

tagStop = cellfun(@(in1,in2,in3) ['stopTrial_{ssd',sprintf('%d',in1),'}_',funGO(in2),'_',funSTOP(in3)],combiCellStop(1,:),combiCellStop(2,:),combiCellStop(3,:),'Uni',0);

% All tags
tagAll = [tagGo,tagStop]';

% All combi cells combines
combiCellAll = [combiCellGo,mat2cell(combiCellStop,size(combiCellStop,1),ones(1,size(combiCellStop,2)))];

% Number of trial categories
nTrialCat = numel([tagGo,tagStop]);
nTrialCatGo = numel(tagGo);
nTrialCatStop = numel(tagStop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. CLASSIFY TRIALS, COMPUTE DESCRIPTIVES, AND SAVE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subj = unique(allObs.subj);
nSubj = length(subj);

for iSubj = 1:nSubj
 
  % 4.1. Pre-allocate arrays for logging
  % =====================================================================
  
  % Dataset array
	obs           = dataset({cell(nTrialCat,1),'trialCat'}, ...
                          {cell(nTrialCat,1),'funGO'}, ...
                          {cell(nTrialCat,1),'funSTOP'}, ...
                          {cell(nTrialCat,1),'onset'}, ...
                          {cell(nTrialCat,1),'duration'}, ...
                          {zeros(nTrialCat,1),'ssd'}, ...
                          {zeros(nTrialCat,1),'nTotal'}, ...
                          {zeros(nTrialCat,1),'nCorr'}, ...
                          {zeros(nTrialCat,1),'nError'}, ...
                          {zeros(nTrialCat,1),'pTotal'}, ...
                          {zeros(nTrialCat,1),'pCorr'}, ...
                          {zeros(nTrialCat,1),'pError'}, ...
                          {cell(nTrialCat,1),'rtCorr'}, ...
                          {cell(nTrialCat,1),'rtError'}, ...
                          {cell(nTrialCat,1),'rtQCorr'}, ...
                          {cell(nTrialCat,1),'rtQError'}, ...
                          {cell(nTrialCat,1),'fCorr'}, ...
                          {cell(nTrialCat,1),'fError'}, ...
                          {cell(nTrialCat,1),'pMassCorr'}, ...
                          {cell(nTrialCat,1),'pMassError'}, ...
                          {cell(nTrialCat,1),'pDefectiveCorr'}, ...
                          {cell(nTrialCat,1),'pDefectiveError'});
                        
	% Some general matrices
  iSsd = [zeros(numel(tagGo),1);cell2mat(combiCellStop(1,:))'];
                        
	for iTrialCat = 1:nTrialCat
    
    obs.trialCat{iTrialCat} = tagAll{iTrialCat};
    
    % If this is a Go trial
    if ~isempty(regexp(tagAll{iTrialCat},'goTrial.*'))
      
      if isequal(signatureGo,[0 0 0]')
        iSelect = find(allObs.subj    == subj(iSubj) & ...
                       allObs.stm2    == 0);
      elseif isequal(signatureGo,[1 0 0]')
        iSelect = find(allObs.subj    == subj(iSubj) & ...
                       allObs.stm2    == 0 & ...
                       allObs.stm1    == combiGo(1,iTrialCat));
      elseif isequal(signatureGo,[0 1 0]')
        iSelect = find(allObs.subj    == subj(iSubj) & ...
                       allObs.stm2    == 0 & ...
                       allObs.rsp1    == combiGo(1,iTrialCat));
      elseif isequal(signatureGo,[0 0 1]')
        iSelect = find(allObs.subj    == subj(iSubj) & ...
                       allObs.stm2    == 0 & ...
                       allObs.cnd     == combiGo(1,iTrialCat));
      elseif isequal(signatureGo,[1 1 0]')
        iSelect = find(allObs.subj    == subj(iSubj) & ...
                       allObs.stm2    == 0 & ...
                       allObs.stm1    == combiGo(1,iTrialCat) & ...
                       allObs.rsp1    == combiGo(2,iTrialCat));
      elseif isequal(signatureGo,[1 0 1]')
        iSelect = find(allObs.subj    == subj(iSubj) & ...
                       allObs.stm2    == 0 & ...
                       allObs.stm1    == combiGo(1,iTrialCat) & ...
                       allObs.cnd     == combiGo(2,iTrialCat));
      elseif isequal(signatureGo,[0 1 1]')
        iSelect = find(allObs.subj    == subj(iSubj) & ...
                       allObs.stm2    == 0 & ...
                       allObs.rsp1    == combiGo(1,iTrialCat) & ...
                       allObs.cnd     == combiGo(2,iTrialCat));
      elseif isequal(signatureGo,[1 1 1]')
        iSelect = find(allObs.subj    == subj(iSubj) & ...
                       allObs.stm2    == 0 & ...
                       allObs.stm1    == combiGo(1,iTrialCat) & ...
                       allObs.rsp1    == combiGo(2,iTrialCat) & ...
                       allObs.cnd     == combiGo(3,iTrialCat));
      end
      
      
    % If this is a Stop trial 
    elseif ~isempty(regexp(tagAll{iTrialCat},'stopTrial.*'))
      
      % Select trials based on GO criteria
      if isequal(signatureGo,[0 0 0]')
        iSelectGo = find(allObs.subj    == subj(iSubj));
      elseif isequal(signatureGo,[1 0 0]')
        iSelectGo = find(allObs.subj    == subj(iSubj) & ...
                         allObs.stm1    == combiCellStop{2,iTrialCat-nTrialCatGo}(1));
      elseif isequal(signatureGo,[0 1 0]')
        iSelectGo = find(allObs.subj    == subj(iSubj) & ...
                         allObs.rsp1    == combiCellStop{2,iTrialCat-nTrialCatGo}(1));
      elseif isequal(signatureGo,[0 0 1]')
        iSelectGo = find(allObs.subj    == subj(iSubj) & ...
                         allObs.cnd     == combiCellStop{2,iTrialCat-nTrialCatGo}(1));
      elseif isequal(signatureGo,[1 1 0]')
        iSelectGo = find(allObs.subj    == subj(iSubj) & ...
                         allObs.stm1    == combiCellStop{2,iTrialCat-nTrialCatGo}(1) & ...
                         allObs.rsp1    == combiCellStop{2,iTrialCat-nTrialCatGo}(2));
      elseif isequal(signatureGo,[1 0 1]')
        iSelectGo = find(allObs.subj    == subj(iSubj) & ...
                         allObs.stm1    == combiCellStop{2,iTrialCat-nTrialCatGo}(1) & ...
                         allObs.cnd     == combiCellStop{2,iTrialCat-nTrialCatGo}(2));
      elseif isequal(signatureGo,[0 1 1]')
        iSelectGo = find(allObs.subj    == subj(iSubj) & ...
                         allObs.rsp1    == combiCellStop{2,iTrialCat-nTrialCatGo}(1)& ...
                         allObs.cnd     == combiCellStop{2,iTrialCat-nTrialCatGo}(2));
      elseif isequal(signatureGo,[1 1 1]')
        iSelectGo = find(allObs.subj    == subj(iSubj) & ...
                         allObs.stm1    == combiCellStop{2,iTrialCat-nTrialCatGo}(1) & ...
                         allObs.rsp1    == combiCellStop{2,iTrialCat-nTrialCatGo}(2) & ...
                         allObs.cnd     == combiCellStop{2,iTrialCat-nTrialCatGo}(3));
      end
      
      % Select trials based on STOP criteria
      if isequal(signatureStop,[0 0 0]')
        iSelectStop = find(allObs.subj    == subj(iSubj) & ...
                           allObs.stm2    > 0 & ...
                           allObs.iSSD    == combiCellStop{1,iTrialCat-nTrialCatGo}(1));
      elseif isequal(signatureStop,[1 0 0]')
        iSelectStop = find(allObs.subj    == subj(iSubj) & ...
                           allObs.stm2    == 1 & ...
                           allObs.iSSD    == combiCellStop{1,iTrialCat-nTrialCatGo}(1) & ...
                           allObs.stm2    == combiCellStop{3,iTrialCat-nTrialCatGo}(1));
      elseif isequal(signatureStop,[0 1 0]')
        iSelectStop = find(allObs.subj    == subj(iSubj) & ...
                           allObs.stm2    > 0 & ...
                           allObs.iSSD    == combiCellStop{1,iTrialCat-nTrialCatGo}(1) & ...
                           allObs.rsp2    == combiCellStop{3,iTrialCat-nTrialCatGo}(1));
      elseif isequal(signatureStop,[0 0 1]')
        iSelectStop = find(allObs.subj    == subj(iSubj) & ...
                           allObs.stm2    > 0 & ...
                           allObs.iSSD    == combiCellStop{1,iTrialCat-nTrialCatGo}(1) & ...
                           allObs.cnd     == combiCellStop{3,iTrialCat-nTrialCatGo}(1));
      elseif isequal(signatureStop,[1 1 0]')
        iSelectStop = find(allObs.subj    == subj(iSubj) & ...
                           allObs.iSSD    == combiCellStop{1,iTrialCat-nTrialCatGo}(1) & ...
                           allObs.stm2    == combiCellStop{3,iTrialCat-nTrialCatGo}(1) & ...
                           allObs.rsp2    == combiCellStop{3,iTrialCat-nTrialCatGo}(2));
      elseif isequal(signatureStop,[1 0 1]')
        iSelectStop = find(allObs.subj    == subj(iSubj) & ...
                           allObs.iSSD    == combiCellStop{1,iTrialCat-nTrialCatGo}(1) & ...
                           allObs.stm2    == combiCellStop{3,iTrialCat-nTrialCatGo}(1) & ...
                           allObs.cnd     == combiCellStop{3,iTrialCat-nTrialCatGo}(2));
      elseif isequal(signatureStop,[0 1 1]')
        iSelectStop = find(allObs.subj    == subj(iSubj) & ...
                           allObs.stm2    > 0 & ...
                           allObs.iSSD    == combiCellStop{1,iTrialCat-nTrialCatGo}(1) & ...
                           allObs.rsp2    == combiCellStop{3,iTrialCat-nTrialCatGo}(1) & ...
                           allObs.cnd     == combiCellStop{3,iTrialCat-nTrialCatGo}(2));
      elseif isequal(signatureStop,[1 1 1]')
        iSelectStop = find(allObs.subj    == subj(iSubj) & ...
                           allObs.iSSD    == combiCellStop{1,iTrialCat-nTrialCatGo}(1) & ...
                           allObs.stm2    == combiCellStop{3,iTrialCat-nTrialCatGo}(1) & ...
                           allObs.rsp2    == combiCellStop{3,iTrialCat-nTrialCatGo}(2) & ...
                           allObs.cnd     == combiCellStop{3,iTrialCat-nTrialCatGo}(3));
      end
      
      % Only keep trials satisfying both criteria
      iSelect = intersect(iSelectGo,iSelectStop);
      
    end
    
    iSelectCorr   = intersect(iSelect, find(allObs.acc == 2));
    iSelectError  = intersect(iSelect,find(allObs.acc ~= 2));
        
    % Number of trials
    obs.nTotal(iTrialCat)   = numel(iSelect);
    obs.nCorr(iTrialCat)    = numel(iSelectCorr);
    obs.nError(iTrialCat)   = numel(iSelectError);
    
    % Probability
    obs.pCorr(iTrialCat)    = numel(iSelectCorr)./numel(iSelect);
    obs.pError(iTrialCat)   = numel(iSelectError)./numel(iSelect);
    
    % Response time
    obs.rtCorr{iTrialCat}   = sort(allObs.rt(iSelectCorr));
    obs.rtError{iTrialCat}  = sort(allObs.rt(iSelectError));
    
    obs.ssd(iTrialCat)      = max([0,unique(nonnans(allObs.ssd(iSelect)))]);
    
    
    if ~isempty(regexp(tagAll{iTrialCat},'goTrial.*')) && obs.nCorr(iTrialCat) > 0
      [obs.rtQCorr{iTrialCat}, ...
       obs.pDefectiveCorr{iTrialCat}, ...
       obs.fCorr{iTrialCat}, ...
       obs.pMassCorr{iTrialCat}]  = sam_bin_data(obs.rtCorr{iTrialCat},obs.pCorr(iTrialCat),obs.nCorr(iTrialCat),qntls,minBinSize);
    end
    
    if obs.nError(iTrialCat) > 0
      [obs.rtQError{iTrialCat}, ...
       obs.pDefectiveError{iTrialCat}, ...
       obs.fError{iTrialCat}, ...
       obs.pMassError{iTrialCat}]  = sam_bin_data(obs.rtError{iTrialCat},obs.pError(iTrialCat),obs.nError(iTrialCat),qntls,minBinSize);
    end
    
    
    if ~isempty(regexp(tagAll{iTrialCat},'goTrial.*'))
      obs.onset{iTrialCat}    = blkdiag(trueNStm{:})*[stmOns(1) 0]';
      obs.duration{iTrialCat} = blkdiag(trueNStm{:})*[stmDur(1) 0]';
    elseif ~isempty(regexp(tagAll{iTrialCat},'stopTrial.*'))
      obs.onset{iTrialCat}    = blkdiag(trueNStm{:})*[stmOns(1) stmOns(1) + obs.ssd(iTrialCat)]';
      obs.duration{iTrialCat} = blkdiag(trueNStm{:})*[stmDur]';
    end
    
%     % 4.1.1. Compute timing diagram
%         % -----------------------------------------------------------------
%                                           % OUTPUT
%         [tStm, ...                        % - Time
%          uStm] ...                        % - Strength of stimulus (t)
%         = sam_spec_timing_diagram ...     % FUNCTION
%          ...                              % INPUT
%         (stimOns(:)', ...                 % - Stimulus onset time
%          stimDur(:)', ...                 % - Stimulus duration
%          [], ...                          % - Strength (default = 1);
%          0, ...                           % - Magnitude of extrinsic noise
%          dt, ...                          % - Time step
%          timeWindow);                     % - Time window
    
  end
    
%   % 4.5. Save data
%   % =======================================================================
%   fName = fullfile(outputDir,sprintf('data_preproc_subj%.2d.mat',subj(iSubj)));
%   save(fName, 'obs');
%   
%   fName = fullfile(outputDir,sprintf('all_trial_data_preproc.mat',subj(iSubj)));
%   save(fName, 'allObs');
end
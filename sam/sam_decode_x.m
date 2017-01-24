function [endoConn,extrMod,exoConn,intrMod,V,ETA,SE,SI,Z0,ZC,T0] = ...
          sam_decode_x(SAM,X,iTrial)
% function [endoConn,extrMod,exoConn,intrMod,V,SE,SI,Z0,ZC,accumOns] = ...
%           sam_decode_x(SAM,X,iTrial)
% SAM_DECODE_X <Synopsis of what this function does> 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% SAM_DECODE_X; 
%  
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Fri 23 Aug 2013 11:47:36 CDT by bram 
% $Modified: Fri 23 Aug 2013 11:47:36 CDT by bram 

 
% CONTENTS 
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS AND SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1.1. Process inputs
% ========================================================================= 

% Column indices
iZ0       = SAM.model.XCat.i.iZ0;
iZc       = SAM.model.XCat.i.iZc;
iV        = SAM.model.XCat.i.iV;
iVe       = SAM.model.XCat.i.iVe;
iEta      = SAM.model.XCat.i.iEta;
iT0       = SAM.model.XCat.i.iT0;
iSe       = SAM.model.XCat.i.iSe;
iSi       = SAM.model.XCat.i.iSi;
iK        = SAM.model.XCat.i.iK;
iWliw     = SAM.model.XCat.i.iWliw;
iWlib     = SAM.model.XCat.i.iWlib;
iWffiw    = SAM.model.XCat.i.iWffiw;

optimScope  = SAM.sim.scope;

nRsp      = SAM.expt.nRsp;
nStm      = SAM.expt.nStm;
nCnd      = SAM.expt.nCnd;

features  = SAM.model.variants.toFit.features;

switch lower(optimScope)
  case 'go'
    iCatClass = SAM.model.variants.toFit.XSpec.i.go.iCatClass;
  case 'stop'
    iCatClass = SAM.model.variants.toFit.XSpec.i.stop.iCatClass;
  case 'all'
    iCatClass = SAM.model.variants.toFit.XSpec.i.all.iCatClass;
end

trialCat    = SAM.optim.obs.trialCat{iTrial};
stmOns      = SAM.optim.obs.onset{iTrial};
stmDur      = SAM.optim.obs.duration{iTrial};

iTarget     = SAM.optim.modelMat.iTarget{iTrial};
iNonTarget  = SAM.optim.modelMat.iNonTarget{iTrial};

% #.#.#. Model matrices
% -------------------------------------------------------------------------------------------------------------------------
endoConn  = SAM.model.mat.endoConn;
exoConn   = SAM.model.mat.exoConn;

% 1.2. Specify dynamic variables
% ========================================================================= 
maxNStm           = max(cell2mat(nStm(:)),[],1);
maxNRsp           = max(cell2mat(nRsp(:)),[],1);

taskFactors       = [maxNStm;maxNRsp;nCnd,nCnd];

trueNStm          = arrayfun(@(x) true(x,1),maxNStm,'Uni',0);
trueNRsp          = arrayfun(@(x) true(x,1),maxNRsp,'Uni',0);

% Parse trial type 
if ~isempty(regexp(trialCat,'^goTrial_', 'once'))
  token = regexp(trialCat,'goTrial_c(\w*)_(\S*)','tokens');
  tagGO   = token{1}{2};
elseif ~isempty(regexp(trialCat,'^stopTrial_', 'once'))
  token = regexp(trialCat,'stopTrial_ssd(\w*)_c(\w*)_(\S*)_(\S*)','tokens');
%   tagSsd  = str2double(token{1}{1});
  tagGO   = token{1}{3};
  tagSTOP = token{1}{4};
end

% 1.2.1. Parameter indices per task factor
% -------------------------------------------------------------------------------------------------------------------------
% Specifies which specific index we need per task factor and accumulator category. A nan means no specific index, a number means a specific index.

indexMat = nan(3,2); 

if ~isempty(regexp(tagGO,'{GO}', 'once'))
elseif ~isempty(regexp(tagGO,'{GO:s+\d}', 'once'))
  indices = regexp(tagGO,'{GO:s(\d+)}','tokens');
  indexMat(1,1) = str2double(indices{1}{1});
elseif ~isempty(regexp(tagGO,'{GO:r+\d}', 'once'))
  indices = regexp(tagGO,'{GO:r(\d+)}','tokens');
  indexMat(2,1) = str2double(indices{1}{1});
elseif ~isempty(regexp(tagGO,'{GO:c+\d}', 'once'))
  indices = regexp(tagGO,'{GO:c(\d+)}','tokens');
  indexMat(3,1) = str2double(indices{1}{1});
elseif ~isempty(regexp(tagGO,'{GO:s+\d,r+\d}', 'once'))
  indices = regexp(tagGO,'{GO:s(\d+),r(\d+)','tokens');
  indexMat(1,1) = str2double(indices{1}{1});
  indexMat(2,1) = str2double(indices{1}{2});
elseif ~isempty(regexp(tagGO,'{GO:s+\d,c+\d}', 'once'))
  indices = regexp(tagGO,'{GO:s(\d+),c(\d+)}','tokens');
  indexMat(1,1) = str2double(indices{1}{1});
  indexMat(3,1) = str2double(indices{1}{2});
elseif ~isempty(regexp(tagGO,'{GO:r+\d,c+\d}', 'once'))
  indices = regexp(tagGO,'{GO:r(\d+),c(\d+)}','tokens');
  indexMat(2,1) = str2double(indices{1}{1});
  indexMat(3,1) = str2double(indices{1}{2});
elseif ~isempty(regexp(tagGO,'{GO:s+\d,r+\d,c+\d}', 'once'))
  indices = regexp(tagGO,'{GO:s(\d+),r(\d+),c(\d+)}','tokens');
  indexMat(1,1) = str2double(indices{1}{1});
  indexMat(2,1) = str2double(indices{1}{2});
  indexMat(3,1) = str2double(indices{1}{3});
end

if ~isempty(regexp(trialCat,'^stopTrial_', 'once'))
  if ~isempty(regexp(tagSTOP,'{STOP}', 'once'))
  elseif ~isempty(regexp(tagSTOP,'{STOP:s+\d}', 'once'))
    indices = regexp(tagSTOP,'{STOP:s(\d+)}','tokens');
    indexMat(1,2) = str2double(indices{1}{1});
  elseif ~isempty(regexp(tagSTOP,'{STOP:r+\d}', 'once'))
    indices = regexp(tagSTOP,'{STOP:r(\d+)}','tokens');
    indexMat(2,2) = str2double(indices{1}{1});
  elseif ~isempty(regexp(tagSTOP,'{STOP:c+\d}', 'once'))
    indices = regexp(tagSTOP,'{STOP:c(\d+)}','tokens');
    indexMat(3,2) = str2double(indices{1}{1});
  elseif ~isempty(regexp(tagSTOP,'{STOP:s+\d,r+\d}', 'once'))
    indices = regexp(tagSTOP,'{STOP:s(\d+),r(\d+)','tokens');
    indexMat(1,2) = str2double(indices{1}{1});
    indexMat(2,2) = str2double(indices{1}{2});
  elseif ~isempty(regexp(tagSTOP,'{STOP:s+\d,c+\d}', 'once'))
    indices = regexp(tagSTOP,'{STOP:s(\d+),c(\d+)}','tokens');
    indexMat(1,2) = str2double(indices{1}{1});
    indexMat(3,2) = str2double(indices{1}{2});
  elseif ~isempty(regexp(tagSTOP,'{STOP:r+\d,c+\d}', 'once'))
    indices = regexp(tagSTOP,'{STOP:r(\d+),c(\d+)}','tokens');
    indexMat(2,2) = str2double(indices{1}{1});
    indexMat(3,2) = str2double(indices{1}{2});
  elseif ~isempty(regexp(tagSTOP,'{STOP:s+\d,r+\d,c+\d}', 'once'))
    indices = regexp(tagSTOP,'{STOP:s(\d+),r(\d+),c(\d+)}','tokens');
    indexMat(1,2) = str2double(indices{1}{1});
    indexMat(2,2) = str2double(indices{1}{2});
    indexMat(3,2) = str2double(indices{1}{3});
  end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. CONVERT X TO INDIVIDUAL PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z0    = get_value_per_xcat(SAM,X,indexMat,iZ0,iCatClass,taskFactors);
zc    = get_value_per_xcat(SAM,X,indexMat,iZc,iCatClass,taskFactors);
v     = get_value_per_xcat(SAM,X,indexMat,iV,iCatClass,taskFactors);
ve    = get_value_per_xcat(SAM,X,indexMat,iVe,iCatClass,taskFactors);
eta   = get_value_per_xcat(SAM,X,indexMat,iEta,iCatClass,taskFactors);
t0    = get_value_per_xcat(SAM,X,indexMat,iT0,iCatClass,taskFactors);
se    = get_value_per_xcat(SAM,X,indexMat,iSe,iCatClass,taskFactors);
si    = get_value_per_xcat(SAM,X,indexMat,iSi,iCatClass,taskFactors);
k     = get_value_per_xcat(SAM,X,indexMat,iK,iCatClass,taskFactors);
wliw  = get_value_per_xcat(SAM,X,indexMat,iWliw,iCatClass,taskFactors);
wlib  = get_value_per_xcat(SAM,X,indexMat,iWlib,iCatClass,taskFactors);
wffiw = get_value_per_xcat(SAM,X,indexMat,iWffiw,iCatClass,taskFactors);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3. ENCODE CONNECTIVITY MATRICES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Consider adding an additional parameter class to distinguish lateral inhibition within and between accumulator classes

% Endogenous connectivity
% =========================================================================
if any(features(2,iK,:))
  endoConnSelf            = endoConn.self * (diag(cell2mat(trueNRsp(:))) * diag(cell2mat(k(:))));
else
  endoConnSelf            = endoConn.self * diag(blkdiag(trueNRsp{:}) * cell2mat(k(:)));
end

if any(features(2,iWliw,:))
  endoConnNonSelfSame     = endoConn.nonSelfSame * (diag(cell2mat(trueNRsp(:))) * diag(cell2mat(wliw(:))));
else
  endoConnNonSelfSame     = endoConn.nonSelfSame * diag(blkdiag(trueNRsp{:}) * cell2mat(wliw(:)));
end

if any(features(2,iWlib,:))
  endoConnNonSelfOther    = endoConn.nonSelfOther * (diag(cell2mat(trueNRsp(:))) * diag(cell2mat(wlib(:))));
else
  endoConnNonSelfOther    = endoConn.nonSelfOther * diag(blkdiag(trueNRsp{:}) * cell2mat(wlib(:)));
end

endoConn                  = endoConnSelf + endoConnNonSelfSame + endoConnNonSelfOther;

% Exogenous connectivity
% =========================================================================
exoConnStimTarget         = exoConn.stimTarget * diag(blkdiag(trueNRsp{:}) * ones(numel(maxNRsp),1));

if any(features(2,iWffiw,:))
  exoConnStimNonTargetSame  = diag(diag(cell2mat(wffiw(:))) * cell2mat(trueNRsp(:))) * exoConn.stimNonTargetSame;
else
  exoConnStimNonTargetSame  = diag(blkdiag(trueNRsp{:}) * cell2mat(wffiw(:))) * exoConn.stimNonTargetSame;
end
exoConn                   = exoConnStimTarget + exoConnStimNonTargetSame;

% Extrinsic and intrinsic modulation: none
% =========================================================================
extrMod = zeros(sum(maxNRsp),sum(maxNRsp),sum(maxNStm));
intrMod = zeros(sum(maxNRsp),sum(maxNRsp),sum(maxNRsp));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 4. ENCODE STARTING POINT MATRIX
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if any(features(2,iZ0,:))
  Z0 = diag(cell2mat(z0)) * cell2mat(trueNRsp(:));
else
  Z0 = blkdiag(trueNRsp{:}) * cell2mat(z0(:));
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 5. ENCODE THRESHOLD MATRIX
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if any(features(2,iZc,:))
  ZC = diag(cell2mat(zc(:))) * cell2mat(trueNRsp(:));
else
  ZC = blkdiag(trueNRsp{:}) * cell2mat(zc(:));
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 6. ACCUMULATION RATES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% if any(features(2,iV,:)) | any(features(2,iVe,:))
%   V   = diag(cell2mat(v(:))) * cell2mat(iTarget(:)) + ...
%         diag(cell2mat(ve(:))) * cell2mat(iNonTarget(:));
% else
  V   = blkdiag(iTarget{:}) * cell2mat(v(:)) + ...
        blkdiag(iNonTarget{:}) * cell2mat(ve);
% end

if any(features(2,iEta,:))
  ETA = diag(cell2mat(eta(:))) * cell2mat(iTarget(:)) + ...
        diag(cell2mat(eta(:))) * cell2mat(iNonTarget(:));
else
  ETA = blkdiag(iTarget{:}) * cell2mat(eta(:)) + ...
        blkdiag(iNonTarget{:}) * cell2mat(eta(:));
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 7. EXTRINSIC AND INTRINSIC NOISE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if any(features(2,iSe,:))
  SE = diag(cell2mat(trueNStm(:))) * diag(cell2mat(se(:)));
else
  SE = diag(blkdiag(trueNStm{:}) * cell2mat(se(:)));
end

if any(features(2,iSi,:))
  SI = diag(cell2mat(trueNRsp(:))) * diag(cell2mat(si(:)));
else
  SI = diag(blkdiag(trueNRsp{:}) * cell2mat(si(:)));
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 8. SPECIFY ONSETS AND DURATIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if any(features(2,iT0,:))
  T0 = diag(cell2mat(t0(:))) * cell2mat(trueNStm(:));
else
  T0 = blkdiag(trueNStm{:}) * cell2mat(t0(:));
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #. NESTED FUNCTIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function xval = get_value_per_xcat(SAM,X,iMat,iCol,iCatClass,taskFactors)
    
    % Default value
    xval = num2cell(SAM.model.XCat.valExcluded(iCol)*ones(2,1));

    % GO parameters
    % =====================================================================

    % Vector of parameter category values in X
    valuesGO      = X(iCatClass{1,iCol});
    
    % Signature
    signatureGO   = SAM.model.variants.toFit.features(:,iCol,1);
    
    % Factorial levels
    if all(signatureGO == [0 0 0]')
      levels        = 1;
    else
      levels        = fullfact(taskFactors(signatureGO,1));
    end
    
    % Number of levels
    nLevel        = size(levels,1);
    
    % Matrix of current indices
    thisIMat      = iMat(signatureGO,1)';
    
    % Identify the values to extract
    if ~isempty(valuesGO)
      if all(signatureGO == [0 0 0]')
        iVal      = 1:numel(valuesGO);
      elseif all(signatureGO == [1 0 0]')
        iVal      = all(levels == repmat(thisIMat,nLevel,1),2);
      elseif all(signatureGO == [0 1 0]')
        iVal      = [thisIMat, setxor(thisIMat, 1:numel(valuesGO))]';
%         iVal      = 1:numel(valuesGO);
      elseif all(signatureGO == [0 0 1]')
        iVal      = all(levels == repmat(thisIMat,nLevel,1),2);
      elseif all(signatureGO == [1 1 0]')
        iVal      = all(levels(:,1) == repmat(thisIMat(:,1),nLevel,1),2);
      elseif all(signatureGO == [1 0 1]')
        iVal      = all(levels == repmat(thisIMat,nLevel,1),2);
      elseif all(signatureGO == [0 1 1]')
        iVal      = all(levels == repmat(thisIMat,nLevel,1),2);
%         iVal      = all(levels(:,2) == repmat(thisIMat(:,2),nLevel,1),2);
      elseif all(signatureGO == [1 1 1]')
        iVal      = all(levels(:,[1,3]) == repmat(thisIMat(:,[1,3]),nLevel,1),2);
      end
      xval{1}   = valuesGO(iVal);
      xval{1}   = xval{1}(:);
    end

    % STOP parameters
    % =====================================================================
    switch lower(SAM.sim.scope)
      case {'stop','all'}
        
        % Vector of parameter category values in X
        valuesSTOP    = X(iCatClass{2,iCol});

        % Signature
        signatureSTOP = SAM.model.variants.toFit.features(:,iCol,2);

        % Factorial levels
        if all(signatureSTOP == [0 0 0]')
          levels        = 1;
        else
          levels        = fullfact(taskFactors(signatureSTOP,2));
        end
        
        % Number of levels
        nLevel        = size(levels,1);

        % Matrix of current indices
        thisIMat      = iMat(signatureSTOP,2)';
        
        % Identify the values to extract
        if ~isempty(valuesSTOP)
          if all(signatureSTOP == [0 0 0]')
            iVal      = 1:numel(valuesSTOP);
          elseif all(signatureSTOP == [1 0 0]')
            iVal      = all(levels == repmat(thisIMat,nLevel,1),2);
          elseif all(signatureSTOP == [0 1 0]')
        iVal      = [thisIMat, setxor(thisIMat, 1:numel(valuesGO))]';
%         iVal      = 1:numel(valuesGO);
          elseif all(signatureSTOP == [0 0 1]')
            iVal      = all(levels == repmat(thisIMat,nLevel,1),2);
          elseif all(signatureSTOP == [1 1 0]')
            iVal      = all(levels(:,1) == repmat(thisIMat(:,1),nLevel,1),2);
          elseif all(signatureSTOP == [1 0 1]')
            iVal      = all(levels == repmat(thisIMat,nLevel,1),2);
          elseif all(signatureSTOP == [0 1 1]')
        iVal      = all(levels == repmat(thisIMat,nLevel,1),2);
%         iVal      = all(levels(:,2) == repmat(thisIMat(:,2),nLevel,1),2);
          elseif all(signatureSTOP == [1 1 1]')
            iVal      = all(levels(:,[1,3]) == repmat(thisIMat(:,[1,3]),nLevel,1),2);
          end
          xval{2}   = valuesSTOP(iVal);
          xval{2}   = xval{2}(:);
        end
        
    end    
  
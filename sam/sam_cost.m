function [cost,altCost,prd] = sam_cost(X,SAM)
% SAM_COST <Synopsis of what this function does> 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% SAM_COST; 
%  
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Sat 21 Sep 2013 19:43:58 CDT by bram 
% $Modified: Sat 21 Sep 2013 19:47:47 CDT by bram

 
% CONTENTS 
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS AND SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimScope  = SAM.sim.scope;
nSim      = SAM.sim.n;
obs       = SAM.optim.obs;
costStat  = SAM.optim.cost.stat.stat;
nTrialCat = size(obs,1);
bic       = cell(nTrialCat,1);
chiSquare = cell(nTrialCat,1);

switch lower(optimScope)
  case 'go'
    nFree     = sum(SAM.model.variants.toFit.XSpec.free.go.free);
  case 'stop'
    nFree     = sum(SAM.model.variants.toFit.XSpec.free.stop.free);
  case 'all'
    nFree     = sum(SAM.model.variants.toFit.XSpec.free.all.free);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. SIMULATE EXPERIMENT
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

prd = sam_sim_expt('optimize',X,SAM);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3. COMPUTE COST
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for iTrialCat = 1:nTrialCat
    
    trialCat    = SAM.optim.obs.trialCat{iTrialCat};

    % Go trial
    % =====================================================================
    if ~isempty(regexp(trialCat,'goTrial.*', 'once'))

        bic{iTrialCat} = zeros(2,1);
        chiSquare{iTrialCat} = zeros(2,1);

        % Correct choice
        % -----------------------------------------------------------------
        if obs.nGoCCorr(iTrialCat) > 0
            [bic{iTrialCat}(1), ...
             chiSquare{iTrialCat}(1)] = compute_cost(prd.rtGoCCorr{iTrialCat}, ...
                                                     obs.rtQGoCCorr{iTrialCat}, ...
                                                     obs.probMassDefectiveGoCCorr{iTrialCat}, ...
                                                     nSim, ...
                                                     obs.nGoCCorr(iTrialCat), ...
                                                     nFree);
        end

        % Choice error
        % -----------------------------------------------------------------
        if obs.nGoCError(iTrialCat) > 0
            [bic{iTrialCat}(2), ...
             chiSquare{iTrialCat}(2)] = compute_cost(prd.rtGoCError{iTrialCat}, ...
                                                     obs.rtQGoCError{iTrialCat}, ...
                                                     obs.probMassDefectiveGoCError{iTrialCat}, ...
                                                     nSim, ...
                                                     obs.nGoCError(iTrialCat), ...
                                                     nFree);
        end

    % Stop trial
    % =====================================================================
    elseif ~isempty(regexp(trialCat,'stopTrial.*', 'once'))

        bic{iTrialCat} = zeros(3,1);
        chiSquare{iTrialCat} = zeros(3,1);
        
        % Successful inhibition
        % -----------------------------------------------------------------
        if obs.nStopICorr(iTrialCat) > 0
            [bic{iTrialCat}(1), ...
             chiSquare{iTrialCat}(1)] = compute_cost(prd.rtStopICorr{iTrialCat}, ...
                                                     obs.rtQStopICorr{iTrialCat}, ...
                                                     obs.probMassDefectiveStopICorr{iTrialCat}, ...
                                                     nSim, ...
                                                     obs.nStopICorr(iTrialCat), ...
                                                     nFree);
        end
        
        % Failed inhibition, correct choice
        % -----------------------------------------------------------------
        if obs.nStopIErrorCCorr(iTrialCat) > 0
            [bic{iTrialCat}(2), ...
             chiSquare{iTrialCat}(2)] = compute_cost(prd.rtStopIErrorCCorr{iTrialCat}, ...
                                                     obs.rtQStopIErrorCCorr{iTrialCat}, ...
                                                     obs.probMassDefectiveStopIErrorCCorr{iTrialCat}, ...
                                                     nSim, ...
                                                     obs.nStopIErrorCCorr(iTrialCat), ...
                                                     nFree);
        end
        
        % Failed inhibition, choice error
        % -----------------------------------------------------------------
        if obs.nStopIErrorCError(iTrialCat) > 0
            [bic{iTrialCat}(3), ...
             chiSquare{iTrialCat}(3)] = compute_cost(prd.rtStopIErrorCError{iTrialCat}, ...
                                                     obs.rtQStopIErrorCError{iTrialCat}, ...
                                                     obs.probMassDefectiveStopIErrorCError{iTrialCat}, ...
                                                     nSim, ...
                                                     obs.nStopIErrorCError(iTrialCat), ...
                                                     nFree);
        end

    end

end

% 3.1. Log all BIC and chi-square values
% =========================================================================
prd.bic         = bic;
prd.chiSquare   = chiSquare;

% 3.3. Compute total
% =========================================================================

allBic          = cell2mat(bic);
allChiSquare    = cell2mat(chiSquare);

switch lower(costStat)
  case 'bic'
    cost    = sum(allBic);
    altCost = sum(allChiSquare);
  case 'chisquare'
    cost    = sum(allChiSquare);
    altCost = sum(allBic);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 4. NESTED FUNCTIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [bic,chiSquare] = compute_cost(rtPrd,rtQObs,probMassDefectiveObs,nSim,nTrialObs,nFree)
  
    % Compute the predicted probability masses
    if isempty(rtPrd)
        probMassDefectivePrd = zeros(size(probMassDefectiveObs));
    else
        probMassPrd = histc(rtPrd,[-Inf,rtQObs,Inf]);
        probMassPrd = probMassPrd(1:end-1);
        probMassDefectivePrd = probMassPrd./nSim;
    end
    
    % Add a small value to bins with a probablity mass of 0 (to prevent
    % division by 0 and hampering optimization)
    probMassDefectivePrd(probMassDefectivePrd == 0) = 1e-4;
    
    % Compute the costs
    chiSquare     = sam_chi_square(probMassDefectiveObs(:),probMassDefectivePrd(:),nTrialObs);
    bic           = sam_bic(probMassDefectiveObs(:),probMassDefectivePrd(:),nTrialObs,nFree);
%     bic           = sam_bic(probMassDefectiveObs(:),probMassDefectivePrd(:),length(rtQObs),nFree);
  
end
end
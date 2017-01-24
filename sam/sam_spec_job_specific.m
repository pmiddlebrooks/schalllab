function SAM = sam_spec_job_specific(SAM,iModel);
% SAM_SPEC_JOB_SPECIFIC <Synopsis of what this function does> 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% SAM_SPEC_JOB_SPECIFIC; 
%  
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Mon 17 Mar 2014 15:05:32 CDT by bram 
% $Modified: Mon 17 Mar 2014 15:05:32 CDT by bram 

% Go to work directory
cd(SAM.io.workDir);

optimScope                    = SAM.sim.scope;

nRsp                          = SAM.expt.nRsp;
nCnd                          = SAM.expt.nCnd;
maxNRsp                       = max(cell2mat(nRsp(:)),[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the model to be fit
SAM.model.variants.toFit      = SAM.model.variants.tree(iModel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Categorize the observations into categories corresponding with the model specifics
% =========================================================================================================================
SAM.optim.obs                 = sam_categorize_data(SAM);

% Specify initial parameters, based on which parameter constraints are defined
% =========================================================================================================================
SAM.optim.x0Base              = sam_get_x0(SAM);

% Specify parameter constraints
% =========================================================================================================================
SAM.optim.constraint          = sam_get_constraint(SAM);

switch lower(optimScope)
  case {'go','stop'}
    
    % Sample 20 uniformly distributed starting points, given constraints
    % =========================================================================================================================
    nStartPoint                   = SAM.optim.nStartPoint;

    LB                            = SAM.optim.constraint.bound.LB;
    UB                            = SAM.optim.constraint.bound.UB;
    A                             = SAM.optim.constraint.linear.A;
    b                             = SAM.optim.constraint.linear.b;
    nonLinCon                     = SAM.optim.constraint.nonlinear.nonLinCon;

    SAM.optim.x0                  = [SAM.optim.x0Base; ...
                                    sam_sample_uniform_constrained_x0(nStartPoint,LB,UB,A,b,nonLinCon,SAM.optim.solver.type)];

  case 'all'
        
    SAM.optim.x0                  = SAM.optim.x0Base;
    
end

% Model predictions
% =========================================================================================================================
nTrialCat                     = size(SAM.optim.obs,1);
SAM.optim.prd                 = dataset({cell(nTrialCat,1),'trialCat'}, ...
                                        {cell(nTrialCat,1),'funGO'}, ...
                                        {cell(nTrialCat,1),'funSTOP'}, ...
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
                                        {cell(nTrialCat,1),'bic'}, ...
                                        {cell(nTrialCat,1),'chiSquare'}, ...
                                        {cell(nTrialCat,1),'modelMat'}, ...
                                        {cell(nTrialCat,1),'dyn'});
SAM.optim.prd.trialCat        = SAM.optim.obs.trialCat;
SAM.optim.prd.onset           = SAM.optim.obs.onset;
SAM.optim.prd.duration        = SAM.optim.obs.duration;
SAM.optim.prd.ssd             = SAM.optim.obs.ssd;

% Maximum number of function evaluations and iterations
% =========================================================================================================================
switch lower(optimScope)
  case 'go'
    nFree = sum(SAM.model.variants.toFit.XSpec.free.go.free);
  case 'stop'
    nFree = sum(SAM.model.variants.toFit.XSpec.free.stop.free);
  case 'all'
    nFree = sum(SAM.model.variants.toFit.XSpec.free.all.free);
end
SAM.optim.solver.opts.MaxFunEvals = 1000 * nFree;
SAM.optim.solver.opts.MaxIter     = 1000 * nFree;
SAM.optim.solver.opts.TolFun      = 1e-5;
SAM.optim.solver.opts.TolX        = 1e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. HACKS (TO BE IMPLEMENTED ELSEWHERE SOON)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dataset array
modelMat                            = dataset({cell(nTrialCat,1),'trialCat'}, ...
                                              {cell(nTrialCat,1),'iTarget'}, ...
                                              {cell(nTrialCat,1),'iNonTarget'});
modelMat.trialCat                   = SAM.optim.obs.trialCat;

for iCnd = 1:nCnd
  
  % Identify trials
  iGoTrial                          = cell2mat(cellfun(@(inp1) ~isempty(regexp(inp1,sprintf('^goTrial.*c%d',iCnd))),modelMat.trialCat,'Uni',0));
  iStopTrial                        = cell2mat(cellfun(@(inp1) ~isempty(regexp(inp1,sprintf('^stopTrial.*c%d',iCnd))),modelMat.trialCat,'Uni',0));
  
  % Identify target accumulators
  iTargetGoTrial                    = cellfun(@(inp1) logical(inp1(:)),{[1,zeros(1,maxNRsp(1)-1)] zeros(1,maxNRsp(2))},'Uni',0);
  iTargetStopTrial                  = cellfun(@(inp1) logical(inp1(:)),{[1,zeros(1,maxNRsp(1)-1)] [1,zeros(1,maxNRsp(2)-1)]},'Uni',0);
  
  modelMat.iTarget(iGoTrial)        = cellfun(@(inp1) iTargetGoTrial,modelMat.iTarget(iGoTrial),'Uni',0);
  modelMat.iTarget(iStopTrial)      = cellfun(@(inp1) iTargetStopTrial,modelMat.iTarget(iStopTrial),'Uni',0);
  
  % Identify non-target accumulators
  nNonTargetGo                      = nRsp{iCnd}(1)-1;
  nNonTargetStop                    = nRsp{iCnd}(2)-1;
  nEmptyGo                          = maxNRsp(1) - nNonTargetGo - 1;
  nEmptyStop                        = maxNRsp(2) - nNonTargetStop - 1;
  
  iNonTargetGoTrial                 = cellfun(@(inp1) logical(inp1(:)),{[0,ones(1,nNonTargetGo),zeros(1,nEmptyGo)] zeros(1,maxNRsp(2))},'Uni',0);
  iNonTargetStopTrial               = cellfun(@(inp1) logical(inp1(:)),{[0,ones(1,nNonTargetGo),zeros(1,nEmptyGo)] [0,ones(1,nNonTargetStop),zeros(1,nEmptyStop)]},'Uni',0);
  
  modelMat.iNonTarget(iGoTrial)     = cellfun(@(inp1) iNonTargetGoTrial,modelMat.iNonTarget(iGoTrial),'Uni',0);
  modelMat.iNonTarget(iStopTrial)   = cellfun(@(inp1) iNonTargetStopTrial,modelMat.iNonTarget(iStopTrial),'Uni',0);
  
end

SAM.optim.modelMat                  = modelMat;
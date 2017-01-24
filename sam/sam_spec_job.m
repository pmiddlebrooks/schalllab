% Script specifying job details
%  
% DESCRIPTION 
% This script contains all the details for the job to run
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Mon 09 Sep 2013 13:07:49 CDT by bram 
% $Modified: Sat 21 Sep 2013 12:24:04 CDT by bram

% Starting parameter index
% -------------------------------------------------------------------------
% There may be a set of starting parameters at which the optimization
% algorithm starts

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS, DEFINE VARIABLES, ADD ACCESS TO PATHS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1.1. Process inputs
% =========================================================================

% Get path string
pathStr   = getenv('pathStr');

% Get optimScope
optimScope  = getenv('optimScope');

% Get model index
iModel    = str2double(getenv('iModel'));

% Get subject index
iSubj     = str2double(getenv('subject'));

% Get starting point index
iStartVal = str2double(getenv('iStartVal'));

% 1.2. Define dynamic variables
% =========================================================================
timeStr   = datestr(now,'yyyy-mm-dd-THHMMSS');

% 1.3. Add paths
% =========================================================================
addpath(genpath('/home/zandbeb/m-files/sam/'));
addpath(genpath('/home/zandbeb/m-files/matlab_code_bbz/'));
addpath(genpath('/home/zandbeb/m-files/matlab_file_exchange_tools/'));
addpath(genpath('/home/zandbeb/m-files/cmtb/'));
addpath(genpath('/scratch/zandbeb/multichoice_stop_model/src/code'));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. RUN AND SAVE JOB
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 2.1. Load the SAM file
% =========================================================================
load(sprintf(pathStr,iSubj,optimScope,iModel));

% 2.2. Add details for logging
% =========================================================================
fNameIterLog                  = sprintf('iterLog_%sTrials_model%.3d_startVal%.3d_started%s.mat',optimScope,iModel,iStartVal,timeStr);
fNameFinalLog                 = sprintf('finalLog_%sTrials_model%.3d_startVal%.3d_started%s.mat',optimScope,iModel,iStartVal,timeStr);

% Iteration log file
fitLog.iterLogFile            = fullfile(SAM.io.workDir,fNameIterLog);

% Iteration lof frequency
fitLog.iterLogFreq            = 50;

% Final log file
fitLog.finalLogFile           = fullfile(SAM.io.workDir,fNameFinalLog);

SAM.optim.log                 = fitLog;

% 2.3. Optimize the fit to the data, starting from parameters corresponding to iStartVal
% =========================================================================
SAM                           = sam_optim(SAM,iStartVal);

% 2.4. Optimize the fit to the data, starting from parameters corresponding to iStartVal
% =========================================================================
fNameSAM                      = sprintf('SAM_%sTrials_model%.3d_exit%d_started%s.mat',optimScope,iModel,SAM.estim.exitFlag,timeStr);
fNameX                        = sprintf('bestX_%sTrials_model%.3d_exit%d_started%s.txt',optimScope,iModel,SAM.estim.exitFlag,timeStr);

save(fullfile(SAM.io.workDir,fNameSAM),'SAM');
save(fullfile(SAM.io.workDir,fNameX),'X','-ascii','-double');
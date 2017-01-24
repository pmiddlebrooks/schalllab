function sam_make_obs_file
% SAM_MAKE_OBS_FILE generates an observation file for each model
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% It is assumed that 
%
%
%
% SYNTAX 
% SAM_MAKE_OBS_FILE; 
%  
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Thu 23 Jan 2014 15:45:35 CST by bram 
% $Modified: Thu 23 Jan 2014 15:45:35 CST by bram 

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Static variables (to be replaced)
% ------------------------------------------------------------------------- 

nStm = [6 1];
nRsp = [6 1];
nCnd = 3;

% How is SSD manipulated?
% - across conditions
% ssdManip = [0 0 1]'
% nSsd = 5;
%
% - across go responses (e.g. a staircase for left and right hands)
% ssdManip = [0 1 0;0 0 0]'
%
% - independent from task factors
%


% 
% ========================================================================= 

 
% 
% ------------------------------------------------------------------------- 


% Loop over models
for iModel = 1:nModel
    
  % Get model features
  features = models(iModel).features;
  
  % Compute number of cells for go trials (correct and commission)
  nGoTrialCells = diag([nStm(1);nRsp(1);nCnd]) * any(features(:,:,1),2);
  nGoTrialCells (nGoTrialCells == 0) = 1;
  nGoTrialCells = prod(nGoTrialCells);
  
  % Compute number of cells for stop trials (correct/inhibit and commission/respond)
  % MAKE THIS DYNAMIC!
  nStopTrialCells = 15;
  
  nCells = nGoTrialCells + nStopTrialCells;
  
end


% Put the data into those cells
% =========================================================================




s1CorrData = cell(nStm(1),nRsp(1),nCnd);
s1CommData = cell(nStm(1),nRsp(1),nCnd);
for iStm = 1:nStm(1)
  for iRsp = 1:nRsp(1)
    for iCnd = 1:nCnd
      
      s1CorrData{iStm,iRsp,iCnd} = find(data.stm1 == iStm & ...
                                    data.stm1 == iRsp & ...
                                    data.cnd1 == iCnd & ...
                                    data.stm2 == 0 & ...
                                    data.acc == 2);
                                  
      s1CommData{iStm,iRsp,iCnd} = find(data.stm1 == iStm & ...
                                    data.stm1 == iRsp & ...
                                    data.cnd1 == iCnd & ...
                                    data.stm2 == 0 & ...
                                    data.acc == 0);
    end
  end
end 
      

fullfact(E(any(features(:,:,1),2)))

goCell = cell(nGoTrialCells,1);



s2CorrData = cell(nStm(1),nRsp(1),nCnd);
s2CommData = cell(nStm(1),nRsp(1),nCnd);
for iStm = 1:nStm(1)
  for iRsp = 1:nRsp(1)
    for iCnd = 1:nCnd
      
      s2CorrData{iStm,iRsp,iCnd} = find(data.stm1 == iStm & ...
                                        data.stm1 == iRsp & ...
                                        data.cnd1 == iCnd & ...
                                        data.stm2 == 0 & ...
                                        data.acc == 2);
                                  
      s2CommData{iStm,iRsp,iCnd} = find(data.stm1 == iStm & ...
                                        data.stm1 == iRsp & ...
                                        data.cnd1 == iCnd & ...
                                        data.stm2 == 0 & ...
                                        data.acc == 0);
    end
  end
end 




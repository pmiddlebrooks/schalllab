function obs = ccm_sam_translate_datafile(subjectID, sessionID, sessionSet)
%
% TRANSLATE_DATAFILE_SAM:   Translates an existing choice countermanding
% datafile into a structure suitable to input into the SAM (Stochastic
% accumulator model)
% toolbox
%
% DESCRIPTION
% Creates a dataset called "obs", with as many rows as there are task
% conditions (e.g. 6 levels of signal strength)
%
% SYNTAX
% SAM_RUN_JOB;
%
% EXAMPLES
%
%
% REFERENCES
%
% .........................................................................
% paul middlebrooks     paul.g.middlebrooks@vanderbilt.edu
% $Created : Sat 21 Sep 2013 12:48:45 CDT by bram
% $Modified: Sat 21 Sep 2013 12:56:52 CDT by bram


% CONTENTS
% 1.PROCESS INPUTS AND SPECIFY VARIABLES
%   1.1. Process inputs
%   1.2. Pre-allocate empty arrays
% 2.FILL DATASET ARRAY
% 3.CREATE DATASET




if nargin < 3
    sessionSet = 'behavior1';
end
subSampleFlag = 1;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS AND SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.1. Process inputs
% =========================================================================

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
pSignalArray = ExtraVar.pSignalArray;
ssdArray = ExtraVar.ssdArray;
% pSignalArray = pSignalArray(pSignalArray ~= .5);
switch lower(subjectID)
    case 'human'
        pSignalArray = [.35 .42 .46 .54 .58 .65];
    case 'broca'
        switch sessionSet
            case 'behavior1'
                pSignalArray = [.41 .45 .48 .52 .55 .59];
            case 'neural1'
                pSignalArray = [.41 .44 .47 .53 .56 .59];
            case 'neural2'
                pSignalArray = [.42 .44 .46 .54 .56 .58];
            otherwise
                pSignalArray = ExtraVar.pSignalArray;
        end
    case 'xena'
        switch sessionSet
            case 'behavior'
                pSignalArray = [.35 .42 .47 .53 .58 .65];
                trialData.targ1CheckerProp(trialData.targ1CheckerProp == .52) = .53;
        end
end


% Get local data path for saving obs file below
[tebDataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessionID);

% Truncate RTs
rtMin                       = 140;
rtMax                       = 1200;
[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTD);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);


% Number of conditions
nCnd        = length(pSignalArray);

% Number of stop-signal delays
nSsd        = length(ssdArray);

% Number of go and stop units
nGoUnit     = 2;
nStopUnit   = 1;



% 1.2. Define and preallocate arrays for the dataset
% =========================================================================
nGo             = nan(nCnd, 1);     % # of total Go trials
nGoCorr         = nan(nCnd, 1);     % # of Go Target trials
nGoComm         = nan(nCnd, 1);     % # of Go Distractor trials
nStop           = nan(nCnd, nSsd); 	% # of stop trials
nStopFailureCorr   	= nan(nCnd, nSsd);	% # of noncanceled stop trials
nStopFailureComm   	= nan(nCnd, nSsd);	% # of noncanceled stop trials
nStopSuccess   	= nan(nCnd, nSsd); 	% # of canceled stop trials
pGoCorr         = nan(nCnd, 1);     % proportion of Go Target Trials
pGoComm         = nan(nCnd, 1);     % proportion of Go Distractor Trials
pStopFailure   	= nan(nCnd, nSsd); 	% inhibition function
pStopFailureCorr   	= nan(nCnd, nSsd); 	% inhibition function
pStopFailureComm   	= nan(nCnd, nSsd); 	% inhibition function
ssd             = nan(nCnd, nSsd); 	% SSDs
inhibFunc   	= nan(nCnd, nSsd);  % inhibition function
rtGoCorr        = cell(nCnd, 1);  	% Go Target RTs
rtGoComm        = cell(nCnd, 1);    % Go Distractor RTs
rtStopFailureCorr   = cell(nCnd, nSsd); % Noncanceled RTs
rtStopFailureComm   = cell(nCnd, nSsd); % Noncanceled RTs
onset           = cell(nCnd, 1 + nSsd); % Response Cue (checker) Onset times
duration      	= cell(nCnd, 1 + nSsd); % Response Cue Duration times






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. FILL DATASET ARRAY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iCnd = 1 : nCnd
    
    % Signal strength/coherence (right target checker color pct)
    iPct = pSignalArray(iCnd) * 100;
    
    % Get default trial selection options
    selectOpt       = ccm_trial_selection;
    selectOpt.rightCheckerPct = iPct;
    
    
    
    % Number of Go Target trials
    selectOpt.ssd       = 'none';
    selectOpt.targDir    	= 'collapse';
    selectOpt.outcome     = {'goCorrectTarget', 'targetHoldAbort'};
    goTargTrial         = ccm_trial_selection(trialData, selectOpt);
    nGoCorr(iCnd)       = length(goTargTrial);
    
    % Number of Go Distractor trials (Errors of Commission)
    selectOpt.outcome     = {'goCorrectDistractor', 'distractorHoldAbort'};
    goDistTrial         = ccm_trial_selection(trialData, selectOpt);
    nGoComm(iCnd)       = length(goDistTrial);
    
    % Number of total Go trials
    nGo(iCnd)           = nGoCorr(iCnd) + nGoComm(iCnd);
    
    % proportion of Go Target Trials
    pGoCorr(iCnd)    	= nGoCorr(iCnd) / nGo(iCnd);
    
    % proportion of Go Distractor Trials
    pGoComm(iCnd)     	= nGoComm(iCnd) / nGo(iCnd);
    
    % Correct and Error RTs
    rtGoCorr{iCnd}  	= trialData.rt(goTargTrial)';
    rtGoComm{iCnd}  	= trialData.rt(goDistTrial)';
    
    % Onset and duration of stimuli (relative to stimulus onset for now)
    onset{iCnd, 1}      = ones(nGoUnit + nStopUnit, 1);
    duration{iCnd, 1}  	= [2250 * ones(nGoUnit, 1); 0];
    
    
    % ************************
    % ************************
    % ************************
    % ************************
    for jSsdInd = 1 : nSsd
        
        % Stop signal delay
        jSsd = ssdArray(jSsdInd);
        selectOpt.ssd       = jSsd;
        
        % ************************************************
        % Differentiate stop to target vs stop to distractor?
        % ************************************************
        %         % Number of Stop Target trials (errors of ommission)
        %         stopTargOutcome         = {'stopIncorrectTarget', 'targetHoldAbort'};
        %         stopTargTrial           = ccm_trial_selection(trialData, stopTargOutcome, pct, ssd, targetHemifield);
        %         nStopCorr(nCnd, jSsd) 	= length(stopTargTrial);
        %
        %         % Number of Stop Distractor trials (Errors of Commission and Ommission)
        %         stopDistOutcome         = {'stopIncorrectDistractor', 'distractorHoldAbort'};
        %         stopDistTrial           = ccm_trial_selection(trialData, stopDistOutcome, pct, ssd, targetHemifield);
        %         nStopComm(nCnd, jSsd) 	= length(stopDistTrial);
        
        % Number of Stop Target trials (errors of ommission)
        selectOpt.outcome       = {'stopIncorrectTarget', 'targetHoldAbort', 'stopIncorrectPreSSDTarget'};
        stopTargTrial           = ccm_trial_selection(trialData, selectOpt);
        selectOpt.outcome       =  {'stopIncorrectDistractor', 'distractorHoldAbort', 'stopIncorrectPreSSDDistractor'};
        stopDistTrial           = ccm_trial_selection(trialData, selectOpt);
        nStopFailureCorr(iCnd, jSsdInd) 	= length(stopTargTrial);
        nStopFailureComm(iCnd, jSsdInd) 	= length(stopDistTrial);
        
        % Number of Stop Distractor trials (Errors of Commission and Ommission)
        selectOpt.outcome       =  {'stopCorrect'};
        stopCorrTrial           = ccm_trial_selection(trialData, selectOpt);
        nStopSuccess(iCnd, jSsdInd) = length(stopCorrTrial);
        
        
        % Number of total Stop trials
        nStop(iCnd, jSsdInd)       = nStopFailureCorr(iCnd, jSsdInd) + nStopFailureComm(iCnd, jSsdInd) + nStopSuccess(iCnd, jSsdInd);
        
        
        % ************************************************
        % Differentiate stop to target vs stop to distractor?
        % ************************************************
        %         % proportion of Stop Target Trials
        %         pStopFailure(iCnd, jSsd) = (nStopCorr(iCnd, jSsd) / nStop(iCnd, jSsd);
        %
        %         % proportion of Stop Distractor Trials
        %         pStopFailure(iCnd, jSsd) = (nStopComm(iCnd, jSsd) / nStop(iCnd, jSsd);
        
        % proportion of Stop Failure Trials
        pStopFailure(iCnd, jSsdInd)    = (nStopFailureCorr(iCnd, jSsdInd) + nStopFailureComm(iCnd, jSsdInd)) / nStop(iCnd, jSsdInd);
        pStopFailureCorr(iCnd, jSsdInd)    = nStopFailureCorr(iCnd, jSsdInd) / nStop(iCnd, jSsdInd);
        pStopFailureComm(iCnd, jSsdInd)    = nStopFailureComm(iCnd, jSsdInd) / nStop(iCnd, jSsdInd);
        inhibFunc(iCnd, jSsdInd)       = pStopFailure(iCnd, jSsdInd);
        
        %         % proportion of Stop Success Trials
        %         pStopSuccess(iCnd, jSsd) 	= nStopSuccess(iCnd, jSsd) / nStop(iCnd, jSsd);
        
        
        % Correct and Error RTs
        rtStopFailureCorr{iCnd, jSsdInd}            = trialData.rt(stopTargTrial)';
        rtStopFailureComm{iCnd, jSsdInd}            = trialData.rt(stopDistTrial)';
        
        
        % Onset and duration of stimuli (relative to stimulus onset for now)
        onset{iCnd, 1+jSsdInd}      = [ones(nGoUnit, 1); jSsd];
        duration{iCnd, 1+jSsdInd}  	= [2250 * ones(nGoUnit, 1); 2250 - jSsd];
        
    end % for jSsd = 1 : nSsd
    
    
    % SSDs in condition iCnd
    ssd(iCnd, :)        = ssdArray;
    
end % for iCnd = 1 : cCnd



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. CREATE DATASET
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obs  = dataset({nGo,'nGo'}, ...
    {nGoCorr,'nGoCorr'}, ...
    {nGoComm,'nGoComm'}, ...
    {nStop,'nStop'}, ...
    {nStopFailureCorr,'nStopFailureCorr'}, ...
    {nStopFailureComm,'nStopFailureComm'}, ...
    {nStopSuccess,'nStopSuccess'}, ...
    {pGoCorr,'pGoCorr'}, ...
    {pGoComm,'pGoComm'}, ...
    {pStopFailure,'pStopFailure'}, ...
    {pStopFailureCorr,'pStopFailureCorr'}, ...
    {pStopFailureComm,'pStopFailureComm'}, ...
    {ssd,'ssd'}, ...
    {inhibFunc,'inhibFunc'}, ...
    {rtGoCorr,'rtGoCorr'}, ...
    {rtGoComm,'rtGoComm'}, ...
    {rtStopFailureCorr,'rtStopFailureCorr'}, ...
    {rtStopFailureComm,'rtStopFailureComm'}, ...
    {onset,'onset'}, ...
    {duration,'duration'});



fName = ['~/schalllab/local_data/sam/', lower(subjectID), '/obs_', sessionID, '_sam.mat'];
save(fName, 'obs', '-mat')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. IF WE WANT TO SUB-SAMPLE THE SSDS (B/C THERE ARE SO MANY OTHERWISE...
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if subSampleFlag
    %     method = 'peaks';
    method = 'fraction';
    fractionThresh = .01;
    
    % Subsample the SSDs for simulation purposes
    %    allStop     = obs.nStopSuccess + obs.nStopFailureCorr + obs.nStopFailureComm;
    allStop     = obs.nStop;
    allStop     = sum(allStop);
    pStop       = allStop ./ sum(allStop);
    
    
    switch method
        case 'fraction'
            subIndices = pStop >= fractionThresh;
        case 'peaks'
            [pks, loc] 	= findpeaks(pStop);
            
            maxSSDs = 7;
            if length(loc) > maxSSDs
                [m, i] = max(allStop(loc(maxSSDs)+2:end))
                i = i + loc(maxSSDs)+1;
                loc = [loc(1:maxSSDs), i];
            end
            loc(loc == 1 | loc == 2) = [];
            
            
            switch lower(subjectID)
                case {'broca','xena'}
                    subIndices = [1 2 loc];  % Use the first 2 ssds and whatever values findpeaks returned
                case {'human'}
                    subIndices = [2:2:20,24];
            end
            %    subPStop    = pStop(subIndices);
            %    normSubPStop = subPStop ./ max(subPStop);
            
    end % switch method
    
    % Go through and subsample the relevant variables in the observation
    % dataset:
    obs.ssd                 = obs.ssd(:, subIndices);
    obs.onset               = obs.onset(:, [1 subIndices + 1]);  % add one to subIndices to account for first trial type being a GO
    obs.duration            = obs.duration(:, [1 subIndices + 1]);  % add one to subIndices to account for first trial type being a GO
    obs.nStop               = obs.nStop(:, subIndices);
    obs.nStopFailureCorr    = obs.nStopFailureCorr(:, subIndices);
    obs.nStopFailureComm    = obs.nStopFailureComm(:, subIndices);
    obs.nStopSuccess        = obs.nStopSuccess(:, subIndices);
    obs.pStopFailure        = obs.pStopFailure(:, subIndices);
    obs.pStopFailureCorr    = obs.pStopFailureCorr(:, subIndices);
    obs.pStopFailureComm    = obs.pStopFailureComm(:, subIndices);
    obs.inhibFunc           = obs.inhibFunc(:, subIndices);
    obs.rtStopFailureCorr   = obs.rtStopFailureCorr(:, subIndices);
    obs.rtStopFailureComm   = obs.rtStopFailureComm(:, subIndices);
    
    fName = ['~/schalllab/local_data/sam/', lower(subjectID), '/obs_', sessionID, '_sam_sub.mat']
    save(fName, 'obs', '-mat')
end





end % function





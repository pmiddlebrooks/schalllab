function ccm_concat_sessions_behavior(subjectID, sessionArray)
%%
% subjectID = 'human';
% subjectID = 'xena';
subjectID = 'broca';

% sessionSet = 'neural2';
sessionSet = 'behavior2';
% sessionSet = 'neural3';

task = 'ccm';
% if nargin < 2
[sessionArray, subjectIDArray] = task_session_array(subjectID, task, sessionSet);
% end
%
% subjectID = 'pm'
% sessionSet = 'allsaccade';
% sessionArray = {'pmAllsaccade'}
% subjectIDArray = {'pm'}
nSession = length(sessionArray);


sessionTag      = [];
responseOnset   = [];
targ1CheckerProp = [];
responseCueOn   = [];
stopSignalOn    = [];
trialOutcome    = [];
targOn          = [];
targAngle       = [];
targAmp         = [];
distAngle       = [];
ssd             = [];
saccAmp        	= [];
saccAngle     	= [];
saccBegin      	= [];
checkerOn      	= [];
saccToTargIndex = [];

for iSession = 1 : nSession
    
    % Load the data
    iSessionID = sessionArray{iSession};
    iSubjectID = subjectIDArray{iSession};
    %     if strcmp('human',subjectID)
    %         iSessionID = ['hu',iSessionID];
    %     end
    variables = [ccm_min_vars, {'targAmp', 'distAngle', 'saccAmp', 'saccBegin'} ];
 [trialData, SessionData, ~] = load_data(iSubjectID, iSessionID, variables);
    
    %     if ~strcmp(SessionData.taskID, 'ccm')
    %         error('Not a choice countermanding session, try again\n')
    %     end
    
    
    iSessionTag         = iSession * ones(size(trialData.trialOutcome, 1), 1);
    sessionTag          = [sessionTag; iSessionTag];
    responseOnset       = [responseOnset; trialData.responseOnset];
    targ1CheckerProp    = [targ1CheckerProp; trialData.targ1CheckerProp];
    responseCueOn       = [responseCueOn; trialData.responseCueOn];
    stopSignalOn        = [stopSignalOn; trialData.stopSignalOn];
    ssd                 = [ssd; trialData.ssd];
    targOn              = [targOn; trialData.targOn];
    trialOutcome        = [trialOutcome; trialData.trialOutcome];
    targAngle           = [targAngle; trialData.targAngle];
    targAmp             = [targAmp; trialData.targAmp];
    distAngle           = [distAngle; trialData.distAngle];
    saccAmp             = [saccAmp; trialData.saccAmp];
    saccAngle           = [saccAngle; trialData.saccAngle];
    saccBegin          	= [saccBegin; trialData.saccBegin];
    checkerOn         	= [checkerOn; trialData.checkerOn];
    saccToTargIndex    	= [saccToTargIndex; trialData.saccToTargIndex];
    
end


trialData = dataset();
trialData.sessionTag        = sessionTag;
trialData.responseOnset     = responseOnset;
trialData.targ1CheckerProp  = targ1CheckerProp;
trialData.responseCueOn     = responseCueOn;
trialData.stopSignalOn      = stopSignalOn;
trialData.ssd               = ssd_session_adjust(ssd);
trialData.trialOutcome      = trialOutcome;
trialData.targAngle         = targAngle;
trialData.targAmp           = targAmp;
trialData.distAngle         = distAngle;
trialData.saccAngle           = saccAngle;
trialData.saccAmp           = saccAmp;
trialData.saccBegin         = saccBegin;
trialData.checkerOn         = checkerOn;
trialData.saccToTargIndex  	= saccToTargIndex;

SessionData.taskID = 'ccm';

% switch lower(subjectID)
%     case 'human'
%         pSignalArray = [.35 .42 .46 .54 .58 .65];
%     case 'broca'
%         switch sessionSet
%             case 'behavior1'
%                 pSignalArray = [.41 .45 .48 .52 .55 .59];
%             case 'neural1'
%                 pSignalArray = [.41 .44 .47 .53 .56 .59];
%             case 'neural2'
%                 pSignalArray = [.42 .44 .46 .54 .56 .58];
%             otherwise
%                 pSignalArray = ExtraVar.pSignalArray;
%         end
%     case 'xena'
%         switch sessionSet
%             case 'behavior'
%                 pSignalArray = [.35 .42 .47 .53 .58 .65];
%                 trialData.targ1CheckerProp(trialData.targ1CheckerProp == .52) = .53;
%         end
% end

% new = trialData;

if strcmp(subjectID, 'human')
    save(['~/schalllab/local_data/human/human_', iSessionID, '.mat'], 'SessionData', 'trialData', '-mat')
else
    save(['~/schalllab/local_data/', subjectID, '/', subjectID, '_', sessionSet, '.mat'], 'trialData', 'SessionData', '-mat')
end



function [SessionData] = cmd_session_behavior(subjectID, sessionID, varargin)
%%

% Set defaults
plotFlag = 1;
printPlot = 0;
for i = 1 : 2 : length(varargin)
    switch varargin{i}
        case 'plotFlag'
            plotFlag = varargin{i+1};
        case 'printPlot'
            printPlot = varargin{i+1};
        otherwise
    end
end

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
ssdArray = ExtraVar.ssdArray;
nTrial = size(trialData.trialOutcome, 1)

if ~strcmp(SessionData.taskID, 'cmd')
    fprintf('Not a simpple countermanding session, try again\n')
    return
end






if plotFlag
    figureHandle = 674;   
end










% ***********************************************************************
% Inhibition Function:
%       &
% SSD vs. Proportion of Response trials
%       &
% SSD vs. Proportion(Correct Choice)
% ***********************************************************************


disp('Stopping data')
dataInh = cmd_inhibition(subjectID, sessionID, 'printPlot', printPlot, 'figureHandle', figureHandle);

% pStopRespond     = dataInh.pStopRespond;
% nStop               = dataInh.nStop;
% nStopStop           = dataInh.nStopStop;
% nStopTarg           = dataInh.nStopTarg;
% nStopDist           = dataInh.nStopDist;
% stopTargetProb      = dataInh.stopTargetProb;
% inhibitionFn        = dataInh.inhibitionFn;
% ssrtGrand                = dataInh.ssrtGrand;
% ssrtMean                = dataInh.ssrtMean;
% ssrtIntegration                = dataInh.ssrtIntegration;
% ssrtIntegrationWeighted                = dataInh.ssrtIntegrationWeighted;
% ssrtIntegrationSimple               = dataInh.ssrtIntegrationSimple;
% pStopRespondGrand = dataInh.pStopRespondGrand;
% inhibitionFnGrand   = dataInh.inhibitionFnGrand;
% ssrtCollapseGrand           = dataInh.ssrtCollapseGrand;
% ssrtCollapseIntegrationWeighted           = dataInh.ssrtCollapseIntegrationWeighted;
% ssrtCollapseIntegration           = dataInh.ssrtCollapseIntegration;
% ssrtCollapseMean           = dataInh.ssrtCollapseMean;
% 
% 













SessionData.ssdArray        = ssdArray;
SessionData.pStopRespond    = dataInh.pStopRespond;
% SessionData.nStop           = nStop;
% SessionData.nStopStop       = nStopStop;
% SessionData.nStopTarg       = nStopTarg;
% SessionData.nStopDist       = nStopDist;
% SessionData.goTargetProb    = goTargetProb;
% SessionData.stopTargetProb  = stopTargetProb;
% SessionData.inhibitionFn    = inhibitionFn;
% SessionData.ssrtGrand       = ssrtGrand;
% SessionData.ssrtMean        = ssrtMean;
% SessionData.ssrtIntegration = ssrtIntegration;
% SessionData.ssrtIntegrationWeighted = ssrtIntegrationWeighted;
% SessionData.ssrtIntegrationSimple = ssrtIntegrationSimple;
% SessionData.ssrtCollapseGrand     	= ssrtCollapseGrand;
% SessionData.ssrtCollapseIntegrationWeighted     	= ssrtCollapseIntegrationWeighted;
% SessionData.ssrtCollapseIntegration     	= ssrtCollapseIntegration;
% SessionData.ssrtCollapseMean     	= ssrtCollapseMean;
% SessionData.pStopRespondGrand = pStopRespondGrand;
% SessionData.inhibitionFnGrand = inhibitionFnGrand;
% SessionData.nGo             = nGo;
% SessionData.nGoRight        = nGoRight;
% SessionData.nStopIncorrect  = nStopIncorrect;
% SessionData.nStopIncorrectRight = nStopIncorrectRight;
% SessionData.goRightLogical  = goRightLogical;
% SessionData.goRightSignalStrength = goRightSignalStrength;
% SessionData.stopRightLogical = stopRightLogical;
% SessionData.stopRightSignalStrength = stopRightSignalStrength;
% SessionData.signalStrengthLeft = signalStrengthLeft;
% SessionData.signalStrengthRight = signalStrengthRight;
% SessionData.goLeftToTarg    = goLeftToTarg;
% SessionData.goRightToTarg   = goRightToTarg;
% SessionData.goLeftToDist    = goLeftToDist;
% SessionData.goRightToDist   = goRightToDist;
% SessionData.stopLeftToTarg  = stopLeftToTarg;
% SessionData.stopRightToTarg = stopRightToTarg;
% SessionData.stopLeftToDist  = stopLeftToDist;
% SessionData.stopRightToDist = stopRightToDist;


        if printPlot
            localFigurePath = local_figure_path;
            print(figureHandle,[localFigurePath, sessionID,'_ccm_behavior'],'-dpdf', '-r300')
        end
% delete(localDataFile);



function Data = ccm_hanes_schall_95_fig7(subjectID, sessionID, options)

%
% function data = ccm_hanes_schall_95_fig7(subjectID, sessionID, options)
%
% Replicate Hanes & Schall figure 7 for Choice Countermanding task, within
% each color coherence
%
%
% If called without any arguments, returns a default options structure.
% If options are input but one is not specified, it assumes default.
%
% input:
%   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%
% Possible options are (default listed first):
%     options.collapseSignal    = Collapse across signal strength (difficulty conditions)?
%            false, true
%     options.collapseTarg 	= collapse angle/directions of the CORRECT
%     TARGET WITHIN a signal strength (so for signal strengths with correct
%     targets on the left, all left targets will be treated as one if set
%     to true
%           false, true
%     options.include50 	= if there is a 50% signal condition, do you
%           want to include it in analyses?
%           false, true
%
%     options.plotFlag       = true, false;
%     options.printPlot       = false, true;
%     options.figureHandle  = optional way to assign the figure to a handle
%
%
% Returns data structure with fields:
%
%   nGo
%   nGoRight
%   nStopIncorrect
%   nStopIncorrectRight
%   goRightLogical
%   goRightSignalStrength
%   stopRightLogical
%   stopRightSignalStrength

%%
subjectID = 'broca';
% subjectID = 'xena';
sessionID = 'xena_behavior1';
sessionID = 'broca_neural2';

ssrtMethod = 'each';
optInh              = ccm_inhibition;
optInh.plotFlag     = false;
optInh.collapseTarg = true;
Data                = ccm_inhibition(subjectID, sessionID, optInh);
nTargPair           = length(Data);

%%
plotFlag        = true;
printPlot       = false;
figureHandle    = 1998;
%{
% Set default options or return a default options structure
if nargin < 3
    options.collapseSignal   	= false;
    options.collapseTarg        = false;
    options.include50           = false;
    
    options.plotFlag            = true;
    options.printPlot           = false;
    options.figureHandle      	= 1995;
    
    % Return just the default options struct if no input
    if nargin == 0
        Data           = options;
        return
    end
end

include50       = options.include50;
plotFlag        = options.plotFlag;
printPlot       = options.printPlot;
figureHandle    = options.figureHandle;

usePreSSD = true;
useCorrectOrAll = 'all';
% ***********************************************************************
% Inhibition Function:
%       &
% SSD vs. Proportion of Response trials
%       &
% SSD vs. Proportion(Correct Choice)
% ***********************************************************************

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
ssdArray = ExtraVar.ssdArray;
pSignalArray = ExtraVar.pSignalArray;
targAngleArray = ExtraVar.targAngleArray;
nTarg = length(targAngleArray);
nTrial = size(trialData.trialOutcome, 1);



% Truncate RTs
MIN_RT = 110;
MAX_RT = 1200;
nSTD   = 3;
[allRT, outlierTrial]   = truncate_rt(trialData.rt, MIN_RT, MAX_RT, nSTD);
trialData = structfun(@(x) x(~outlierTrial,:), trialData, 'uni', false);
allRT(outlierTrial) = [];




ssdArrayRaw = trialData.stopSignalOn - trialData.responseCueOn;

% If there weren't stop trials, skip all stop-related analyses
if isempty(ssdArray)
    Data = [];
    disp('ccm_inhibition.m: No stop trials or stop trial analyses not requested');
    return
end


if ~include50
    pSignalArray(pSignalArray == .5) = [];
end

if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a choice countermanding session, try again\n')
    return
end



% Which Signal Strength levels to analyze?
switch options.collapseSignal
    case true
        nSignal = 2;
    case false
        nSignal = length(pSignalArray);
end


% If collapsing into all left and all right need to note here that there are "2" angles to deal with
% (important for calling ccm_trial_selection.m)
leftTargInd = (targAngleArray < -89) & (targAngleArray > -270) | ...
    (targAngleArray > 90) & (targAngleArray < 269);
rightTargInd = ~leftTargInd;
if options.collapseTarg
    nTargPair = 1;
else
    nTargPair = sum(rightTargInd);
    % do nothing, all target angles will be considered separately
end











%}


% Loop through targets (or collapse them if desired) and
% account for all target pairs if the session had more than one target
% pair
for kTarg = 1 : nTargPair
    ssdArray = Data(kTarg).ssdArray;
    pSignalArray = Data(kTarg).pSignalArray;
    nSSD = length(Data(kTarg).ssdArray);
    nSignal = length(Data(kTarg).pSignalArray);
    stopTargRTObs = cell(nSignal, nSSD);
    stopDistRTObs = cell(nSignal, nSSD);
    stopTargRTPrd = cell(nSignal, nSSD);
    stopDistRTPrd = cell(nSignal, nSSD);
    
    
    % ============================================================
    % Set up graphs
    if plotFlag
        figureHandle = figureHandle + 1;
        nRow = 2;
        nColumn = nSignal;
        screenOrSave = 'save';
        if printPlot
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
        else
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
        end
        choicePlotXMargin = .03;
        clf
    end
    % ============================================================
    
    
    
    for iPropIndex = 1 : nSignal
        
        
        % ============================================================
        % PLOTTING: Set up graphs
        % axes names
        axCorr = 1;
        axErr = 2;
        
        % Set up plot axes
        % inhibition function: grand (collapsed over all signal strengths)
        ax(axCorr, iPropIndex) = axes('units', 'centimeters', 'position', [xAxesPosition(axCorr, iPropIndex) yAxesPosition(axCorr, iPropIndex) axisWidth axisHeight]);
        cla
        hold(ax(axCorr, iPropIndex), 'on')
        ylim(ax(axCorr, iPropIndex), [100 500]);
          xlim(ax(axCorr, iPropIndex), [0 max(ssdArray)+20]);
       
        % inhibition function for each signal strength
        ax(axErr, iPropIndex) = axes('units', 'centimeters', 'position', [xAxesPosition(axErr, iPropIndex) yAxesPosition(axErr, iPropIndex) axisWidth axisHeight]);
        cla
        hold(ax(axErr, iPropIndex), 'on')
         ylim(ax(axErr, iPropIndex), [100 500]);
         xlim(ax(axErr, iPropIndex), [0 max(ssdArray)+20]);
       % ============================================================
        
        
        % ============================================================
        % DATA:
        
        
        % Go Trial data for this checker proportion
        % -----------------------------------
        
        % Correct and Error Go Trial RTs
        iGoTargRT = Data(kTarg).goTargRT{iPropIndex};
        iGoDistRT = Data(kTarg).goDistRT{iPropIndex};
        
        
        % Which ssrt should we use for calculating predicted stop RTs?
        switch ssrtMethod
            case 'each'
                iSSRT = Data(kTarg).ssrtIntegrationWeighted(iPropIndex);
            case 'collapsed'
                iSSRT = Data(kTarg).ssrtCollapseIntegrationWeighted;
        end
        
        % Stop Trial data for this checker proportion, within each SSD
        % -----------------------------------
        for jSSDIndex = 1 : nSSD
            ijSSD = Data(kTarg).ssdArray(jSSDIndex);
            
            stopTargRTObs{iPropIndex, jSSDIndex} = Data(kTarg).stopTargRT{iPropIndex, jSSDIndex};
            stopDistRTObs{iPropIndex, jSSDIndex} = Data(kTarg).stopDistRT{iPropIndex, jSSDIndex};
            
            stopTargRTPrd{iPropIndex, jSSDIndex} = iGoTargRT(iGoTargRT <= ijSSD + iSSRT);
            stopDistRTPrd{iPropIndex, jSSDIndex} = iGoDistRT(iGoDistRT <= ijSSD + iSSRT);
            
            
        end % iSSDIndex
        
        
        
       % Take the mean of the observed and predicted RTs for this color
       % coherence, to plot below
       iStopTargRTObs = cellfun(@nanmean, (stopTargRTObs(iPropIndex, :)));
       iStopTargRTPrd = cellfun(@nanmean, (stopTargRTPrd(iPropIndex, :)));
       iStopDistRTObs = cellfun(@nanmean, (stopDistRTObs(iPropIndex, :)));
       iStopDistRTPrd = cellfun(@nanmean, (stopDistRTPrd(iPropIndex, :)));
       % Get rid of predicted values for SSDs with not stopTarg/Dist trials
        iStopTargRTPrd(isnan(iStopTargRTObs)) = nan;
       iStopDistRTPrd(isnan(iStopDistRTObs)) = nan;
        % ============================================================
        % PLOTTING: Set up graphs
        if plotFlag
            % Determine color to use for plot based on which checkerboard color
            % proportion being used.
                 cMap = ccm_colormap(pSignalArray);
%            if options.collapseSignal
%                 cMap = ccm_colormap([0 1]);
%             else
%                 cMap = ccm_colormap(pSignalArray);
%             end
            inhColor = cMap(iPropIndex,:);
            
            
            %             plot(ax(axInhEach), ssdTimePoints, inhibitionFn{iPropIndex}, 'color', inhColor, 'linewidth', 2)
            plot(ax(axCorr, iPropIndex), Data(kTarg).ssdArray, iStopTargRTObs, '.', 'color', inhColor, 'markersize', 35)
            plot(ax(axCorr, iPropIndex), Data(kTarg).ssdArray, iStopTargRTPrd, 'o', 'color', inhColor, 'markersize', 15)

                        plot(ax(axErr, iPropIndex), Data(kTarg).ssdArray, iStopDistRTObs, '.', 'color', inhColor, 'markersize', 35)
            plot(ax(axErr, iPropIndex), Data(kTarg).ssdArray, iStopDistRTPrd, 'o', 'color', inhColor, 'markersize', 15)

        end
        
        
        
        
        
    end % iPropIndex
    
    
    
    
end % for kTarg = 1 : nTargPair

end % function













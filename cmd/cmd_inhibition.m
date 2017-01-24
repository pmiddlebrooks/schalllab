function Data = cmd_inhibition(subjectID, sessionID, options)
%%
%
% function data = ccm_inhibition(subjectID, sessionID, plotFlag, figureHandle)
%
% Response inhibition analyses for choice countermanding task.
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

% Set defaults
% Set default options or return a default options structure
if nargin < 3
    options.collapseTarg        = true;
    options.include50           = false;
    options.USE_PRE_SSD         = true;
    
    options.plotFlag            = true;
    options.printPlot           = false;
    options.figureHandle      	= 7575;
    
    % Return just the default options struct if no input
    if nargin == 0
        Data           = options;
        return
    end
end
plotFlag        = options.plotFlag;
printPlot       = options.printPlot;
figureHandle    = options.figureHandle;


% Load the data
[trialData, SessionData, ExtraVar] = load_data(subjectID, sessionID);
ssdArray        = ExtraVar.ssdArray;
timePoints                 = ssdArray(1) : ssdArray(end);
targAngleArray  = ExtraVar.targAngleArray;
nTarg           = length(targAngleArray);


if ~strcmp(SessionData.taskID, 'cmd')
    fprintf('Not a countermanding session, try again\n')
    return
end



% Truncate RTs
MIN_RT = 100;
MAX_RT = 2000;
STD_MULTIPLE    = 3;
[allRT, rtOutlierTrial] = truncate_rt(trialData.rt, MIN_RT, MAX_RT, STD_MULTIPLE);
trialData.rt(rtOutlierTrial) = nan;



if plotFlag
    nRow = 2;
    nColumn = nTarg;
    screenOrSave = 'save';
    if printPlot
        [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
    else
        [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
    end
    clf
    choicePlotXMargin = .03;
    inhColor = [0 0 0];
    stopColor = [1 0 0];
    ssdMargin = 20;
    
    % axes names
    axCollapse = nTarg + 1;
end














% ***********************************************************************
% Inhibition Function:
% ***********************************************************************

% If there weren't stop trials, skip all stop-related analyses
if plotFlag
    
    % Make one axis for each target angle
    for i = 1 : nTarg
        % left targets inhibition function
        ax(i) = axes('units', 'centimeters', 'position', [xAxesPosition(1, i) yAxesPosition(1, i) axisWidth axisHeight]);
        cla
        hold(ax(i), 'on')
    end
    
    % all targets inhibition function
    ax(axCollapse) = axes('units', 'centimeters', 'position', [mean(xAxesPosition(2, :)) yAxesPosition(2, 1) axisWidth axisHeight]);
    cla
    hold(ax(axCollapse), 'on')
end




% ----------------------------------------------------------------
% -----         Data collpased across all targets            -----

nStop               = nan(1, length(ssdArray));
nStopStop               = nan(1, length(ssdArray));
nStopTarg               = nan(1, length(ssdArray));
pStopRespond        = nan(1, length(ssdArray));


selectOpt = cmd_trial_selection;
for jSSDIndex = 1 : length(ssdArray)
    
    selectOpt.ssdRange       = ssdArray(jSSDIndex);
    
    % stop correct trials
    selectOpt.outcomeArray   = {'stopCorrect'};
    stopStopTrial         = cmd_trial_selection(trialData,  selectOpt);
    
    % stop incorrect trials
    selectOpt.outcomeArray 	= {'targetHoldAbort', 'stopIncorrect', 'stopIncorrectTarget', 'stopIncorrectPreSSDTarget'};
    stopIncorrectTrial       = cmd_trial_selection(trialData,  selectOpt);
    
    % Inhibition function data points:
    nStopStop(jSSDIndex)        	= length(stopStopTrial);
    nStopTarg(jSSDIndex)        	= length(stopIncorrectTrial);
    nStop(jSSDIndex)        	= length(stopIncorrectTrial) + length(stopStopTrial);
    pStopRespond(jSSDIndex) 	= length(stopIncorrectTrial) / nStop(jSSDIndex);
    
end % iSSDIndex


% Need go trials cumulative RTs for calculation of SSRT
selectOpt.outcomeArray   = {'goCorrectTarget','goCorrect'};
selectOpt.ssdRange       = 'none';
goTargetTrial         	= cmd_trial_selection(trialData, selectOpt);
goRT        = trialData.rt(goTargetTrial);


[fitParametersG, lowestSSE] = Weibull(ssdArray, pStopRespond, nStop);
inhibitionFn                = weibull_curve(fitParametersG, timePoints);

% SSRTs
ssrt        = get_ssrt(ssdArray, pStopRespond, nStop, goRT, fitParametersG);
        ssrtGrand = ssrt.grand;
        ssrtIntegrationWeighted = ssrt.integrationWeighted;
        ssrtIntegrationSimple = ssrt.integrationSimple;
        ssrtIntegration = ssrt.integration;
        ssrtMean = ssrt.mean;

fprintf('SSRT Collapsed:  %.2f\n', ssrt.integrationWeighted)






if plotFlag
    plot(ax(axCollapse), timePoints, inhibitionFn, 'color', inhColor, 'linewidth', 2)
    plot(ax(axCollapse), ssdArray, pStopRespond, '.k', 'markersize', 25)
    
    
    % inhibition function
    xlim(ax(axCollapse), [ssdArray(1)-ssdMargin ssdArray(end)+ssdMargin])
    set(ax(axCollapse), 'xtick', ssdArray)
    set(ax(axCollapse), 'xtickLabel', ssdArray)
    ylim(ax(axCollapse), [0 1]);
    set(get(ax(axCollapse), 'ylabel'), 'String', 'p(Respond | stop)')
    
end % if plotflag


Data.goRT               = goRT;
Data.ssrtGrand               = ssrtGrand;
Data.ssrtMean               = ssrtMean;
Data.ssrtIntegration               = ssrtIntegration;
Data.ssrtIntegrationWeighted               = ssrtIntegrationWeighted;
Data.ssrtIntegrationSimple               = ssrtIntegrationSimple;
Data.pStopRespond    = pStopRespond;
Data.ssdArray           = ssdArray;
Data.nStop              = nStop;
Data.nStopStop              = nStopStop;
Data.nStopTarg             = nStopTarg;
Data.inhibitionFn       = inhibitionFn;
Data.nTrial             = sum(nStop) + length(goTargetTrial);












% ----------------------------------------------------------------
% -----             Data separated by target                 -----

if nTarg > 1
    nStopStopAngle          = nan(nTarg, length(ssdArray));
    nStopTargAngle          = nan(nTarg, length(ssdArray));
    nStopAngle          = nan(nTarg, length(ssdArray));
    pStopRespondAngle   = nan(nTarg, length(ssdArray));
    inhibitionFnAngle   = cell(nTarg, 1);
    ssrtAngle           = cell(nTarg, 1);
    selectOpt = cmd_trial_selection;
    
    ssrtGrandAngle          = nan(nTarg, 1);
    ssrtIntegrationWeightedAngle          = nan(nTarg, 1);
    ssrtIntegrationSimpleAngle          = nan(nTarg, 1);
    ssrtIntegrationAngle          = cell(nTarg, 1);
    ssrtMeanAngle          = nan(nTarg, 1);

    
    for iTarg = 1 : nTarg
        iAngle = targAngleArray(iTarg);
        
        selectOpt.targAngle = iAngle;
        
        for jSSDIndex = 1 : length(ssdArray)
            selectOpt.ssdRange       = ssdArray(jSSDIndex);
            
            % stop correct trials
            selectOpt.outcomeArray   = {'stopCorrect'};
            stopStopTrialAngle    = cmd_trial_selection(trialData,  selectOpt);
            
            % stop incorrect trials
            selectOpt.outcomeArray	= {'targetHoldAbort', 'stopIncorrect', 'stopIncorrectTarget', 'stopIncorrectPreSSDTarget'};
            stopIncorrectTrialAngle  = cmd_trial_selection(trialData,  selectOpt);
            
            nStopStopAngle(iTarg, jSSDIndex)               = length(stopStopTrialAngle);
            nStopTargAngle(iTarg, jSSDIndex)               = length(stopIncorrectTrialAngle);
            nStopAngle(iTarg, jSSDIndex)               = length(stopStopTrialAngle) + length(stopIncorrectTrialAngle);
            pStopRespondAngle(iTarg, jSSDIndex)     = length(stopIncorrectTrialAngle) / nStopAngle(iTarg, jSSDIndex);
            
        end % iSSDIndex
        
        
        
        % Need go trials cumulative RTs for calculation of SSRT
        selectOpt.outcomeArray   = {'goCorrectTarget','goCorrect'};
        selectOpt.ssdRange       = 'none';
        iGoTargetTrial         	= cmd_trial_selection(trialData, selectOpt);
        iGoRT   = trialData.rt(iGoTargetTrial);
        
        
        % % Inhibition function data
        % keepSSDR                = ~isnan(pStopRespondAngle(iTarg, jSSDIndex));
        % pStopRespondAngle    = pStopRespondAngle(iTarg, keepSSDR);
        % nStopAngle              = nStopAngle(iTarg, keepSSDR);
        % ssdArrayRight           = ssdArray(keepSSDR);
        
        % Inhibition functions
        [fitParameters, lowestSSE] = Weibull(ssdArray, pStopRespondAngle(iTarg, :), nStopAngle(iTarg, :));
        inhibitionFnAngle{iTarg}           = weibull_curve(fitParameters, timePoints);
        
        % SSRTs
        ssrtAngle{iTarg}   = get_ssrt(ssdArray, pStopRespondAngle(iTarg, :), nStopAngle(iTarg, :), iGoRT, fitParameters);
        ssrtGrandAngle(iTarg) = ssrtAngle{iTarg}.grand;
        ssrtIntegrationWeightedAngle(iTarg) = ssrtAngle{iTarg}.integrationWeighted;
        ssrtIntegrationSimpleAngle(iTarg) = ssrtAngle{iTarg}.integrationSimple;
        ssrtIntegrationAngle{iTarg} = ssrtAngle{iTarg}.integration;
        ssrtMeanAngle(iTarg) = ssrtAngle{iTarg}.mean;
        
        
        fprintf('SSRTs:  Angle %d:\t%.2f\n', targAngleArray(iTarg), ssrtAngle{iTarg}.integrationWeighted)
        
        
        
        
        
        
        
        if plotFlag
            plot(ax(iTarg), timePoints, inhibitionFnAngle{iTarg}, 'color', inhColor, 'linewidth', 2)
            plot(ax(iTarg), ssdArray, pStopRespondAngle(iTarg, :), '.k', 'markersize', 25)
            
            % inhibition function
            xlim(ax(iTarg), [ssdArray(1)-ssdMargin ssdArray(end)+ssdMargin])
            set(ax(iTarg), 'xtick', ssdArray)
            set(ax(iTarg), 'xtickLabel', ssdArray)
            ylim(ax(iTarg), [0 1]);
            yL = sprintf('%d deg  p(Respond | stop)', targAngleArray(iTarg));
            set(get(ax(iTarg), 'ylabel'), 'String', yL)
            
        end % if plotflag
        
    end % iTarg
    
    
Data.ssrtGrandAngle               = ssrtGrandAngle;
Data.ssrtMeanAngle               = ssrtMeanAngle;
Data.ssrtIntegrationAngle               = ssrtIntegrationAngle;
Data.ssrtIntegrationWeightedAngle               = ssrtIntegrationWeightedAngle;
Data.ssrtIntegrationSimpleAngle               = ssrtIntegrationSimpleAngle;
    Data.pStopRespondAngle    = pStopRespondAngle;
    Data.nStopAngle              = nStopAngle;
    Data.nStopStopAngle              = nStopStopAngle;
    Data.nStopTargAngle              = nStopTargAngle;
    Data.inhibitionFnAngle       = inhibitionFnAngle;
    
end





if printPlot
    print(gcf, ['~/schalllab/tempfigures/',sessionID, '_Behavior'], '-dpdf')
end




% Data.ssd                = ssd;
% Data.goRTMinusSSD       = goRTMinusSSD;
% Data.nStopStop          = nStopStop;
% Data.nStopTarg          = nStopTarg;
% Data.nStopDist          = nStopDist;
% Data.stopTargetProb     = stopTargetProb;
% % Data.ssrt               = ssrt;
% % Data.ssrtMean               = ssrtMean;
% % Data.ssrtIntegration               = ssrtIntegration;
% Data.pStopRespondGrand = pStopRespondGrand;
% Data.inhibitionFnGrand  = inhibitionFnGrand;
% % Data.ssrtGrand          = ssrtGrand;
% % Data.ssrtMeanGrand          = ssrtMeanGrand;
% % Data.ssrtIntegrationGrand          = ssrtIntegrationGrand;
%
% Data.stopRespondRTPredictGrand = stopRespondRTPredictGrand;
% Data.stopRespondRTObserveGrand = stopRespondRTObserveGrand;
% Data.stopRespondRTPredict = stopRespondRTPredict;
% Data.stopRespondRT = stopRespondRT;
%
%

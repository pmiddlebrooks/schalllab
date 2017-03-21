function Data = ccm_inhibition_rt(subjectID, sessionID, options)

%
% function data = ccm_inhibition(subjectID, sessionID, plotFlag, figureHandle)
%
% Response inhibition analyses for choice countermanding task.
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
% sessionID = 'bp093n02';
% Set default options or return a default options structure
if nargin < 3
    options.collapseSignal   	= false;
    options.collapseTarg        = true;
    options.include50           = false;
    options.USE_PRE_SSD         = true;
    options.USE_TWO_COLORS         = false;
    
    options.plotFlag            = true;
    options.printPlot           = true;
    options.figureHandle      	= 333;
    
    % Return just the default options struct if no input
    if nargin == 0
        Data           = options;
        return
    end
end

plotFlag        = options.plotFlag;
printPlot       = options.printPlot;
figureHandle    = options.figureHandle;



optInh                  = options;
optInh.plotFlag         = false;
optInh.printPlot        = false;
optInh.collapseSignal   = options.collapseSignal;
Data                    = ccm_inhibition(subjectID, sessionID, optInh);

% Which Signal Strength levels to analyze
switch options.collapseSignal
    case true
        nSignal = 2;
    case false
        nSignal = length(Data.pSignalArray);
end

% Plotting parameters:
inhColor = [0 0 0];
inhLineW = 2;
goLineW = 4;
stopColor = [1 0 0];
stopLineW = 4;

switch subjectID
    case 'xena'
        xlimSubj = [50 600];
    otherwise
        xlimSubj = [50 850];
end


stepSize = .01;
nStep = 1 / stepSize;
rtDiff = cell(nSignal, 1);
rtDiffTarg = cell(nSignal, 1);



ssdTimePoints = Data.ssdArray(1) : Data.ssdArray(end);
goRT            = cell(nSignal, 1);
stopRT          = cell(nSignal, 1);
propGoRT        = cell(nSignal, 1);
propStopRT      = cell(nSignal, 1);
goTargRT        = cell(nSignal, 1);
stopTargRT      = cell(nSignal, 1);
propGoTargRT    = cell(nSignal, 1);
propStopTargRT  = cell(nSignal, 1);

if options.collapseSignal
    cMap = ccm_colormap([0 1]);
else
    cMap = ccm_colormap(Data.pSignalArray);
end



if plotFlag
    nRow = 3;
    nColumn = length(Data.pSignalArray);
    if printPlot
        [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
    else
        [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
    end
    clf
end

for iCoh = 1 : nSignal
    
    iGoColor = cMap(iCoh,:);
    %     if Data.pSignalArray(iCoh) < .5
    %         iGoColor = cMap(1,:);
    %     else
    %         iGoColor = cMap(end,:);
    %     end
    
    % Set up plot axes
    % Ond=e plot for each choice difficulty
    if plotFlag
        
        % Correct error GO and STOP RT plots
        ax(1, iCoh) = axes('units', 'centimeters', 'position', [xAxesPosition(1, iCoh) yAxesPosition(1, iCoh) axisWidth axisHeight]);
        cla
        hold(ax(1, iCoh), 'on')
        set(ax(1, iCoh), 'Xlim', xlimSubj)
        
        % Correct-only GO and STOP RT plots
        ax(2, iCoh) = axes('units', 'centimeters', 'position', [xAxesPosition(2, iCoh) yAxesPosition(2, iCoh) axisWidth axisHeight]);
        cla
        hold(ax(2, iCoh), 'on')
        set(ax(2, iCoh), 'Xlim', xlimSubj)
        
        if iCoh == 1
            ylabel(ax(1, iCoh), 'All RTs')
            ylabel(ax(2, iCoh), 'Correct RTs')
        end
        if iCoh > 1
            set(ax(1, iCoh), 'yticklabel', [])
            set(ax(2, iCoh), 'yticklabel', [])
        end
        
    end
    
    % CORRECT AND ERROR CHOICE RTS
    
    % Cumulative GO RTs
    goRT{iCoh} = sort(nonnans(Data.goTotalRT{iCoh, 1}));
    propGoRT{iCoh} = cmtb_edf(goRT{iCoh},goRT{iCoh});
    
    % Cumulative STOP RTs
    stopRT{iCoh} = sort(nonnans([cell2mat(Data.stopTargRT(iCoh,:)');cell2mat(Data.stopDistRT(iCoh,:)')]));
    propStopRT{iCoh} = cmtb_edf(stopRT{iCoh},stopRT{iCoh});
    % Make STOP RTs defective w.r.t. GO RTs
    propStopRT{iCoh} =  propStopRT{iCoh} * (length(stopRT{iCoh}) / (length(stopRT{iCoh}) + length(goRT{iCoh})));
    
    
    
    
    % CORRECT CHOICE RTS
    
    % Cumulative GO RTs
    goTargRT{iCoh} = sort(nonnans(Data.goTargRT{iCoh, 1}));
    propGoTargRT{iCoh} = cmtb_edf(goTargRT{iCoh},goTargRT{iCoh});
    
    % Cumulative STOP RTs
    stopTargRT{iCoh} = sort(nonnans(cell2mat(Data.stopTargRT(iCoh,:)')));
    propStopTargRT{iCoh} = cmtb_edf(stopTargRT{iCoh},stopTargRT{iCoh});
    % Make STOP RTs defective w.r.t. GO RTs
    propStopTargRT{iCoh} =  propStopTargRT{iCoh} * (length(stopTargRT{iCoh}) / (length(stopTargRT{iCoh}) + length(goTargRT{iCoh})));
    if plotFlag
        % Plot Cumulative GO RTs
        plot(ax(1, iCoh), goRT{iCoh}, propGoRT{iCoh}, 'color', iGoColor, 'linewidth', goLineW)
        % Plot Cumulative STOP RTs
        plot(ax(1, iCoh), stopRT{iCoh}, propStopRT{iCoh}, 'color', stopColor, 'linewidth', stopLineW)
        % Plot Inhibition Function
        plot(ax(1, iCoh), Data.ssdArray, Data.stopRespondProb(iCoh,:), '.k', 'markersize', 25);
        plot(ax(1, iCoh), ssdTimePoints, Data.inhibitionFn{iCoh}, 'color', inhColor, 'linewidth', inhLineW)
        
        
        
        % Plot Cumulative GO RTs
        plot(ax(2, iCoh), goTargRT{iCoh}, propGoTargRT{iCoh}, 'color', iGoColor, 'linewidth', goLineW)
        % Plot Cumulative STOP RTs
        plot(ax(2, iCoh), stopTargRT{iCoh}, propStopTargRT{iCoh}, 'color', stopColor, 'linewidth', stopLineW)
        % Plot Inhibition Function
        plot(ax(2, iCoh), Data.ssdArray, Data.stopRespondProb(iCoh,:), '.k', 'markersize', 25);
        plot(ax(2, iCoh), ssdTimePoints, Data.inhibitionFn{iCoh}, 'color', inhColor, 'linewidth', inhLineW)
    end
    
    
    
    
    % goRT - stopRT cumulative distribution Fns
    iRtDiff = nan(nStep, 1);
    iRtDiffTarg = nan(nStep, 1);
    for i = stepSize : stepSize : 1
        iGoInd = ceil(length(goRT{iCoh}) * i);
        iStopInd = ceil(length(stopRT{iCoh}) * i);
        
        iRtDiff(round(i * 1/stepSize)) = goRT{iCoh}(iGoInd) - stopRT{iCoh}(iStopInd);
        
        iGoTargInd = ceil(length(goTargRT{iCoh}) * i);
        iStopTargInd = ceil(length(stopTargRT{iCoh}) * i);
        
        if iStopTargInd > 0
            iRtDiffTarg(round(i * 1/stepSize)) = goTargRT{iCoh}(iGoTargInd) - stopTargRT{iCoh}(iStopTargInd);
        end
    end
    rtDiff{iCoh} = iRtDiff;
    rtDiffTarg{iCoh} = iRtDiffTarg;
    
end




if plotFlag
    % LEFT COLOR TRIALS: correct and error choices:  goRT - stopRT cumulative distribution Fns
    ax(3, 1) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 1) yAxesPosition(3, 1) axisWidth axisHeight]);
    hold(ax(3, 1), 'on')
    set(ax(3, 1), 'ylim', [0 1], 'xlim', [-100 350])
    
    % RIGHT COLOR TRIALS: correct and error choices: goRT - stopRT cumulative distribution Fns
    ax(3, 2) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 2) yAxesPosition(3, 2) axisWidth axisHeight]);
    hold(ax(3, 2), 'on')
    set(ax(3, 2), 'ylim', [0 1], 'xlim', [-100 350])
    
    % LEFT COLOR TRIALS: correct-only choices:  goRT - stopRT cumulative distribution Fns
    ax(3, 3) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 3) yAxesPosition(3, 3) axisWidth axisHeight]);
    hold(ax(3, 3), 'on')
    set(ax(3, 3), 'ylim', [0 1], 'xlim', [-100 350])
    
    % RIGHT COLOR TRIALS: correct-only choices; goRT - stopRT cumulative distribution Fns
    ax(3, 4) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 4) yAxesPosition(3, 4) axisWidth axisHeight]);
    hold(ax(3, 4), 'on')
    set(ax(3, 4), 'ylim', [0 1], 'xlim', [-100 350])
end
for iCoh = 1 : nSignal/2
    iPropIndexL = nSignal/2 + 1 - iCoh;
    iPropIndexR = iCoh + nSignal/2;  % Reverse order of plotting to keep color overlays similar between left and right target
    
    if plotFlag
        plot(ax(3, 1), rtDiff{iPropIndexL}, [stepSize : stepSize : 1], '-', 'color', cMap(iPropIndexL,:), 'lineWidth', 4)
        plot(ax(3, 2), rtDiff{iPropIndexR}, [stepSize : stepSize : 1], '-', 'color', cMap(iPropIndexR,:), 'lineWidth', 4)
        plot(ax(3, 3), rtDiffTarg{iPropIndexL}, [stepSize : stepSize : 1], '-', 'color', cMap(iPropIndexL,:), 'lineWidth', 4)
        plot(ax(3, 4), rtDiffTarg{iPropIndexR}, [stepSize : stepSize : 1], '-', 'color', cMap(iPropIndexR,:), 'lineWidth', 4)
    end
end

if plotFlag
    plot(ax(3, 1), [0 0], [0 1], '--k')
    ylabel(ax(3, 1), 'LEFT: goRT - StopRT ')
    
    plot(ax(3, 2), [0 0], [0 1], '--k')
    ylabel(ax(3, 2), 'RIGHT: goRT - StopRT ')
    set(ax(3, 2), 'yticklabel', [])
    
    plot(ax(3, 3), [0 0], [0 1], '--k')
    ylabel(ax(3, 3), 'Correct LEFT: goRT - StopRT ')
    set(ax(3, 3), 'yticklabel', [])
    
    plot(ax(3, 4), [0 0], [0 1], '--k')
    ylabel(ax(3, 4), 'Correct RIGHT: goRT - StopRT ')
    set(ax(3, 4), 'yticklabel', [])
    
    
    
    h=axes('Position', [0 0 1 1], 'Visible', 'Off');
    titleString = sprintf('\n\n%s \t %s', subjectID,sessionID);
    text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')
    
    
end















if printPlot && plotFlag
    print(figureHandle,fullfile(local_figure_path, subjectID, [sessionID,'_ccm_inhibition_rt.pdf']),'-dpdf', '-r300')
end


Data.stopRT         = stopRT;
Data.goRT           = goRT;
Data.propStopRT     = propStopRT;
Data.propGoRT       = propGoRT;
Data.stopTargRT     = stopTargRT;
Data.goTargRT       = goTargRT;
Data.propStopTargRT = propStopTargRT;
Data.propGoTargRT   = propGoTargRT;
end % function













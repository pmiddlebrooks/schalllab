function Data = ccm_nelson_etal_2010_fig8(subjectID, sessionID, options)



%%
selectionVar = {'trialOutcome', 'targ1CheckerProp', 'ssd', 'targAngle', 'saccToTargIndex', 'saccAngle', 'rt'};
tdSelect = trialData(:, selectionVar);


% no-stop  -  no-stop  -  no-stop

% Get default options structure:
selectOpt = ccm_trial_selection;

seqOpt              = selectOpt;
seqOpt(1).outcome   = {'goCorrectTarget'};
seqOpt(1).ssd       = 'none';

seqOpt(2)           = selectOpt;
seqOpt(2).outcome   = {'goCorrectTarget'};
seqOpt(2).ssd       = 'none';

seqOpt(3)           = selectOpt;
seqOpt(3).outcome   = {'goCorrectTarget'};
seqOpt(3).ssd       = 'none';

goGoGo              = ccm_trial_sequence(tdSelect, seqOpt);



% no-stop  -  canceled  -  no-stop

% Get default options structure:
selectOpt = ccm_trial_selection;

seqOpt              = selectOpt;
seqOpt(1).outcome   = {'goCorrectTarget'};
seqOpt(1).ssd       = 'none';

seqOpt(2)           = selectOpt;
seqOpt(2).outcome   = {'stopCorrect'};

seqOpt(3)           = selectOpt;
seqOpt(3).outcome   = {'goCorrectTarget'};
seqOpt(3).ssd       = 'none';

goCancGo            = ccm_trial_sequence(tdSelect, seqOpt);



% no-stop  -  noncanceled  -  no-stop

% Get default options structure:
selectOpt = ccm_trial_selection;

seqOpt              = selectOpt;
seqOpt(1).outcome   = {'goCorrectTarget'};
seqOpt(1).ssd       = 'none';

seqOpt(2)           = selectOpt;
seqOpt(2).outcome   = {'stopIncorrectTarget'};

seqOpt(3)           = selectOpt;
seqOpt(3).outcome   = {'goCorrectTarget'};
seqOpt(3).ssd       = 'none';

goNonCancGo         = ccm_trial_sequence(tdSelect, seqOpt);

%% Set up plotting
plotFlag = true;
printPlot = true;
    if plotFlag
        figureHandle = 84;
        nRow = 2;
        nColumn = 2;
        screenOrSave = 'save';
        if printPlot
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
        else
            [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
        end
        clf
        choicePlotXMargin = .03;
        ssdMargin = 20;
        ylimArray = [];
        
        % axes names
        axTriplet = 1;
%         axInhEach = 2;
%         axInhRTSSD = 6;
%         SSDvPCorrect = 3;
%         ssrtPRight = 4;
%         SSDvSigStrength = 5;
        
        
        % Set up plot axes
        % inhibition function: grand (collapsed over all signal strengths)
        ax(axTriplet) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
        cla
        hold(ax(axTriplet), 'on')
        
%         % inhibition function for each signal strength
%         ax(axInhEach) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
%         cla
%         hold(ax(axInhEach), 'on')
%         
%         % goRTMean - SSD inhibition function for each signal strength
%         ax(axInhRTSSD) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 3) yAxesPosition(1, 3) axisWidth axisHeight]);
%         cla
%         hold(ax(axInhRTSSD), 'on')
%         
%         % p(right) vs ssrt
%         ax(ssrtPRight) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 3) yAxesPosition(2, 3) axisWidth axisHeight]);
%         cla
%         hold(ax(ssrtPRight), 'on')
%         
%         % SSD vs p(correct)
%         ax(SSDvPCorrect) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 1) yAxesPosition(3, 1) axisWidth axisHeight]);
%         cla
%         hold(ax(SSDvPCorrect), 'on')
%         
%         % SSD vs p(correct)
%         ax(SSDvSigStrength) = axes('units', 'centimeters', 'position', [xAxesPosition(3, 3) yAxesPosition(3, 3) axisWidth axisHeight]);
%         cla
%         hold(ax(SSDvSigStrength), 'on')
    end
 

    plot(ax(axTriplet), [1 2], [tdSelect.rt(goGoGo) tdSelect.rt(goGoGo+2)], '-k')
    plot(ax(axTriplet), [3 4], [tdSelect.rt(goCancGo) tdSelect.rt(goCancGo+2)], '-r')
    plot(ax(axTriplet), [5 6], [tdSelect.rt(goNonCancGo) tdSelect.rt(goNonCancGo+2)], '--r')
    
    [h,p,ci,stats] = ttest2(tdSelect.rt(goGoGo), tdSelect.rt(goGoGo+2))
    [h,p,ci,stats] = ttest2(tdSelect.rt(goCancGo), tdSelect.rt(goCancGo+2))
    [h,p,ci,stats] = ttest2(tdSelect.rt(goNonCancGo), tdSelect.rt(goNonCancGo+2))
    
    mean([tdSelect.rt(goGoGo), tdSelect.rt(goGoGo+2)])
    mean([tdSelect.rt(goCancGo), tdSelect.rt(goCancGo+2)])
    mean([tdSelect.rt(goNonCancGo), tdSelect.rt(goNonCancGo+2)])

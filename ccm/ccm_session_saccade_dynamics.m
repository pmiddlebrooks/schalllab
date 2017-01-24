function Data = ccm_session_saccade_dynamics(subject, session, Options)

if nargin < 3
    Options = ccm_options;
    if nargin == 0
        Data           = Options;
        return
    end
end
%%
% Options = ccm_options;
% subject = 'broca';
% session = 'bp093n02';

% Load the data
[trialData, SessionData, ExtraVar] = load_data(subject, session);
colorCohArray    = ExtraVar.pSignalArray;
targAngleArray	= ExtraVar.targAngleArray;
ssdArray        = ExtraVar.ssdArray;
nSSD            = length(ssdArray);
colorCohArray(colorCohArray == .5) = [];

if ~strcmp(SessionData.taskID, 'ccm')
    fprintf('Not a chioce countermanding saccade session, try again\n')
    return
end

% Make color coherence array a 2-row matrix of left and right values
origCohArray = colorCohArray;
colorCohArray = reshape(colorCohArray, length(colorCohArray)/2, 2)';
colorCohArray(1,:) = fliplr(colorCohArray(1,:));
%%







% Initialize cells (direction X color coherence)
goTrial = cell(2,size(colorCohArray, 2));
goBegin = cell(2,size(colorCohArray, 2));
goEnd = cell(2,size(colorCohArray, 2));
goEyeX = cell(2,size(colorCohArray, 2));
goEyeY = cell(2,size(colorCohArray, 2));
goDeg = cell(2,size(colorCohArray, 2));
goVel = cell(2,size(colorCohArray, 2));


stopTrial = cell(2,size(colorCohArray, 2));
stopBegin = cell(2,size(colorCohArray, 2));
stopEnd = cell(2,size(colorCohArray, 2));
stopEyeX = cell(2,size(colorCohArray, 2));
stopEyeY = cell(2,size(colorCohArray, 2));
stopDeg = cell(2,size(colorCohArray, 2));
stopVel = cell(2,size(colorCohArray, 2));








%   Get relevant trials for each condition
% ___________________________________________

Opt = ccm_options;

% Loop through color coherence values
for i = 1 : size(colorCohArray, 2)
    
    
    % Loop through directions (left to right)
    for d = 1 : 2;
        
        
        Opt.rightCheckerPct = colorCohArray(d,i)*100;
        
        % Go trials
        % -------------------------
        Opt.outcome = {'goCorrectTarget'};
        Opt.ssd = 'none';
        
        goTrial{d, i} = ccm_trial_selection(trialData, Opt);
        goTrial{d, i}(isnan(trialData.saccToTargIndex(goTrial{d, i}))) = [];

        % Get beginning and ending times of saccades
        goBegin{d, i} = cellfun(@(x,y) x(y), trialData.saccBegin(goTrial{d, i}), num2cell(trialData.saccToTargIndex(goTrial{d, i})), 'uni', false);
        goEnd{d, i} = cellfun(@(x,y,z) x+y(z), goBegin{d, i}, trialData.saccDuration(goTrial{d, i}), num2cell(trialData.saccToTargIndex(goTrial{d, i})), 'uni', false);
        
        % Get horiz(X) and vert (Y) components of saccades
        goEyeX{d, i} = cellfun(@(x,y,z) x(y:z), trialData.eyeX(goTrial{d, i}), goBegin{d, i}, goEnd{d, i}, 'uni', false);
        goEyeY{d, i} = cellfun(@(x,y,z) x(y:z), trialData.eyeY(goTrial{d, i}), goBegin{d, i}, goEnd{d, i}, 'uni', false);
        
        % Use pythagorean theorom to calculate saccade tragjectory in degrees
        goDeg{d, i} = cellfun(@(x,y) sqrt(x.^2 + y.^2), goEyeX{d, i}, goEyeY{d, i}, 'uni', false);
        if d == 1
            goDeg{d, i} = cellfun(@(x) -x, goDeg{d, i}, 'uni', false);
        end
        goVel{d, i} = cellfun(@(x) [0; abs(diff(x(:)))], goDeg{d, i}, 'uni', false);
        
        
        
        % Stop trials
        % -------------------------
        Opt.outcome = {'stopIncorrectTarget', 'targetHoldAbort'};
        Opt.ssd = 'collapse'; % collapse SSDs
        
        stopTrial{d, i} = ccm_trial_selection(trialData, Opt);
        stopTrial{d, i}(isnan(trialData.saccToTargIndex(stopTrial{d, i}))) = [];
        
        % Get beginning and ending times of saccades
        stopBegin{d, i} = cellfun(@(x,y) x(y), trialData.saccBegin(stopTrial{d, i}), num2cell(trialData.saccToTargIndex(stopTrial{d, i})), 'uni', false);
        stopEnd{d, i} = cellfun(@(x,y,z) x+y(z), stopBegin{d, i}, trialData.saccDuration(stopTrial{d, i}), num2cell(trialData.saccToTargIndex(stopTrial{d, i})), 'uni', false);
        
        % Get horiz(X) and vert (Y) components of saccades
        stopEyeX{d, i} = cellfun(@(x,y,z) x(y:z), trialData.eyeX(stopTrial{d, i}), stopBegin{d, i}, stopEnd{d, i}, 'uni', false);
        stopEyeY{d, i} = cellfun(@(x,y,z) x(y:z), trialData.eyeY(stopTrial{d, i}), stopBegin{d, i}, stopEnd{d, i}, 'uni', false);
        
        % Use pythagorean theorom to calculate saccade tragjectory in degrees
        stopDeg{d, i} = cellfun(@(x,y) sqrt(x.^2 + y.^2), stopEyeX{d, i}, stopEyeY{d, i}, 'uni', false);
        if d == 1
            stopDeg{d, i} = cellfun(@(x) -x, stopDeg{d, i}, 'uni', false);
        end
        stopVel{d, i} = cellfun(@(x) [0; abs(diff(x(:)))], stopDeg{d, i}, 'uni', false);
        
        
        
    end
end










%      set up plot
% _____________________________
if Options.plotFlag
    figureHandle = 6;
    
    % Get colors for difference color coherence values
    cMap = ccm_colormap(origCohArray);
    % For this function, flip part of the map upside down (to start with
    % hard colors for both response sides and end with easy colors)
    cMap(1:size(colorCohArray, 2), :) = flipud(cMap(1:size(colorCohArray, 2), :));
    
    markerSize = 10;
    lineW = 2;
    
    nCol = 2;
    nRow = 4;
    colGo = 1;
    colStop = 2;
    rowTrace = 1;
    rowVelo = 2;
    rowRT = 3;
    rowAmp = 4;
    
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nRow, nCol, 'portrait', figureHandle);
    clf
    h=axes('Position', [0 0 1 1], 'Visible', 'Off');
    set(gcf, 'Name', 'Saccade Dynamics','NumberTitle', 'off')
    %         titleString = sprintf('%s \t %s', sessionID, Unit(kUnitIndex).name);
    %         text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')
    
    
    xMaxTrace = max([abs(goEyeX{1,1}{1}(end)) abs(goEyeX{2,1}{1}(end))]) * 1.2;
    yMaxTrace = xMaxTrace;
    
    yMaxVelo = max([max(goVel{1,1}{1}) max(goVel{2,1}{1})]) * 1.3;
    yMinVelo = min([min(goVel{1,1}{1}) min(goVel{2,1}{1})]) * .7;
    
    % Create axes
    % Loop through color coherence values
    for j = 1 : nCol
        
        %      Raw eye traces to targets
        % _____________________________
        ax(rowTrace, j) = axes('units', 'centimeters', 'position', [xAxesPosition(rowTrace, j) yAxesPosition(rowTrace, j) axisWidth axisHeight]);
        set(ax(rowTrace, j), 'ylim', [-yMaxTrace yMaxTrace], 'xlim', [-xMaxTrace xMaxTrace])
        cla
        hold(ax(rowTrace, j), 'on')
        
        %      Velocity as a function of degrees
        % _____________________________
        ax(rowVelo, j) = axes('units', 'centimeters', 'position', [xAxesPosition(rowVelo, j) yAxesPosition(rowVelo, j) axisWidth axisHeight]);
        set(ax(rowVelo, j), 'ylim', [0 yMaxVelo], 'xlim', [-xMaxTrace xMaxTrace])
        cla
        hold(ax(rowVelo, j), 'on')
        
        %      Peak velocity as a function of degrees
        % _____________________________
        ax(rowRT, j) = axes('units', 'centimeters', 'position', [xAxesPosition(rowRT, j) yAxesPosition(rowRT, j) axisWidth axisHeight]);
        %                 set(ax(rowRT, j), 'ylim', [yMinVelo yMaxVelo], 'xlim', [-xMaxTrace xMaxTrace])
        cla
        hold(ax(rowRT, j), 'on')
        
        
    end
    
    
    
    
    % Loop through color coherence values
    for i = 1 : size(colorCohArray, 2)
        
        % color index to draw correct color from cMap
        colorInd = i;
        
        % Loop through directions (left to right)
        for d = 1 : 2;
            
            % color index to draw correct color from cMap
            if d == 2
                colorInd = i + size(colorCohArray, 2);
            end
            
            %      Raw eye traces to targets
            % _____________________________
            axes(ax(rowTrace, colGo))
            cellfun(@(x,y) plot(x,y, 'color', cMap(colorInd,:)), goEyeX{d, i}, goEyeY{d, i})
            axes(ax(rowTrace, colStop))
            cellfun(@(x,y) plot(x,y, 'color', cMap(colorInd,:)), stopEyeX{d, i}, stopEyeY{d, i})
            
            %      Velocity as a function of degrees
            % _____________________________
            axes(ax(rowVelo, colGo))
            cellfun(@(x,y) plot(x,y, 'color', cMap(colorInd,:)), goDeg{d, i}, goVel{d, i})
            axes(ax(rowVelo, colStop))
            cellfun(@(x,y) plot(x,y, 'color', cMap(colorInd,:)), stopDeg{d, i}, stopVel{d, i})
            
            %      Peak velocity as a function of degrees
            % _____________________________
            axes(ax(rowRT, colGo))
%             plot(trialData.rt(goTrial{d, i}), cellfun(@max, goVel{d, i}),'.', 'markerSize', markerSize, 'markerFacecolor', cMap(colorInd,:), 'markerEdgecolor', cMap(colorInd,:))
            plot(trialData.rt(goTrial{d, i}), cellfun(@max, goVel{d, i}),'.', 'markerSize', markerSize, 'color', cMap(colorInd,:))
            axes(ax(rowRT, colStop))
%             plot(trialData.rt(stopTrial{d, i}), cellfun(@max, stopVel{d, i}),'.', 'markerSize', markerSize, 'markerFacecolor', cMap(colorInd,:), 'markerEdgecolor', cMap(colorInd,:))
            plot(trialData.rt(stopTrial{d, i}), cellfun(@max, stopVel{d, i}),'.', 'markerSize', markerSize, 'color', cMap(colorInd,:))
            
            
            
        end
    end
    
    
    if Options.printPlot
        print(figureHandle,fullfile(local_figure_path, subject, [session, '_ccm_saccade_dynamics.pdf']),'-dpdf', '-r300')
    end
    
end
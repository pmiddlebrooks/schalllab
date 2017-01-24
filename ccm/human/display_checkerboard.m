function display_checkerboard

magentaTargSide = 'right';

secPerStim = 5;
squareWidthPixel         = 2; % pixels
    choiceStimulusAmplitude = 0;
    choiceStimulusAngle = 0;

% rightTargPropArray   = [0 1 ];
rightTargPropArray     = [.41 .45 .48 .5 .52 .55 .59];
rightTargetRate = .5; % How often should right side be target?
fiftyPercentRate = .6; % how often should 50% signal strength be presented RELATIVE TO OTHER proportions?

% Checkered Stimulus constants
nCheckerColumn     = 20;
nCheckerRow        = 10;
nSquare            = nCheckerColumn * nCheckerRow;
pixelPerDegree = 40;
squareWidth = squareWidthPixel / pixelPerDegree;

% Initialize variables
stimulusSquaresArray    = zeros(4, nSquare);
iStimulusColorsArray    = zeros(3, nSquare);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Set up Experiment Variables, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KbName('UnifyKeyNames');
stopKey=KbName('ESCAPE');

% Get distance values and conversions for stimuli presentation
metersFromScreen    = .66;
screenWidthMeters   = .4;
theta               = asind(screenWidthMeters / 2 / sqrt(((screenWidthMeters / 2)^2) + metersFromScreen^2));
screenSize          = get(0, 'ScreenSize');
screenWidthPixel  	= screenSize(3);
screenHeightPixel 	= screenSize(4);
matlabCenterX       = screenWidthPixel/2;
matlabCenterY       = screenHeightPixel/2;
pixelPerDegree     = screenWidthPixel / (2*theta);

stimulusAmplitudePixel = choiceStimulusAmplitude * pixelPerDegree;



whichScreen     = 0;
backGround      = [60 60 60];
[window, centerPoint] = Screen('OpenWindow', whichScreen, backGround);
priorityLevel   = MaxPriority(window);
Priority(priorityLevel);

commandwindow;
flipTime        = Screen('GetFlipInterval',window);% get the flip rate of current monitor.
preFlipBuffer   = flipTime / 2;
% a='Resolution & Refresh ok';
% if timing>.012 || timing<.011 %make sure we're running in 1024*768 at 85 Hz, else stop
%     clear a;
%     display('Please change screen to 1024x768 and 85Hz');
% end
% display a;
% frameFrequency = 85;
frameFrequency = 1 / flipTime;


% ***************************************************************
%     Generate a new random number stream based on clock
% ***************************************************************
%        Replace the default stream with a stream whose seed is based on CLOCK, so
%        RAND will return different values in different MATLAB sessions.  NOTE: It
%        is usually not desirable to do this more than once per MATLAB
%        session.
s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setDefaultStream(s);







% ------------------------------
%     CHECKERED STIMULUS
% ------------------------------
% First determine the center of the 10X10 stimulus in x,y coordinates w.r.t the center of the
% screen
stimulusEyeLinkX     = stimulusAmplitudePixel * cosd(choiceStimulusAngle);
stimulusEyeLinkY     = stimulusAmplitudePixel * sind(choiceStimulusAngle);
stimulusMatlabX     = matlabCenterX + stimulusEyeLinkX;
stimulusMatlabY     = matlabCenterY - stimulusEyeLinkY;
stimulusWindowPixel      = squareWidthPixel * nCheckerColumn; % For now, make it window the same size as the stimulus
choiceStimulusWindow = stimulusWindowPixel / pixelPerDegree;
stimulusWindow       = [stimulusEyeLinkY - stimulusWindowPixel/2, stimulusEyeLinkX - stimulusWindowPixel/2, stimulusEyeLinkY + stimulusWindowPixel/2, stimulusEyeLinkX + stimulusWindowPixel/2];
% Get the positions of each checker square
for iRow = 1 : nCheckerRow
    for iColumn = 1 : nCheckerColumn
        stimulusIndex   = iColumn + (nCheckerColumn * (iRow - 1));
        % Get the centers of each square
        iSquareCenterX  = stimulusMatlabX + squareWidthPixel * ((iColumn-1) - (nCheckerColumn/2 - .5));
        iSquareCenterY  = stimulusMatlabY + squareWidthPixel * ((nCheckerRow/2 - .5) - (iRow-1));
        
        % Get the upper left and lower right corners of each square
        iSquareLeft     = iSquareCenterX - squareWidthPixel/2;
        iSquareTop      = iSquareCenterY + squareWidthPixel/2;
        iSquareRight    = iSquareCenterX + squareWidthPixel/2;
        iSquareBottom   = iSquareCenterY - squareWidthPixel/2;
        
        stimulusSquaresArray(1:4, stimulusIndex) = [iSquareLeft; iSquareTop; iSquareRight; iSquareBottom];
    end
end




% ***********************************************************************
%                       START THE TRIAL LOOP
% ***********************************************************************
% Cue the user to begin the task by pressing the space bar
DrawFormattedText(window, 'Press the space to begin.', 'center', 'center', [0 200 0]);
Screen('Flip', window);
junk=NaN;
beforeStart = GetSecs;
while isnan(junk) && GetSecs - beforeStart < 30;
    [kd, sec, kc] = KbCheck;
    if (kd==1) && (kc(32)==1 || kc(44) == 1)
        junk = 1;
    elseif (kd==1) && (kc(81)==1)
        clear junk
    end
end
Screen('FillRect', window, backGround);
[vbl SOT] = Screen('Flip', window);



runningTask = 1;

% for iTrial = 0 : 2
while runningTask
    
    
        %   GENERATE THE CHECKERED STIMULUS:
        %   This algorithm ensures colors will be vertically evenly
        %   dispersed across the checkerboard (a vertical anti-cluster
        %   algorithm)
        % -------------------------------
        acceptTargetFlag = 0;
        randomThreshold = rand;
        while ~acceptTargetFlag
            % Randomize the colors with some proportion of blue squares.
            randomizedRightTarg         = int16(randperm(length(rightTargPropArray)));
            iPropRightTargIndex    = randomizedRightTarg(1);
            iRandomRightP    = rightTargPropArray(iPropRightTargIndex);
            
            if iRandomRightP == .5 && randomThreshold < fiftyPercentRate
                acceptTargetFlag = 1;
            elseif iRandomRightP < .5 && randomThreshold > rightTargetRate
                acceptTargetFlag = 1;
            elseif iRandomRightP > .5 && randomThreshold <= rightTargetRate
                acceptTargetFlag = 1;
            end
        end
        
        
        %   WHICH IS THE TARGET?
        % -------------------------------
        if iRandomRightP > .5
            rightTargetTrial = 1;
            majorityColor = 0;
            minorityColor = 1;
            iMinorityP = 1 - iRandomRightP;
        elseif iRandomRightP < .5
            rightTargetTrial = 0;
            majorityColor = 1;
            minorityColor = 0;
            iMinorityP = iRandomRightP;
        elseif iRandomRightP == .5
            rightTargetTrial = randi(2)-1;  % assign the target randomly
            majorityColor = 0;  % arbitrary
            minorityColor = 1;
            iMinorityP = iRandomRightP;
        end
        
        
        
        
        %         randomRightIndices       = randperm(nSquare);
        %         randomRight              = randomRightIndices(1 : nRight);
        
        
        iCheckerboardArray      = majorityColor + zeros(nSquare, 1);  % Initialize all checkers to majority (target) color
        nRowRemain = nCheckerRow;
        nMinority                   = int16(iMinorityP * nSquare);
        % Initialize a maximum and minimum allowed minority checkers for the first row
        nMax 	= ceil(nMinority / nRowRemain);
        nMin 	= floor(nMinority / nRowRemain);
        
        for iRow = 1 : nCheckerRow
            
            % Choose the number of minorityColor checkers to insert in the row
            randInsert = randi(2)-1;
            if randInsert == 0
                nInsert = nMax;
            elseif randInsert == 1
                nInsert = nMin;
            end
            % Fill row with as many zeros as there will be Targ1 squares
                iSquare 				= (iRow-1) * nCheckerColumn + 1;
                iCheckerboardArray(iSquare : iSquare + nInsert-1) = minorityColor;
            % randomly shuffle the ones and zeros within the row
            for iColumn = 1 : nCheckerColumn
                iSquare 				= (iRow-1) * nCheckerColumn + iColumn;
                tempIndex 				= (iRow-1) * nCheckerColumn + randi(nCheckerColumn);
                tempColor 				= iCheckerboardArray(tempIndex);
                iCheckerboardArray(tempIndex) 	= iCheckerboardArray(iSquare);
                iCheckerboardArray(iSquare) 		= tempColor;
            end
            % Update the varialbes used to calculate how many minority colors
            % checkers will go into the next row
            nRowRemain = nRowRemain - 1;
            nMinority = nMinority - nInsert;
            nMax 	= ceil(nMinority / nRowRemain);
            nMin 	= floor(nMinority / nRowRemain);
        end
        randomRight = iCheckerboardArray == 0;
        
        
        
        
        cyanGun = 174;
        magentaGun = 255;
        switch magentaTargSide
            case 'right'
                leftTargetCheckerColor = [0 cyanGun cyanGun];
                rightTargetCheckerColor = [magentaGun 0 magentaGun];
            case 'left'
                rightTargetCheckerColor = [0 cyanGun cyanGun];
                leftTargetCheckerColor = [magentaGun 0 magentaGun];
        end
        
        
        % Make all squares left target color as default
        iStimulusColorsArray(1, :) = leftTargetCheckerColor(1);
        iStimulusColorsArray(2, :) = leftTargetCheckerColor(2);
        iStimulusColorsArray(3, :) = leftTargetCheckerColor(3);
        % Then add the right target
        iStimulusColorsArray(1, randomRight) = rightTargetCheckerColor(1);
        iStimulusColorsArray(2, randomRight) = rightTargetCheckerColor(2);
        iStimulusColorsArray(3, randomRight) = rightTargetCheckerColor(3);
        
    
        
        
        
       
         Screen('FillRect', window, iStimulusColorsArray, stimulusSquaresArray);
       [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
       tic
       while toc < secPerStim
       end
   Screen('FillRect', window, backGround);
[vbl SOT] = Screen('Flip', window);
 
    
    % Allow option to quit task by pressing escape at the end of a trial
    tic
    while toc < 2
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        if keyIsDown && find(keyCode) == stopKey
            runningTask = 0;
        end
    end
    
    
end % trial loop

cleanup

function cleanup
Screen('CloseAll');
ShowCursor;
end


end


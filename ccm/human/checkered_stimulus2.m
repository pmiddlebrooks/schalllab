function checkered_stimulus2(SSDArrayScreenFlips, stimulusEccentricity, stimulusAngle, targetEccentricity, target1Angle, saveFlag)

% stimulusEccentricity: in some units (e.g. pixels), the distance from the
%       center of the screen to the center of the stimulus
% stimulusAngle: in degrees, the angle from the center of the screen to the
%       center of the stimulus

% Format (numbering) of the squares in the stimulus (a 3 X 3 example):
%             0   1   2
%             3   4   5
%             6   7   8

% example SSDArrayScreenFlips = [22 26 30 34]

initials    = input('Enter initials of subject            ', 's');
clockVector = clock;
session     = [num2str(year(now)), '_', num2str(month(now)), '_', num2str(day(now)), '_', num2str(clockVector(4)), '_', num2str(clockVector(5))];
saveFileName = ['choiceStop_', initials, '_', session];

whichScreen     = 0;
window          = Screen('OpenWindow', whichScreen, [0 0 0]);
commandwindow;
frameFrequency  = Screen('Nominalframerate', window);
    priorityLevel= MaxPriority(window);
    Priority(priorityLevel);

%     if frameFrequency == 0
%     frameFrequency = 60; %Hz
% else
%     disp('Current screen refresh is not 60Hz- need to change timing parameters in the script');
%     Screen('CloseAll');
%     return
% end;
    timing=Screen('GetFlipInterval',window);% get the flip rate of current monitor.
    a='Resolution & Refresh ok';
    if timing>.012 || timing<.011 %make sure we're running in 1024*768 at 85 Hz, else stop
        clear a;
        display('Please change screen to 1024x768 and 85Hz');
    end
    display a;



% ***************************************************************
%     Generate a new random number stream based on clock
% ***************************************************************
%        Replace the default stream with a stream whose seed is based on CLOCK, so
%        RAND will return different values in different MATLAB sessions.  NOTE: It
%        is usually not desirable to do this more than once per MATLAB
%        session.
s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setDefaultStream(s);


% ***************************************************************
%     CONSTANTS
% ***************************************************************
stopTrialProportion     = 1/3;%.3333;
screenRefreshRate       = 1 / frameFrequency;
SSDArray                = SSDArrayScreenFlips * (screenRefreshRate - .002); %subtracting 2 ms each cycle to enusre the SSD occurs before the next screen refresh
maxGoWait               = 1.5;
maxResponseWait         = 1.7;
stimulusScaleConversion = 10 / 100; % for now 10 pixels per every 100 pixels from fixation
feedbackTime            = 1.5;

% define the units you're using:
scrsz                   = get(0, 'ScreenSize');
screenWidth             = scrsz(3);
screenHeight            = scrsz(4);

% Target selection keys
rightTargetKey          = KbName('RightArrow');
leftTargetKey           = Kbname('LeftArrow');
escapeKey               = KbName('ESCAPE');


% Fixation Spot constants
fixationWidth       = 16;
fixationColor       = [200 200 200];
fixationHoldTime    = 1.5;

% Go signal constants
goSignalWidth       = fixationWidth - 2;
goSignalColor       = [50 50 50];

% Stop signal constants
stopSignalWidth     = fixationWidth - 2;
stopSignalColor     = [20 255 20];

% Target constants
rightTarget         = 1;
leftTarget          = 2;

targetWidth         = targetEccentricity * stimulusScaleConversion;
distractorWidth     = targetWidth;
% distractorAngle   = targetAngle + 180;
target2Angle        = target1Angle + 180;
targetColor         = [200 200 200];
distractorColor     = targetColor;

% Checkered Stimulus constants
stimulusColumns     = 10;
stimulusRows        = 10;
nSquares            = stimulusColumns * stimulusRows;
squareWidth         = 5;
% proportionBlues   = [0 1];
proportionBlues     = [.44 .48 .52 .56];



% ***********************************************************************
%       INITIALIZE STUFF
% ***********************************************************************
% Variables that will get filled each trial
trialOnset              = [];
trialDuration           = [];
fixationSpotOnset       = [];
fixationSpotDuration    = [];
targetOnset             = [];
targetDuration          = [];
distractorOnset         = [];
distractorDuration      = [];
stopOnset               = [];
stopDuration            = [];
choiceStimulusOnset     = [];
choiceStimulusDuration  = [];
responseCueOnset        = [];
responseOnset           = [];
feedbackOnset           = [];
feedbackDuration        = [];

choiceStimulusColor     = {};
choiceStimulusColorProportion = [];
targetPosition          = [];
distractorPosition      = [];

SSD                     = [];
realSSD                 = []; % calculated from Screen Flips
iSSDIndex               = 1;
nSSD                    = length(SSDArray);
lastStopOutcome         = 'stopIncorrectTarget';


% Initialize variables
stimulusSquaresArray    = zeros(4, nSquares);
iStimulusColorsArray    = zeros(3, nSquares);



% ***********************************************************************
%     STIMULI POSITION INFORMATION
% ***********************************************************************

% ------------------------------
%     FIXATION STIMULUS
% ------------------------------
fixationOffsetX     = 0;
fixationOffsetY     = 0;
fixationLeft        = screenWidth/2 - fixationWidth/2;
fixationTop         = screenHeight/2 - fixationWidth/2;
fixationRight       = screenWidth/2 + fixationWidth/2;
fixationBottom      = screenHeight/2 + fixationWidth/2;
fixationLocation    = [fixationLeft fixationTop fixationRight fixationBottom];

% ------------------------------
%     GO SIGNAL STIMULUS
% ------------------------------
goSignalLeft        = screenWidth/2 - goSignalWidth/2;
goSignalTop         = screenHeight/2 - goSignalWidth/2;
goSignalRight       = screenWidth/2 + goSignalWidth/2;
goSignalBottom      = screenHeight/2 + goSignalWidth/2;
goSignalLocation    = [goSignalLeft goSignalTop goSignalRight goSignalBottom];

% ------------------------------
%     STOP SIGNAL STIMULUS
% ------------------------------
stopSignalLeft      = screenWidth/2 - stopSignalWidth/2;
stopSignalTop       = screenHeight/2 - stopSignalWidth/2;
stopSignalRight     = screenWidth/2 + stopSignalWidth/2;
stopSignalBottom    = screenHeight/2 + stopSignalWidth/2;
stopSignalLocation  = [stopSignalLeft stopSignalTop stopSignalRight stopSignalBottom];

% ------------------------------
%     TARGET STIMULUI
% ------------------------------
target1OffsetX      = screenWidth/2 + targetEccentricity * cosd(target1Angle);
target1OffsetY      = screenHeight/2 - targetEccentricity * sind(target1Angle);

target1Left         = target1OffsetX - targetWidth/2;
target1Top          = target1OffsetY - targetWidth/2;
target1Right        = target1OffsetX + targetWidth/2;
target1Bottom       =   target1OffsetY + targetWidth/2;
target1Location     = [target1Left, target1Top, target1Right, target1Bottom];

target2OffsetX      = screenWidth/2 + targetEccentricity * cosd(target2Angle);
target2OffsetY      = screenHeight/2 - targetEccentricity * sind(target2Angle);

target2Left         = target2OffsetX - targetWidth/2;
target2Top          = target2OffsetY - targetWidth/2;
target2Right        = target2OffsetX + targetWidth/2;
target2Bottom       = target2OffsetY + targetWidth/2;
target2Location     = [target2Left, target2Top, target2Right, target2Bottom];

% ------------------------------
%     CHECKERED STIMULUS
% ------------------------------
% First determine the center of the 10X10 stimulus in x,y coordinates w.r.t the center of the
% screen
stimulusOffsetX     = screenWidth/2 + stimulusEccentricity * cosd(stimulusAngle);
stimulusOffsetY     = screenHeight/2 - stimulusEccentricity * sind(stimulusAngle);
% Get the positions of each checker square
for y = 1 : stimulusRows
    for x = 1 : stimulusColumns
        stimulusIndex   = x + (stimulusColumns * (y - 1));
        % Get the centers of each square
        iSquareCenterX  = stimulusOffsetX + squareWidth * ((x-1) - (stimulusColumns/2 - .5));
        iSquareCenterY  = stimulusOffsetY + squareWidth * ((stimulusRows/2 - .5) - (y-1));
        
        % Get the upper left and lower right corners of each square
        iSquareLeft     = iSquareCenterX - squareWidth/2;
        iSquareTop      = iSquareCenterY + squareWidth/2;
        iSquareRight    = iSquareCenterX + squareWidth/2;
        iSquareBottom   = iSquareCenterY - squareWidth/2;
        
        stimulusSquaresArray(1:4, stimulusIndex) = [iSquareLeft; iSquareTop; iSquareRight; iSquareBottom];
    end
end




% ***********************************************************************
%     START THE TRIAL LOOP
% ***********************************************************************
runningTask = 1;
iTrial = 0;
taskStartTime = GetSecs;
while runningTask
    iTrialOnset             = GetSecs - taskStartTime;
    iTrialOnsetComputerTime = GetSecs;
    iTrial                  = iTrial + 1;
    % Initialize variables that may or may not get filled
    iStopOnset              = nan;
    iStopDuration           = nan;
    iResponseOnset          = nan;
    iRealSSD                = nan;
    iSSD                    = nan;
    
    
    
    %   GENERATE THE CHECKERED STIMULUS
    % -------------------------------
    % Randomize the colors with some proportion of blue squares.
    randomizedBlues         = int16(randperm(length(proportionBlues)));
    proportionBluesIndex    = randomizedBlues(1);
    randomBlueProportion    = proportionBlues(proportionBluesIndex);
    nBlue                   = int16(randomBlueProportion * nSquares);
    randomBlueIndices       = randperm(nSquares);
    randomBlue              = randomBlueIndices(1 : nBlue);
    
    
    % *** NEED TO MAKE COLORS ISOLUMINANT  ***
    % Make all squares red as default
    iStimulusColorsArray(1, :) = 255;
    iStimulusColorsArray(2, :) = 0;
    iStimulusColorsArray(3, :) = 255;
    % Then add the blue square indices
    iStimulusColorsArray(1, randomBlue) = 0;
    iStimulusColorsArray(2, randomBlue) = 255;
    iStimulusColorsArray(3, randomBlue) = 255;
    
    
    %   WHICH IS THE TARGET?
    % -------------------------------
    if randomBlueProportion > .5
        iTrialTarget        = rightTarget;
        iTargetPosition     = [target1OffsetX - screenWidth/2, screenHeight/2 - target1OffsetY];
        iDistractorPosition = [target2OffsetX - screenWidth/2, screenHeight/2 - target2OffsetY];
    else
        iTrialTarget        = leftTarget;
        iTargetPosition     = [target2OffsetX - screenWidth/2, screenHeight/2 - target2OffsetY];
        iDistractorPosition = [target1OffsetX - screenWidth/2, screenHeight/2 - target1OffsetY];
    end
    
    %     %   WHICH IS THE TARGET?----- FOR RED AND BLUE TARGETS (MATCH TO
    %     SAMPLE VERSION)
    %     % -------------------------------
    %     if randomBlueProportion > .5 && blueTargetSide == 1
    %         iTrialTarget = rightTarget;
    %     elseif randomBlueProportion > .5 && blueTargetSide == 2
    %         iTrialTarget = leftTarget;
    %     elseif randomBlueProportion < .5 && blueTargetSide == 1
    %         iTrialTarget = leftTarget;
    %     elseif randomBlueProportion < .5 && blueTargetSide == 2
    %         iTrialTarget = rightTarget;
    %     end
    
    
    %    GO OR STOP TRIAL?
    % -------------------------------
    randomProportion = rand;
    if randomProportion > stopTrialProportion
        iTrialType = 'go';
    else
        iTrialType = 'stop';
    end
    
    
    
    % **********************************
    %     FIXATION SPOT
    % **********************************
    Screen('FillRect', window, fixationColor, fixationLocation);
    [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
    iFixationSpotOnset = StimulusOnsetTime - iTrialOnsetComputerTime;
    
    
    % **********************************
    %  FIXATION + TARGETS + CHECKERED STIMULUS + GO
    % **********************************
    Screen('FillRect', window, targetColor, target1Location);
    Screen('FillRect', window, targetColor, target2Location);
    Screen('FillRect', window, iStimulusColorsArray, stimulusSquaresArray);
    Screen('FillRect', window, fixationColor, fixationLocation);
    Screen('FillRect', window, goSignalColor, goSignalLocation);
    
    % some script, like:
    % "when eye in window, hold time, then move on
    waitSecs(fixationHoldTime);
    [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
    iFixationSpotDuration   = StimulusOnsetTime - iTrialOnsetComputerTime - iFixationSpotOnset;
    iTargetOnset            = StimulusOnsetTime - iTrialOnsetComputerTime;
    iDistractorOnset        = StimulusOnsetTime - iTrialOnsetComputerTime;
    iChoiceStimulusOnset    = StimulusOnsetTime - iTrialOnsetComputerTime;
    iResponseCueOnset       = StimulusOnsetTime - iTrialOnsetComputerTime;
    
    switch iTrialType
        case 'go'
            responded = 0;
            tic
            while toc < maxGoWait && ~responded
                [ keyIsDown, seconds, keyCode ] = KbCheck;
                if keyIsDown
                    iResponseOnset = iResponseCueOnset + toc;
                    if iTrialTarget == rightTarget && find(keyCode) == rightTargetKey
                        % Flag it was correct right
                        iTrialOutcome = 'goCorrectTarget';
                    elseif iTrialTarget == rightTarget && find(keyCode) == leftTargetKey
                        % Flag as distractor left
                        iTrialOutcome = 'goCorrectDistractor';
                    elseif iTrialTarget == leftTarget && find(keyCode) == leftTargetKey
                        % Flag as correct left
                        iTrialOutcome = 'goCorrectTarget';
                    elseif iTrialTarget == leftTarget && find(keyCode) == rightTargetKey
                        % Flag as distractor right
                        iTrialOutcome = 'goCorrectDistractor';
                    elseif find(keyCode) ~= leftTargetKey && find(keyCode) ~= rightTargetKey
                        % Flag as error key response
                        iTrialOutcome = 'errorKeyPress';
                    end
                    responded = 1;
                    WaitSecs(.5);
                end
            end
            if ~keyIsDown
                iTrialOutcome = 'goIncorrect';
            end
            
        case 'stop'
            iSSDIndex   = staircase(lastStopOutcome, iSSDIndex, nSSD);
            iSSD        = SSDArray(iSSDIndex);
            %             randomizedSSDs = int16(randperm(length(SSDArray)));
            %             iSSD = SSDArray(randomizedSSDs(1));
            
            % **********************************
            %  FIXATION + TARGETS + CHECKERED STIMULUS + STOP
            % **********************************
            Screen('FillRect', window, targetColor, target1Location);
            Screen('FillRect', window, targetColor, target2Location);
            Screen('FillRect', window, iStimulusColorsArray, stimulusSquaresArray);
            Screen('FillRect', window, fixationColor, fixationLocation);
            Screen('FillRect', window, stopSignalColor, stopSignalLocation);
            
            % Wait for the SSD to turn on the stop signal
            responded = 0;
            tic
            while toc < iSSD && ~responded
                [ keyIsDown, seconds, keyCode ] = KbCheck;
                if keyIsDown
                    % End the current trial and record the response made
                    iResponseOnset = iResponseCueOnset + toc;
                    if iTrialTarget == rightTarget && find(keyCode) == rightTargetKey
                        % Flag it would've been correct right but RT was before SSD
                        iTrialOutcome = 'stopIncorrectPreSSDTarget';
                    elseif iTrialTarget == rightTarget && find(keyCode) == leftTargetKey
                        % Flag it would've been distractor left but RT was before SSD
                        iTrialOutcome = 'stopIncorrectPreSSDDistractor';
                    elseif iTrialTarget == leftTarget && find(keyCode) == leftTargetKey
                        % Flag it would've been correct left but RT was before SSD
                        iTrialOutcome = 'stopIncorrectPreSSDTarget';
                    elseif iTrialTarget == leftTarget && find(keyCode) == rightTargetKey
                        % Flag it would've been distractor right but RT was before SSD
                        iTrialOutcome = 'stopIncorrectPreSSDDistractor';
                    elseif find(keyCode) ~= leftTargetKey && find(keyCode) ~= rightTargetKey
                        % Flag as error key response
                        iTrialOutcome = 'errorKeyPress';
                    end
                    responded = 1;
                    WaitSecs(.5);
                end
            end
            [~, stopOnsetTime, ~, ~, ~] = Screen('Flip', window);
            iRealSSD    = stopOnsetTime - iTargetOnset;
            iStopOnset  = stopOnsetTime - iTrialOnsetComputerTime;
            
            while toc < maxResponseWait && ~responded
                [ keyIsDown, seconds, keyCode ] = KbCheck;
                if keyIsDown
                    % End the current trial and record the response made
                    iResponseOnset = iResponseCueOnset + toc;
                    if iTrialTarget == rightTarget && find(keyCode) == rightTargetKey
                        % Flag it incorrected stop to correct right
                        iTrialOutcome = 'stopIncorrectTarget';
                    elseif iTrialTarget == rightTarget && find(keyCode) == leftTargetKey
                        % Flag it incorrected stop to distractor left
                        iTrialOutcome = 'stopIncorrectDistractor';
                    elseif iTrialTarget == leftTarget && find(keyCode) == leftTargetKey
                        % Flag it incorrected stop to correct leftTarget
                        iTrialOutcome = 'stopIncorrectTarget';
                    elseif iTrialTarget == leftTarget && find(keyCode) == rightTargetKey
                        % Flag it incorrected stop to distractor right
                        iTrialOutcome = 'stopIncorrectDistractor';
                    elseif find(keyCode) ~= leftTargetKey && find(keyCode) ~= rightTargetKey
                        % Flag as error key response
                        iTrialOutcome = 'errorKeyPress';
                    end
                    responded = 1;
                    WaitSecs(.5);
                end
            end
            if ~keyIsDown
                % Flag it as a correctful stop
                iTrialOutcome = 'stopCorrect';
            end
            lastStopOutcome = iTrialOutcome;
        otherwise
            disp('something wrong- neither a sotp or a go trial?')
    end
    
    % Allow option to quit task by pressing escape at the end of a trial
    tic
    while toc < 1
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        if keyIsDown && find(keyCode) == escapeKey
            runningTask = 0;
        end
    end
    [window, timeout] = feedback(window, iTrialOutcome);
    [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
    iTargetDuration         = StimulusOnsetTime - iTrialOnsetComputerTime + iTargetOnset;
    iDistractorDuration     = StimulusOnsetTime - iTrialOnsetComputerTime + iDistractorOnset;
    iChoiceStimulusDuration = StimulusOnsetTime - iTrialOnsetComputerTime + iChoiceStimulusOnset;
    iFeedbackOnset          = StimulusOnsetTime - iTrialOnsetComputerTime;
    iFeedbackDuration       = feedbackTime + timeout;
    if strcmp(iTrialType, 'stop')
        iStopDuration = StimulusOnsetTime - iTrialOnsetComputerTime + iStopOnset;
    end
    
    waitSecs(feedbackTime + timeout);
    iTrialDuration = GetSecs - iTrialOnsetComputerTime;
    
    
    
    % Add the trial's variables to the data set.
    
    % Event Timing
    % ------------
    trialOnset              = [trialOnset; iTrialOnset];
    trialDuration           = [trialDuration; iTrialDuration];
    fixationSpotOnset       = [fixationSpotOnset; iFixationSpotOnset];
    fixationSpotDuration    = [fixationSpotDuration; iFixationSpotDuration];
    targetOnset             = [targetOnset; iTargetOnset];
    targetDuration          = [targetDuration; iTargetDuration];
    distractorOnset         = [distractorOnset; iDistractorOnset];
    distractorDuration      = [distractorDuration; iDistractorDuration];
    choiceStimulusOnset     = [choiceStimulusOnset; iChoiceStimulusOnset];
    choiceStimulusDuration  = [choiceStimulusDuration; iChoiceStimulusDuration];
    stopOnset               = [stopOnset; iStopOnset];
    stopDuration            = [stopDuration; iStopDuration];
    responseCueOnset        = [responseCueOnset; iResponseCueOnset];
    responseOnset           = [responseOnset; iResponseOnset];
    feedbackOnset           = [feedbackOnset; iFeedbackOnset];
    feedbackDuration        = [feedbackDuration; iFeedbackDuration];
    
    
    % Stimulus Properties
    % ------------
    choiceStimulusColor     = [choiceStimulusColor; iStimulusColorsArray];
    choiceStimulusColorProportion = [choiceStimulusColorProportion; randomBlueProportion];
    targetPosition          = [targetPosition; iTargetPosition];
    distractorPosition      = [distractorPosition; iDistractorPosition];
    
    %     SSD = [SSD; iSSD];
    %     realSSD = [realSSD; iRealSSD];
    
end


nTrial = length(trialOnset);




% ********************************************************************
%                    Trial Data
% ********************************************************************


% Event Timing
% ------------




% Event Properties
% ---------------------------------------------------------------
stopTrialXXXX             = ones(nTrial, 1) .* stopTrialProportion;



% Location of Stimuli
% ---------------------------------------------------------------

fixationPosition          = [ones(nTrial, 1) * fixationOffsetX, ones(nTrial, 1) * fixationOffsetY];
choiceStimulusPosition    = [ones(nTrial, 1) * (stimulusOffsetX - screenWidth/2), ones(nTrial, 1) * (screenHeight/2 - stimulusOffsetY)];

fixationSize              = ones(nTrial, 1) * fixationWidth;
targetSize                = ones(nTrial, 1) * targetWidth;
distractorSize            = ones(nTrial, 1) * distractorWidth;
choiceStimulusSize        = ones(nTrial, 1) * squareWidth * stimulusColumns;
fixationColor             = ones(nTrial, 1) * fixationColor;
targetColor               = ones(nTrial, 1) * targetColor;
distractorColor           = ones(nTrial, 1) * distractorColor;



trialData = dataset(...
    {trialOnset,            'trialOnset'},...
    {trialDuration,         'trialDuration'},...
    {fixationSpotOnset,     'fixationSpotOnset'},...
    {fixationSpotDuration,  'fixationSpotDuration'},...
    {targetOnset,           'targetOnset'},...
    {targetDuration,        'targetDuration'},...
    {distractorOnset,       'distractorOnset'},...
    {distractorDuration,    'distractorDuration'},...
    {choiceStimulusOnset,   'choiceStimulusOnset'},...
    {choiceStimulusDuration, 'choiceStimulusDuration'},...
    {stopOnset,             'stopOnset'},...
    {stopDuration,          'stopDuration'},...
    {responseCueOnset,      'responseCueOnset'},...
    {responseOnset,         'responseOnset'},...
    {feedbackOnset,         'feedbackOnset'},...
    {feedbackDuration,      'feedbackDuration'},...
    {stopTrialXXXX,         'stopTrialXXXX'},...
    {choiceStimulusColorProportion, 'choiceStimulusColorProportion'},...
    {fixationPosition,      'fixationPosition'},...
    {targetPosition,        'targetPosition'},...
    {distractorPosition,    'distractorPosition'},...
    {choiceStimulusPosition, 'choiceStimulusPosition'},...
    {fixationSize,          'fixationSize'},...
    {targetSize,            'targetSize'},...
    {distractorSize,        'distractorSize'},...
    {choiceStimulusSize,    'choiceStimulusSize'},...
    {fixationColor,         'fixationColor'},...
    {targetColor,           'targetColor'},...
    {distractorColor,       'distractorColor'},...
    {choiceStimulusColor,   'choiceStimulusColor'});
    

% ********************************************************************
% Session Data
% ********************************************************************


sessionData.task.effector = 'keypress';

sessionData.timing.totalDuration = trialOnset(end) + trialDuration(end) - trialOnset(1); % seconds


sessionData.subjectID = initials;
sessionData.sessionID = session;












% ---------- Window Cleanup ----------

% Closes all windows.
Screen('CloseAll');

% Restores the mouse cursor.
ShowCursor;


% % Fill the data set with the data
% trial = (1 : iTrial)';
% experimentDataSet = dataset({trial, 'trial'}, {trialOnset, 'trialOnset'}, {trialType, 'trialType'}, {trialOutcome, 'trialOutcome'}, ...
%     {reactionTime, 'reactionTime'}, {SSD, 'SSD'}, {realSSD, 'realSSD'}, {blueProportion, 'blueProportion'});

if saveFlag
    save(saveFileName, 'trialData', 'sessionData');
    %     save(saveFileName, 'experimentDataSet', 'experimentStartTime', 'stimulusEccentricity', 'stimulusAngle', 'targetEccentricity', ...
    %         'target1Angle', 'stopTrialProportion', 'SSDArray', 'maxGoWait', 'maxResponseWait', 'stimulusScaleConversion', ...
    %         'fixationWidth', 'fixationColor', 'fixationHoldTime', 'goSignalWidth', 'goSignalColor', 'stopSignalWidth', 'stopSignalColor', ...
    %         'targetWidth', 'targetColor', 'stimulusColumns', 'stimulusRows', 'squareWidth', 'proportionBlues', 'stimulusColorsArray');
end



% *******************************************************************
function [window, timeout] = feedback(window, iTrialOutcome)


scrsz = get(0, 'ScreenSize');
screenWidth = scrsz(3);
screenHeight = scrsz(4);

incorrectTextColor = [250, 50, 50];
correctTextColor = [50, 220, 50];
Screen('TextFont', window, 'Times');
Screen('TextSize', window, 70);
Screen('TextStyle', window, 1);
if strcmp(iTrialOutcome, 'goCorrectTarget') || strcmp(iTrialOutcome, 'stopCorrect')
    [nx, ny, bbox] = DrawFormattedText(window, 'Nice job', 'center', 'center', correctTextColor);
    timeout = 0;
elseif strcmp(iTrialOutcome, 'goCorrectDistractor')
    DrawFormattedText(window, 'Wrong Target', 'center', 'center', incorrectTextColor);
    timeout = .5;
elseif strcmp(iTrialOutcome, 'goIncorrect')
    DrawFormattedText(window, 'You should have responded', 'center', 'center', incorrectTextColor);
    timeout = 1;
elseif strcmp(iTrialOutcome, 'stopIncorrectTarget')
    DrawFormattedText(window, 'You Should Have Stopped', 'center', 'center', incorrectTextColor);
    timeout = .5;
elseif strcmp(iTrialOutcome, 'stopIncorrectDistractor')
    DrawFormattedText(window, 'You Should Have Stopped', 'center', 'center', incorrectTextColor);
    timeout = .5;
elseif strcmp(iTrialOutcome, 'errorKeyPress')
    DrawFormattedText(window, 'Careful which key you press', 'center', 'center', [190, 190, 190]);
    timeout = .5;
elseif strcmp(iTrialOutcome, 'stopIncorrectPreSSDTarget') || strcmp(iTrialOutcome, 'stopIncorrectPreSSDDistractor')
    DrawFormattedText(window, 'Wow, that was a quck response', 'center', 'center', [190, 190, 190]);
    timeout = .5;
end

% Screen('Flip', window);
% waitSecs(feedbackTime + timeout)
return


% *******************************************************************
function newSSDIndex = staircase(lastStopOutcome, lastSSDIndex, nSSD)

maxStepSize = 3;
iStepSize = randi(maxStepSize);

if strcmp(lastStopOutcome, 'stopCorrect')
    newSSDIndex = min(iStepSize + lastSSDIndex, nSSD);
elseif strcmp(lastStopOutcome, 'stopIncorrectTarget') || strcmp(lastStopOutcome, 'stopIncorrectDistractor') || ...
        strcmp(lastStopOutcome, 'stopIncorrectPreSSDTarget') || strcmp(lastStopOutcome, 'stopIncorrectPreSSDDistractor')
    newSSDIndex = max(iStepSize - lastSSDIndex, 1);
end
return

function choice_stop_task_saccade_practice



% stimulusAmplitude: in some units (e.g. pixels), the distance from the
%       center of the screen to the center of the stimulus
% stimulusAngle: in degrees, the angle from the center of the screen to the
%       center of the stimulus

% Format (numbering) of the squares in the stimulus (a 3 X 3 example):
%             0   1   2
%             3   4   5
%             6   7   8

% example SSDArrayScreenFlips = [22 26 30 34]

if nargin == 0
    SSDArrayScreenFlips = [6 : 6 : 6 + 6*10];
    choiceStimulusAmplitude = 3;
    choiceStimulusAngle = 90;
    targetAmplitude = 10;
    rightTargetAngle = 0;
    plotFlag = 1;
    saveFlag = 1;
end

% rightTargProportionArray   = [0 1];
% rightTargProportionArray     = [.41 .45 .48 .5 .52 .55 .59];
% rightTargProportionArray     = [.35 .42 .46 .5 .54 .58 .65];
rightTargProportionArray     = [.35 .45 .55 .65];
rightTargetProportion = .5; % How often should right side be target?
fiftyPercentRate = .6; % how often should 50% signal strength be presented RELATIVE TO OTHER proportions?

nProportion             = length(rightTargProportionArray);
% INITIAL_SSD_INDEX       = 6;
% iSSDIndexArray          = INITIAL_SSD_INDEX * ones(nProportion, 1);
INITIAL_SSD_INDEX       = [6 7 8 8 8 7 6];
iSSDIndexArray          = INITIAL_SSD_INDEX;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Set up Experiment Variables, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get distance values and conversions for stimuli presentation
metersFromScreen    = .54;
screenWidthMeters   = .4;
theta               = asind(screenWidthMeters / 2 / sqrt(((screenWidthMeters / 2)^2) + metersFromScreen^2));
screenSize          = get(0, 'ScreenSize');
screenWidthPixel  	= screenSize(3);
screenHeightPixel 	= screenSize(4);
matlabCenterX       = screenWidthPixel/2;
matlabCenterY       = screenHeightPixel/2;
pixelPerDegree     = screenWidthPixel / (2*theta)

targetAmplitudePixel = targetAmplitude * pixelPerDegree;
stimulusAmplitudePixel = choiceStimulusAmplitude * pixelPerDegree;


% Get info for saving a data file.
subjectID     = input('Enter subject number            ', 's');
clockVector     = clock;
session         = [datestr(now,'mmdd'),  'saccade'];
saveFileName    = ['data/', subjectID, session];


whichScreen     = 0;
backGround      = [50 50 50];
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


% ***************************************************************
%     CONSTANTS
% ***************************************************************
stopTrialProportion     = .43;
trialsPerBlock          = 10;
BLOCKS_TO_RUN           = 10;
totalTrial              = trialsPerBlock * BLOCKS_TO_RUN;
SSDArray                = SSDArrayScreenFlips * (flipTime); %subtracting 2 ms each cycle to enusre the SSD occurs before the next screen refresh
graceResponse           = 1.7; % seconds allowed to make a saccade
graceObtainFixation     = 2;  % seconds program will wait for user to obtain fixation
stopHoldDuration        = graceResponse * 1;
graceSaccadeDuration    = .1; % seconds allowed intra saccade time
postSaccadeHoldDuration = .4;  % duration to hold post-saccade fixation
stimulusScaleConversion = 1/ 10; % for now 10 pixels per every 100 pixels from fixation
feedbackTime            = .25;
dummymode               = 0;



% Keyboard assignments
KbName('UnifyKeyNames');
rightTargetKey          = KbName('m');
leftTargetKey           = Kbname('z');
escapeKey               = KbName('ESCAPE');


% Fixation Spot constants
fixationAngle       = 0;
fixationAmplitude   = 0;
fixationWidth       = .5;
fixationWidthPixel  = fixationWidth * pixelPerDegree;
fixWindowScale      = 3;  % fix window is 4 times size of fix point
fixationWindow      = fixationWidth*fixWindowScale;
fixationWindowPixel	= [-fixationWidthPixel*fixWindowScale, -fixationWidthPixel*fixWindowScale, fixationWidthPixel*fixWindowScale, fixationWidthPixel*fixWindowScale];
fixationWindowPixel      = CenterRect(fixationWindowPixel, centerPoint);
fixationColor       = [200 200 200];
fixationHoldBase    = .5;
fixationHoldAdd     = (0 : 10 : 500)  / 1000;


% Go signal constants
goSignalWidthPixel       = fixationWidthPixel - 2;
goSignalColor       = [50 50 50];

% Stop signal constants
stopSignalWidthPixel     = fixationWidthPixel - 2;
stopSignalColor     = [255 150 0];


% targetWidth         = .5 + (targetAmplitude * stimulusScaleConversion);
targetWidth         = 1;
targetWidthPixel   	= targetWidth * pixelPerDegree;
targetWindowScale   = 3;
targetWindow        = targetWindowScale * targetWidth;
targetWindowPixel 	= targetWindow * pixelPerDegree;
distractorWidth     = targetWidth*.5;
distractorWidthPixel     = distractorWidth * pixelPerDegree;
distractorWindow    = targetWindow;
distractorWindowPixel = targetWindowPixel;
% distractorAngle   = targetAngle + 180;
leftTargetAngle    	= rightTargetAngle + 180;
targetColor         = [200 200 200];
distractorColor     = targetColor;

% Checkered Stimulus constants
stimulusColumns     = 10;
stimulusRows        = 10;
nSquare            = stimulusColumns * stimulusRows;
squareWidthPixel         = 3; % pixels
squareWidth = squareWidthPixel / pixelPerDegree;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Set up Eyelink
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup Eyelink System- must have eyelink computer started and running
%eyelink
init = Eyelink('Initialize');
el = EyelinkInitDefaultsNC(window);
edffilename = [subjectID,  datestr(now,'mmdd'), 'em',  '.edf'];
openEDF = Eyelink('OpenFile', edffilename);

HideCursor;

EyelinkDoTrackerSetup(el, 'c');
%Drift Correction
drift = EyelinkDoDriftCorrect(el);
DrawFormattedText(window, 'Which eye? 0=Left; 1=Right.', 'center', 'center');
Screen('Flip', window);
junk = NaN;
while isnan(junk)
    [kde, sece, kce] = KbCheck;
    if (kde==1) && ((kce(48)==1) || (kce(49)==1) || (kce(96)==1) || (kce(97)==1))
        junk=1;
    elseif (kde==1) && (kce(81)==1)
        clear junk
    end
end
if ((kce(48)==1) || (kce(96)==1) ); eyeRecorded=0;
elseif ( (kce(49)==1) || (kce(97)==1));eyeRecorded=1;
else clear eyeRecorded;
end

Screen('FillRect', window, backGround);
[vbl, SOT] = Screen('Flip', window);

Eyelink('Command', 'sample_rate = 500');
if eyeRecorded == 0;
    Eyelink('Command', 'file_sample_filter=LEFT, GAZE, AREA, STATUS');
    Eyelink('Command', 'file_event_filter=LEFT,  FIXATION, SACCADE, BLINK, MESSAGE');
    Eyelink('Command', 'link_event_filter = LEFT, SACCADE,BLINK, MESSAGE');
else
    Eyelink('Command', 'file_sample_filter=RIGHT, GAZE, AREA, STATUS');
    Eyelink('Command', 'file_event_filter= RIGHT, FIXATION, SACCADE, BLINK, MESSAGE');
    Eyelink('Command', 'link_event_filter = RIGHT,SACCADE,BLINK, MESSAGE');
end

screenpixstring = sprintf('screen_pixel_coords= %f,%f,%f,%f', centerPoint);      %you don't want eyelink to take those from physical.ini, because that means changing physical.ini for every screen resolution setting
Eyelink('Command', 'clear_screen 0');









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             BEGIN TASK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% ***********************************************************************
%       INITIALIZE STUFF
% ***********************************************************************
% Variables that will get filled each trial
trialOutcome            = {};
trialOnset              = [];
trialDuration           = [];
abortOnset              = [];
fixationSpotOnset       = [];
fixationOnset           = [];
fixationSpotDuration    = [];
targetOnset             = [];
targetDuration          = [];
distractorOnset         = [];
distractorDuration      = [];
choiceStimulusOnset     = [];
choiceStimulusDuration  = [];
stopOnset               = [];
stopDuration            = [];
SSD                     = [];
realSSD                 = []; % calculated from Screen Flips
responseCueOnset        = [];
responseOnset           = [];
targetWindowEntered     = [];
feedbackOnset           = [];
feedbackDuration        = [];

targetAngle          = [];
distractorAngle      = [];

target1CheckerProportion = [];
choiceStimulusColor     = {};
checkerboardArray       = {};

nSSD                    = length(SSDArray);
lastStopOutcomeArray    = cell(nProportion, 1);

% Initialize variables
stimulusSquaresArray    = zeros(4, nSquare);
iStimulusColorsArray    = zeros(3, nSquare);


% ***********************************************************************
%     STIMULI POSITION INFORMATION
% ***********************************************************************

% ------------------------------
%     FIXATION STIMULUS
% ------------------------------
fixationEyeLinkX     = 0;
fixationEyeLinkY     = 0;
fixationLeft        = matlabCenterX - fixationWidthPixel/2;
fixationTop         = matlabCenterY - fixationWidthPixel/2;
fixationRight       = matlabCenterX + fixationWidthPixel/2;
fixationBottom      = matlabCenterY + fixationWidthPixel/2;
fixationSquare    = [fixationLeft fixationTop fixationRight fixationBottom];

% ------------------------------
%     GO SIGNAL STIMULUS
% ------------------------------
goSignalLeft        = matlabCenterX - goSignalWidthPixel/2;
goSignalTop         = matlabCenterY - goSignalWidthPixel/2;
goSignalRight       = matlabCenterX + goSignalWidthPixel/2;
goSignalBottom      = matlabCenterY + goSignalWidthPixel/2;
goSignalSquare    = [goSignalLeft goSignalTop goSignalRight goSignalBottom];

% ------------------------------
%     STOP SIGNAL STIMULUS
% ------------------------------
stopSignalLeft      = matlabCenterX - stopSignalWidthPixel/2;
stopSignalTop       = matlabCenterY - stopSignalWidthPixel/2;
stopSignalRight     = matlabCenterX + stopSignalWidthPixel/2;
stopSignalBottom    = matlabCenterY + stopSignalWidthPixel/2;
stopSignalLocation  = [stopSignalLeft stopSignalTop stopSignalRight stopSignalBottom];

% ------------------------------
%     TARGET STIMULUI
% ------------------------------
rightTargetEyeLinkX      = targetAmplitudePixel * cosd(rightTargetAngle);
rightTargetEyeLinkY      = targetAmplitudePixel * sind(rightTargetAngle);
rightTargetMatlabX      = matlabCenterX + rightTargetEyeLinkX;
rightTargetMatlabY      = matlabCenterY - rightTargetEyeLinkY;

rightTargetLeft         = rightTargetMatlabX - targetWidthPixel/2;
rightTargetTop          = rightTargetMatlabY - targetWidthPixel/2;
rightTargetRight        = rightTargetMatlabX + targetWidthPixel/2;
rightTargetBottom       = rightTargetMatlabY + targetWidthPixel/2;
rightTargetSquare     = [rightTargetLeft, rightTargetTop, rightTargetRight, rightTargetBottom];

leftTargetEyeLinkX      = targetAmplitudePixel * cosd(leftTargetAngle);
leftTargetEyeLinkY      = targetAmplitudePixel * sind(leftTargetAngle);
leftTargetMatlabX      = matlabCenterX + leftTargetEyeLinkX;
leftTargetMatlabY      = matlabCenterY - leftTargetEyeLinkY;

leftTargetLeft         = leftTargetMatlabX - targetWidthPixel/2;
leftTargetTop          = leftTargetMatlabY - targetWidthPixel/2;
leftTargetRight        = leftTargetMatlabX + targetWidthPixel/2;
leftTargetBottom       = leftTargetMatlabY + targetWidthPixel/2;
leftTargetSquare     = [leftTargetLeft, leftTargetTop, leftTargetRight, leftTargetBottom];

% ------------------------------
%     CHECKERED STIMULUS
% ------------------------------
% First determine the center of the 10X10 stimulus in x,y coordinates w.r.t the center of the
% screen
stimulusEyeLinkX     = stimulusAmplitudePixel * cosd(choiceStimulusAngle);
stimulusEyeLinkY     = stimulusAmplitudePixel * sind(choiceStimulusAngle);
stimulusMatlabX     = matlabCenterX + stimulusEyeLinkX;
stimulusMatlabY     = matlabCenterY - stimulusEyeLinkY;
stimulusWindowPixel      = squareWidthPixel * stimulusColumns; % For now, make it window the same size as the stimulus
choiceStimulusWindow = stimulusWindowPixel / pixelPerDegree;
stimulusWindow       = [stimulusEyeLinkY - stimulusWindowPixel/2, stimulusEyeLinkX - stimulusWindowPixel/2, stimulusEyeLinkY + stimulusWindowPixel/2, stimulusEyeLinkX + stimulusWindowPixel/2];
% Get the positions of each checker square
for iRow = 1 : stimulusRows
    for iColumn = 1 : stimulusColumns
        stimulusIndex   = iColumn + (stimulusColumns * (iRow - 1));
        % Get the centers of each square
        iSquareCenterX  = stimulusMatlabX + squareWidthPixel * ((iColumn-1) - (stimulusColumns/2 - .5));
        iSquareCenterY  = stimulusMatlabY + squareWidthPixel * ((stimulusRows/2 - .5) - (iRow-1));
        
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
DrawFormattedText(window, 'Press the space to begin block 1.', 'center', 'center', [0 200 0]);
Screen('Flip', window);
junk=NaN;
beforeStart = GetSecs;
while isnan(junk) && GetSecs - beforeStart < 30;
    [kd, sec, kc] = KbCheck;
    if (kd==1) && (kc(32)==1)
        junk = 1;
    elseif (kd==1) && (kc(81)==1)
        clear junk
    end
end
Screen('FillRect', window, backGround);
[vbl SOT] = Screen('Flip', window);



runningTask = 1;
iTrial = 1;
taskStartTime = GetSecs;
newTrialVariables = 1; % A flag, set to zero if an abort occurs and we want all trial variables to stay the same
iBlock = 1;

% for iTrial = 0 : 2
while runningTask
    
    
    %Take a break every "trialsPerBlock" trials
    if (mod(iTrial, trialsPerBlock) == 1) && (iTrial ~= 1)
        iBlock = iBlock + 1;
msg = sprintf('Press the space bar to begin block %d', iBlock);
        DrawFormattedText(window, msg, 'center', 'center', [0 200 0]);
        Screen('Flip', window);
        junky = NaN;
        WaitSecs(2);
        takeBreak = GetSecs;
        while isnan(junky) && GetSecs - takeBreak < 60;
            [kd, sec, kc] = KbCheck;
            if (kd==1) && (kc(32)==1)
                junky = 1;
            elseif (kd==1) && (kc(81)==1)
                clear junky
            end
%         if kd == 1 && find(kc) == escapeKey
%             runningTask = 0;
%         end
        end
    end
    
    
    
    
    
    
    
    
    iTrialOutcome          	= nan;
    % Initialize variables that may or may not get filled
    iTrialDuration       	= nan;
    iAbortOnset             = nan;
    fixationAddIndex        = randperm(length(fixationHoldAdd));
    iPreTargetFixDuration   = .4 + fixationHoldAdd(fixationAddIndex(1));
    iPostTargetFixDuration  = fixationHoldBase + fixationHoldAdd(fixationAddIndex(2));
    iFixationSpotOnset      = nan;
    iFixationOnset          = nan;
    iFixationSpotDuration   = nan;
    iTargetOnset            = nan;
    iTargetDuration         = nan;
    iDistractorOnset        = nan;
    iDistractorDuration     = nan;
    iChoiceStimulusOnset    = nan;
    iChoiceStimulusDuration = nan;
    iStopOnset              = nan;
    iStopDuration           = nan;
    iSSD                    = nan;
    iRealSSD                = nan;
    iResponseCueOnset       = nan;
    iResponseOnset          = nan;
    iTargetWindowEntered    = nan;
    iFeedbackOnset          = nan;
    iFeedbackDuration       = nan;
    
    
    
    
    
    
    
    
    
    % If the previous trial was aborted before it began (user did not
    % obtain fixation), then keep all the variables the same. Only update
    % them if it's a new trial
    if newTrialVariables
        
        %   GENERATE THE CHECKERED STIMULUS
        % -------------------------------
        acceptTargetFlag = 0;
        randomThreshold = rand;
        while ~acceptTargetFlag
            % Randomize the colors with some proportion of blue squares.
            randomizedRightTarg         = int16(randperm(length(rightTargProportionArray)));
            iProportionRightTargIndex    = randomizedRightTarg(1);
            randomRightProportion    = rightTargProportionArray(iProportionRightTargIndex);
            
            if randomRightProportion == .5 && randomThreshold < fiftyPercentRate
                acceptTargetFlag = 1;
            elseif randomRightProportion < .5 && randomThreshold > rightTargetProportion
                acceptTargetFlag = 1;
            elseif randomRightProportion > .5 && randomThreshold <= rightTargetProportion
                acceptTargetFlag = 1;
            end
        end
        iProportionRightTarg = randomRightProportion;
        nRight                   = int16(randomRightProportion * nSquare);
        randomRightIndices       = randperm(nSquare);
        randomRight              = randomRightIndices(1 : nRight);
        
        
        
        iCheckerboardArray      = ones(nSquare, 1);  % Initiate all checkers to left target (ones)
        cyanGun = 175;
        magentaGun = 255;
        leftTargetCheckerColor = [0 cyanGun cyanGun];
        rightTargetCheckerColor = [magentaGun 0 magentaGun];
        
        % *** NEED TO MAKE COLORS ISOLUMINANT  ***
        % Make all squares left target color as default
        iStimulusColorsArray(1, :) = leftTargetCheckerColor(1);
        iStimulusColorsArray(2, :) = leftTargetCheckerColor(2);
        iStimulusColorsArray(3, :) = leftTargetCheckerColor(3);
        % Then add the right target
        iStimulusColorsArray(1, randomRight) = rightTargetCheckerColor(1);
        iStimulusColorsArray(2, randomRight) = rightTargetCheckerColor(2);
        iStimulusColorsArray(3, randomRight) = rightTargetCheckerColor(3);
        iCheckerboardArray(randomRight)      = 0;  % Change right target checkers to zeros
        
        %   WHICH IS THE TARGET?
        % -------------------------------
        if randomRightProportion > .5
            rightTargetTrial = 1;
        elseif randomRightProportion < .5
            rightTargetTrial = 0;
        elseif randomRightProportion == .5
            rightTargetTrial = randi(2)-1;
        end
        if rightTargetTrial
            iTargetAngle        = rightTargetAngle;
            iDistractorAngle   	= leftTargetAngle;
            iTargetSquare       = rightTargetSquare;
            iDistractorSquare   = leftTargetSquare;
            iTargetLocation     = [rightTargetMatlabX - matlabCenterX, matlabCenterY - rightTargetMatlabY];
            iDistractorLocation = [leftTargetMatlabX - matlabCenterX, matlabCenterY - leftTargetMatlabY];
            iTargetWindow       = [rightTargetMatlabX - targetWindowPixel, rightTargetMatlabY - targetWindowPixel, rightTargetMatlabX + targetWindowPixel, rightTargetMatlabY + targetWindowPixel];
            iDistractorWindow   = [leftTargetMatlabX - distractorWindowPixel, leftTargetMatlabY - distractorWindowPixel, leftTargetMatlabX + distractorWindowPixel, leftTargetMatlabY + distractorWindowPixel];
        else
            iTargetAngle        = leftTargetAngle;
            iDistractorAngle   	= rightTargetAngle;
            iTargetSquare       = leftTargetSquare;
            iDistractorSquare   = rightTargetSquare;
            iTargetLocation     = [leftTargetMatlabX - matlabCenterX, matlabCenterY - leftTargetMatlabY];
            iDistractorLocation = [rightTargetMatlabX - matlabCenterX, matlabCenterY - rightTargetMatlabY];
            iTargetWindow       = [leftTargetMatlabX - targetWindowPixel, leftTargetMatlabY - targetWindowPixel, leftTargetMatlabX + targetWindowPixel, leftTargetMatlabY + targetWindowPixel];
            iDistractorWindow   = [rightTargetMatlabX - distractorWindowPixel, rightTargetMatlabY - distractorWindowPixel, rightTargetMatlabX + distractorWindowPixel, rightTargetMatlabY + distractorWindowPixel];
        end
    end
    
    
    %    GO OR STOP TRIAL?
    % -------------------------------
    randomProportion = rand;
    if randomProportion > stopTrialProportion
        iTrialType = 'go';
    else
        iTrialType = 'stop';
    end
    
    
    
    % ****************************************************************************************
    %            BEGIN STAGES
    % ****************************************************************************************
    
    
    % Initialize stage logical variables each trial
    addTrialDataFlag    =  1;  % Gets set to false for pre-fixation aborts
    stagePreFixation    = 1;
    stageFixation       = 0;
    stageTargetOn       = 0;
    stageChoiceOn       = 0;
    stageStopSignalOn   = 0;
    stageInFlight       = 0;
    stageOnTarget       = 0;
    stageOnDistractor   = 0;
    stageFeedback       = 0;
    
    preStopSignalResponse = 0;
    
    % ****************************************************************************************
    %            PRE-FIXATION STAGE
    % ****************************************************************************************
    % Turn on the fixation spot and wait forunit subject to start trial by
    % fixating and pressing space bar (which initiates a drift correction
    % and subsequently the trial
    
    
    Screen('FillRect', window, backGround);
    Screen('FillRect', window, fixationColor, fixationSquare);
    %     Screen('FrameRect', window, targetColor, fixationWindowPixel);
    %         DrawFormattedText(window, 'Press the space to drift correct.', 'center', 'center', [0 200 0]);
    [~, fixationOnsetTime, ~, ~, ~] = Screen('Flip', window);
    
    
    
    
    %Start Recording Eye for current trial
    status = ['record_status_message "Trial ' num2str(iTrial) '"'];
    Eyelink('Command', status);
    Eyelink('Message', '%s%d', 'Trial=', iTrial);
    
    
    
    drift = 1;
    tic
    while drift
        drift = EyelinkDoDriftCorrect(el, matlabCenterX, matlabCenterY, 0, 1);
        % drift = Eyelink('Command', 'do_drift_correct', 0, 0, 0, 0)
        if toc > 5
            drift = 0;
        end
    end
    
    iTrialOnset             = GetSecs - taskStartTime;
    iTrialOnsetComputerTime = GetSecs;
    Eyelink('StartRecording');
    
    
    
    
    
    
    iFixationSpotOnset = fixationOnsetTime - iTrialOnsetComputerTime;
    tic
    while stagePreFixation
        
        % Unit code for checking eye position
        %-----------------------------------------------------------------------------------------------------
        if dummymode == 0
            error = Eyelink('CheckRecording');
            if error
                fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
%                     runningTask = 0;
%                     cleanup;
iTrialOutcome = 'eyelinkError';
stageFeedback = 1;
break
            end
            
            if Eyelink( 'NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                event = Eyelink( 'NewestFloatSample');
                event.gx;
                event.gy;
                if eyeRecorded ~= -1 % do we know which eye to use yet?
                    % if we do, get current gaze position from sample
                    x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                    y = event.gy(eyeRecorded+1);
                    % do we have valid data and is the pupil visible?
                    if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                        eyeX = x;
                        eyeY = y;
                        
                        %                         % Show where gaze is (troubleshooting)
                        %                         gazeRect = [ x-5 y-5 x+6 y+6];
                        %                         colour = round(rand(3,1)*255); % coloured dot
                        %                         Screen('FillOval', window,colour, gazeRect);
                        %                         Screen('FillRect', window, fixationColor, fixationSquare);
                        %                         Screen('FrameRect', window, targetColor, fixationWindow);
                        %                         [~, fixationOnsetTime, ~, ~, ~] = Screen('Flip', window);
                        
                    end
                end
            end
        else
            
            % Query current mouse cursor position (our "pseudo-eyetracker") -
            % (mx,my) is our gaze position.
            [eyeX, eyeY, ~] = GetMouse(window); %#ok<*NASGU>
        end
        %-----------------------------------------------------------------------------------------------------
        
        
        if inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc <= graceObtainFixation
            Eyelink('Message', 'Fixation_Start');
            iFixationOnset      = GetSecs - iTrialOnsetComputerTime;
            stagePreFixation    = 0;
            stageFixation       = 1;
        elseif GetSecs > fixationOnsetTime + graceObtainFixation
            % If subject aborted, exit this stage and start a new trial
            % with all other parameters the same
            stagePreFixation    = 0;
            newTrialVariables 	= 0;
            addTrialDataFlag = 0;
        end
        % We wait 1 ms each loop-iteration so that we
        % don't overload the system in realtime-priority:
        WaitSecs(0.001);
    end
    
    
    
    
    % ****************************************************************************************
    %            FIXATION STAGE
    % ****************************************************************************************
    % Subject has begun fixating, can either hold fixation until the
    % targets come on and advance to stageTargetOn, or can abort the trial
    % and start a new trial by not maintaining fixation
    if stageFixation
        
        iTrial                  = iTrial + 1;
        % Prepare the screen for target onsets
        Screen('FillRect', window, fixationColor, fixationSquare);
        Screen('FillRect', window, targetColor, rightTargetSquare);
        Screen('FillRect', window, targetColor, leftTargetSquare);
        %                 Screen('FrameRect', window, targetColor, iTargetWindow);
        %                 Screen('FrameRect', window, targetColor, iDistractorWindow);
        %                 Screen('FrameRect', window, targetColor, fixationWindowPixel);
        tic
        while stageFixation
            
            % Unit code for checking eye position
            %-----------------------------------------------------------------------------------------------------
            if dummymode == 0
                error = Eyelink('CheckRecording');
                if error
                    fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
%                     runningTask = 0;
%                     cleanup;
iTrialOutcome = 'eyelinkError';
stageFeedback = 1;
break
                end
                
                if Eyelink( 'NewFloatSampleAvailable') > 0
                    % get the sample in the form of an event structure
                    event = Eyelink( 'NewestFloatSample');
                    event.gx;
                    event.gy;
                    if eyeRecorded ~= -1 % do we know which eye to use yet?
                        % if we do, get current gaze position from sample
                        x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                        y = event.gy(eyeRecorded+1);
                        % do we have valid data and is the pupil visible?
                        if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                            eyeX = x;
                            eyeY = y;
                            
                            %                             % Show where gaze is (troubleshooting)
                            %                             gazeRect = [ x-5 y-5 x+6 y+6];
                            %                             colour = round(rand(3,1)*255); % coloured dot
                            %                             Screen('FillOval', window,colour, gazeRect);
                            %                             Screen('FillRect', window, fixationColor, fixationSquare);
                            %                             Screen('FrameRect', window, targetColor, fixationWindow);
                            %                             [~, fixationOnsetTime, ~, ~, ~] = Screen('Flip', window);
                        end
                    end
                end
            else
                
                % Query current mouse cursor position (our "pseudo-eyetracker") -
                % (mx,my) is our gaze position.
                [eyeX, eyeY, ~] = GetMouse(window); %#ok<*NASGU>
            end
            %-----------------------------------------------------------------------------------------------------
            
            if ~inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc < iPreTargetFixDuration
                Eyelink('Message', 'Fixation_Abort');
                iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
                iTrialOutcome       = 'fixationAbort';
                stageFixation       = 0;
                stageFeedback       = 1;
                newTrialVariables            = 0;
            elseif inAcceptWindow(eyeX, eyeY, fixationWindowPixel) &&  toc >= iPreTargetFixDuration
                stageFixation       = 0;
                stageTargetOn       = 1;
            end
        end
    end
    
    
    
    % ****************************************************************************************
    %            TARGETS ON STAGE
    % ****************************************************************************************
    % Subject has obtained fixation, and the targets appeared,
    % can either conitue to fixate until choice stimulus comes on, or can abort the trial
    % and start a new trial by not maintaining fixation
    if stageTargetOn
        [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
        Eyelink('Message', 'Targets_On');
        iTargetOnset            = StimulusOnsetTime - iTrialOnsetComputerTime;
        iDistractorOnset        = StimulusOnsetTime - iTrialOnsetComputerTime;
        
        % Prepare the screen for choice stimulus onset
        Screen('FillRect', window, targetColor, rightTargetSquare);
        Screen('FillRect', window, targetColor, leftTargetSquare);
        Screen('FillRect', window, iStimulusColorsArray, stimulusSquaresArray);
        Screen('FillRect', window, fixationColor, fixationSquare);
        Screen('FillRect', window, goSignalColor, goSignalSquare);
        %        Screen('FrameRect', window, targetColor, iTargetWindow);
        %         Screen('FrameRect', window, targetColor, iDistractorWindow);
        %         Screen('FrameRect', window, targetColor, fixationWindowPixel);
        tic
        while stageTargetOn
            
            % Unit code for checking eye position
            %-----------------------------------------------------------------------------------------------------
            if dummymode == 0
                error = Eyelink('CheckRecording');
                if error
                    fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
%                     runningTask = 0;
%                     cleanup;
iTrialOutcome = 'eyelinkError';
stageFeedback = 1;
break
                end
                
                if Eyelink( 'NewFloatSampleAvailable') > 0
                    % get the sample in the form of an event structure
                    event = Eyelink( 'NewestFloatSample');
                    event.gx;
                    event.gy;
                    if eyeRecorded ~= -1 % do we know which eye to use yet?
                        % if we do, get current gaze position from sample
                        x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                        y = event.gy(eyeRecorded+1);
                        % do we have valid data and is the pupil visible?
                        if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                            eyeX = x;
                            eyeY = y;
                            
                            %                             % Show where gaze is (troubleshooting)
                            %                             gazeRect = [ x-5 y-5 x+6 y+6];
                            %                             colour = round(rand(3,1)*255); % coloured dot
                            %                             Screen('FillOval', window,colour, gazeRect);
                            %                             Screen('FillRect', window, fixationColor, fixationSquare);
                            %                             Screen('FrameRect', window, targetColor, fixationWindow);
                            %                             [~, fixationOnsetTime, ~, ~, ~] = Screen('Flip', window);
                            
                        end
                    end
                end
            else
                
                % Query current mouse cursor position (our "pseudo-eyetracker") -
                % (mx,my) is our gaze position.
                [eyeX, eyeY, ~] = GetMouse(window); %#ok<*NASGU>
            end
            %-----------------------------------------------------------------------------------------------------
            
            
            if ~inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc < iPostTargetFixDuration
                Eyelink('Message', 'Fixation_Abort');
                iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
                iTrialOutcome       = 'fixationAbort';
                stageTargetOn       = 0;
                stageFeedback       = 1;
                newTrialVariables            = 0;
            elseif inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc >= iPostTargetFixDuration
                stageTargetOn       = 0;
                stageChoiceOn       = 1;
            end
        end
    end
    
    
    
    
    if stageChoiceOn
        [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
        iFixationSpotDuration   = StimulusOnsetTime - iTrialOnsetComputerTime - iFixationSpotOnset;
        iChoiceStimulusOnset    = StimulusOnsetTime - iTrialOnsetComputerTime;
        iResponseCueOnset       = StimulusOnsetTime - iTrialOnsetComputerTime;
        Eyelink('Message', 'Choice_Stimulus_On');
        
        preStopSignalSaccade = 0;
        tic
        switch iTrialType
            case 'go'
                % **********************************************************************************************
                %         TARGETS AND STIMULI ON STAGE
                % **********************************************************************************************
                % All stimuli come on at once (in this version at least). In the case
                % of a go trial, wait until subject makes a saccade to one of the
                % targets or aborts by timing out. In the case of a stop trial, wait
                % until subject makes a saccade (an error) or
                while stageChoiceOn
                    
                    
                    % Unit code for checking eye position
                    %-----------------------------------------------------------------------------------------------------
                    if dummymode == 0
                        error = Eyelink('CheckRecording');
                        if error
                            fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
%                     runningTask = 0;
%                     cleanup;
iTrialOutcome = 'eyelinkError';
stageFeedback = 1;
break
                        end
                        
                        if Eyelink( 'NewFloatSampleAvailable') > 0
                            % get the sample in the form of an event structure
                            event = Eyelink( 'NewestFloatSample');
                            event.gx;
                            event.gy;
                            if eyeRecorded ~= -1 % do we know which eye to use yet?
                                % if we do, get current gaze position from sample
                                x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                                y = event.gy(eyeRecorded+1);
                                % do we have valid data and is the pupil visible?
                                if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                                    eyeX = x;
                                    eyeY = y;
                                    
                                    %                                     % Show where gaze is (troubleshooting)
                                    %                                     gazeRect = [ x-5 y-5 x+6 y+6];
                                    %                                     colour = round(rand(3,1)*255); % coloured dot
                                    %                                     Screen('FillOval', window,colour, gazeRect);
                                    %                                     Screen('FillRect', window, fixationColor, fixationSquare);
                                    %                                     Screen('FrameRect', window, targetColor, fixationWindow);
                                    %                                     [~, fixationOnsetTime, ~, ~, ~] = Screen('Flip', window);
                                    
                                end
                            end
                        end
                    else
                        
                        % Query current mouse cursor position (our "pseudo-eyetracker") -
                        % (mx,my) is our gaze position.
                        [eyeX, eyeY, ~] = GetMouse(window); %#ok<*NASGU>
                    end
                    %-----------------------------------------------------------------------------------------------------
                    
                    
                    % Made a saccade
                    if ~inAcceptWindow(eyeX, eyeY, fixationWindowPixel)
                        Eyelink('Message', 'Response_Onset');
                        iResponseOnset = GetSecs - iTrialOnsetComputerTime;
                        stageChoiceOn         = 0;
                        stageInFlight       = 1;
                        
                        % Waited too long to make a saccade
                    elseif inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc > graceResponse
                        Eyelink('Message', 'Response_Timed_Out');
                        iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
                        iTrialOutcome       = 'goIncorrect';
                        stageChoiceOn         = 0;
                        stageFeedback       = 1;
                    end
                end
                
                
            case 'stop'
                stopSignalOn = 0;
                iSSDIndex   = staircase(lastStopOutcomeArray{iProportionRightTargIndex}, iSSDIndexArray(iProportionRightTargIndex), nSSD);
                iSSD        = SSDArray(iSSDIndex);
                iSSDIndexArray(iProportionRightTargIndex) = iSSDIndex;   % Update the array to staircase SSDs in each discriminability level
                
                % **********************************************************************************************
                %         STOP SIGNAL STAGE
                % **********************************************************************************************
                Screen('FillRect', window, targetColor, rightTargetSquare);
                Screen('FillRect', window, targetColor, leftTargetSquare);
                Screen('FillRect', window, iStimulusColorsArray, stimulusSquaresArray);
                Screen('FillRect', window, fixationColor, fixationSquare);
                Screen('FillRect', window, stopSignalColor, stopSignalLocation);
                
                while stageChoiceOn
                    
                    % Unit code for checking eye position
                    %-----------------------------------------------------------------------------------------------------
                    if dummymode == 0
                        error = Eyelink('CheckRecording');
                        if error
                            fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
%                     runningTask = 0;
%                     cleanup;
iTrialOutcome = 'eyelinkError';
stageFeedback = 1;
break
                        end
                        
                        if Eyelink( 'NewFloatSampleAvailable') > 0
                            % get the sample in the form of an event structure
                            event = Eyelink( 'NewestFloatSample');
                            event.gx;
                            event.gy;
                            if eyeRecorded ~= -1 % do we know which eye to use yet?
                                % if we do, get current gaze position from sample
                                x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                                y = event.gy(eyeRecorded+1);
                                % do we have valid data and is the pupil visible?
                                if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                                    eyeX = x;
                                    eyeY = y;
                                    
                                    
                                end
                            end
                        end
                    else
                        
                        % Query current mouse cursor position (our "pseudo-eyetracker") -
                        % (mx,my) is our gaze position.
                        [eyeX, eyeY, ~] = GetMouse(window); %#ok<*NASGU>
                    end
                    %-----------------------------------------------------------------------------------------------------
                    
                    
                    % turn on the Stop Signal when appropriate
                    if ~stopSignalOn && toc >= iSSD
                        [~, stopOnsetTime, ~, ~, ~] = Screen('Flip', window);
                        stopSignalOn = 1;
                        Eyelink('Message', 'Stop_Signal_On');
                        iRealSSD    = stopOnsetTime - iTrialOnsetComputerTime - iChoiceStimulusOnset;
                        iStopOnset  = stopOnsetTime - iTrialOnsetComputerTime;
                        stageChoiceOn = 0;
                        stageStopSignalOn = 1;
                    elseif ~stopSignalOn && ~inAcceptWindow(eyeX, eyeY, fixationWindowPixel)
                        Eyelink('Message', 'Response_Onset');
                        iResponseOnset = GetSecs - iTrialOnsetComputerTime;
                        stageChoiceOn         = 0;
                        stageInFlight       = 1;
                        preStopSignalResponse = 1;
                        while toc < iSSD
                        end
                        [~, stopOnsetTime, ~, ~, ~] = Screen('Flip', window);
                        stopSignalOn = 1;
                        Eyelink('Message', 'Stop_Signal_On');
                        iRealSSD    = stopOnsetTime - iTrialOnsetComputerTime - iChoiceStimulusOnset;
                        iStopOnset  = stopOnsetTime - iTrialOnsetComputerTime;
                    end
                end
                
                while stageStopSignalOn
                    
                    
                    % Unit code for checking eye position
                    %-----------------------------------------------------------------------------------------------------
                    if dummymode == 0
                        error = Eyelink('CheckRecording');
                        if error
                            fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
%                     runningTask = 0;
%                     cleanup;
iTrialOutcome = 'eyelinkError';
stageFeedback = 1;
break
                        end
                        
                        if Eyelink( 'NewFloatSampleAvailable') > 0
                            % get the sample in the form of an event structure
                            event = Eyelink( 'NewestFloatSample');
                            event.gx;
                            event.gy;
                            if eyeRecorded ~= -1 % do we know which eye to use yet?
                                % if we do, get current gaze position from sample
                                x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                                y = event.gy(eyeRecorded+1);
                                % do we have valid data and is the pupil visible?
                                if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                                    eyeX = x;
                                    eyeY = y;
                                    
                                    
                                end
                            end
                        end
                    else
                        
                        % Query current mouse cursor position (our "pseudo-eyetracker") -
                        % (mx,my) is our gaze position.
                        [eyeX, eyeY, ~] = GetMouse(window); %#ok<*NASGU>
                    end
                    %-----------------------------------------------------------------------------------------------------
                    
                    
                    % Made a saccade
                    if ~inAcceptWindow(eyeX, eyeY, fixationWindowPixel)
                        Eyelink('Message', 'Response_Onset');
                        iResponseOnset = GetSecs - iTrialOnsetComputerTime;
                        stageStopSignalOn         = 0;
                        stageInFlight       = 1;
                    elseif inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc > iSSD + stopHoldDuration
                        iTrialOutcome       = 'stopCorrect';
                        lastStopOutcomeArray{iProportionRightTargIndex} = iTrialOutcome;
                        stageStopSignalOn         = 0;
                        stageFeedback       = 1;
                    end
                end
        end
        newTrialVariables = 1;
    end
    
    
    
    
    
    % **********************************************************************************************
    %         IN FLIGHT STAGE
    % **********************************************************************************************
    tic
    while stageInFlight
        
        
        % Unit code for checking eye position
        %-----------------------------------------------------------------------------------------------------
        if dummymode == 0
            error = Eyelink('CheckRecording');
            if error
                fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
%                     runningTask = 0;
%                     cleanup;
iTrialOutcome = 'eyelinkError';
stageFeedback = 1;
break
            end
            
            if Eyelink( 'NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                event = Eyelink( 'NewestFloatSample');
                event.gx;
                event.gy;
                if eyeRecorded ~= -1 % do we know which eye to use yet?
                    % if we do, get current gaze position from sample
                    x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                    y = event.gy(eyeRecorded+1);
                    % do we have valid data and is the pupil visible?
                    if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                        eyeX = x;
                        eyeY = y;
                        
                        
                    end
                end
            end
        else
            
            % Query current mouse cursor position (our "pseudo-eyetracker") -
            % (mx,my) is our gaze position.
            [eyeX, eyeY, ~] = GetMouse(window); %#ok<*NASGU>
        end
        %-----------------------------------------------------------------------------------------------------
        
        
        if toc > graceSaccadeDuration
            iTrialOutcome       = 'saccadeAbort';
            iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
            stageInFlight       = 0;
            stageFeedback       = 1;
        elseif inAcceptWindow(eyeX, eyeY, iTargetWindow) && toc <= graceSaccadeDuration
            iTargetWindowEntered = GetSecs - iTrialOnsetComputerTime;
            stageInFlight       = 0;
            stageOnTarget       = 1;
        elseif inAcceptWindow(eyeX, eyeY, iDistractorWindow) && toc <= graceSaccadeDuration
            iTargetWindowEntered = GetSecs - iTrialOnsetComputerTime;
            stageInFlight       = 0;
            stageOnDistractor       = 1;
        elseif inAcceptWindow(eyeX, eyeY, stimulusWindow) && toc <= graceSaccadeDuration
            % Flag it was go, but to the choice stimulus
            iTrialOutcome = 'choiceStimulusAbort';
            stageInFlight       = 0;
            stageFeedback       = 1;
        end
    end
    
    
    
    % **********************************************************************************************
    %         ON TARGET STAGE
    % **********************************************************************************************
    tic
    while stageOnTarget
        
        
        % Unit code for checking eye position
        %-----------------------------------------------------------------------------------------------------
        if dummymode == 0
            error = Eyelink('CheckRecording');
            if error
                fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
%                     runningTask = 0;
%                     cleanup;
iTrialOutcome = 'eyelinkError';
stageFeedback = 1;
break
            end
            
            if Eyelink( 'NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                event = Eyelink( 'NewestFloatSample');
                event.gx;
                event.gy;
                if eyeRecorded ~= -1 % do we know which eye to use yet?
                    % if we do, get current gaze position from sample
                    x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                    y = event.gy(eyeRecorded+1);
                    % do we have valid data and is the pupil visible?
                    if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                        eyeX = x;
                        eyeY = y;
                        
                        
                    end
                end
            end
        else
            
            % Query current mouse cursor position (our "pseudo-eyetracker") -
            % (mx,my) is our gaze position.
            [eyeX, eyeY, ~] = GetMouse(window); %#ok<*NASGU>
        end
        %-----------------------------------------------------------------------------------------------------
        
        
        switch iTrialType
            case 'go'
                if inAcceptWindow(eyeX, eyeY, iTargetWindow) && toc >= postSaccadeHoldDuration
                    % Flag it was correct to the target
                    iTrialOutcome = 'goCorrectTarget';
                    stageOnTarget = 0;
                    stageFeedback = 1;
                elseif ~inAcceptWindow(eyeX, eyeY, iTargetWindow) && toc < postSaccadeHoldDuration
                    % Flag it was target aborted
                    iTrialOutcome = 'targetHoldAbort';
                    stageOnTarget = 0;
                    stageFeedback = 1;
                end
            case 'stop'
                if inAcceptWindow(eyeX, eyeY, iTargetWindow) && toc >= postSaccadeHoldDuration
                    % Flag it was correct to the target
                    if ~preStopSignalResponse
                        iTrialOutcome = 'stopIncorrectTarget';
                    elseif preStopSignalResponse
                        iTrialOutcome = 'stopIncorrectPreSSDTarget';
                    end
                    stageOnTarget = 0;
                    stageFeedback = 1;
                    lastStopOutcomeArray{iProportionRightTargIndex} = iTrialOutcome;
                elseif ~inAcceptWindow(eyeX, eyeY, iTargetWindow) && toc < postSaccadeHoldDuration
                    % Flag it was target aborted
                    iTrialOutcome = 'targetHoldAbort';
                    stageOnTarget = 0;
                    stageFeedback = 1;
                    lastStopOutcomeArray{iProportionRightTargIndex} = iTrialOutcome;
                end
        end
    end
    
    
    
    % **********************************************************************************************
    %         ON DISTRACTOR STAGE
    % **********************************************************************************************
    tic
    while stageOnDistractor
        
        
        % Unit code for checking eye position
        %-----------------------------------------------------------------------------------------------------
        if dummymode == 0
            error = Eyelink('CheckRecording');
            if error
                fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
%                     runningTask = 0;
%                     cleanup;
iTrialOutcome = 'eyelinkError';
stageFeedback = 1;
break
            end
            
            if Eyelink( 'NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                event = Eyelink( 'NewestFloatSample');
                event.gx;
                event.gy;
                if eyeRecorded ~= -1 % do we know which eye to use yet?
                    % if we do, get current gaze position from sample
                    x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                    y = event.gy(eyeRecorded+1);
                    % do we have valid data and is the pupil visible?
                    if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                        eyeX = x;
                        eyeY = y;
                        
                        
                    end
                end
            end
        else
            
            % Query current mouse cursor position (our "pseudo-eyetracker") -
            % (mx,my) is our gaze position.
            [eyeX, eyeY, ~] = GetMouse(window); %#ok<*NASGU>
        end
        %-----------------------------------------------------------------------------------------------------
        
        
        switch iTrialType
            case 'go'
                if inAcceptWindow(eyeX, eyeY, iDistractorWindow) && toc >= postSaccadeHoldDuration
                    % Flag it was correct to the target
                    iTrialOutcome = 'goCorrectDistractor';
                    stageOnDistractor = 0;
                    stageFeedback = 1;
                elseif ~inAcceptWindow(eyeX, eyeY, iDistractorWindow) && toc < postSaccadeHoldDuration
                    % Flag it was target aborted
                    iTrialOutcome = 'distractorHoldAbort';
                    stageOnDistractor = 0;
                    stageFeedback = 1;
                end
            case 'stop'
                if inAcceptWindow(eyeX, eyeY, iDistractorWindow) && toc >= postSaccadeHoldDuration
                    % Flag it was correct to the target
                    if ~preStopSignalResponse
                        iTrialOutcome = 'stopIncorrectDistractor';
                    elseif preStopSignalResponse
                        iTrialOutcome = 'stopIncorrectPreSSDDistractor';
                    end
                    stageOnDistractor = 0;
                    stageFeedback = 1;
                    lastStopOutcomeArray{iProportionRightTargIndex} = iTrialOutcome;
                elseif ~inAcceptWindow(eyeX, eyeY, iDistractorWindow) && toc < postSaccadeHoldDuration
                    % Flag it was target aborted
                    iTrialOutcome = 'distractorHoldAbort';
                    stageOnDistractor = 0;
                    stageFeedback = 1;
                    lastStopOutcomeArray{iProportionRightTargIndex} = iTrialOutcome;
                end
        end
    end
    
    
    
    % **********************************************************************************************
    %         FEEDBACK STAGE
    % **********************************************************************************************
    while stageFeedback
        if iTrial > 20
            feedbackTime = .2;
        end
        [window, timeout] = feedback(window, iTrialOutcome, iTrialType);
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
        stageFeedback = 0;
    end
    
    iTrialDuration = GetSecs - iTrialOnsetComputerTime;
    Eyelink('StopRecording');
    
    
    
    % Allow option to quit task by pressing escape at the end of a trial
    tic
    while toc < 1
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        if keyIsDown && find(keyCode) == escapeKey
            runningTask = 0;
        end
    end
    
    
    
    if  addTrialDataFlag
        % Add the trial's variables to the data set.
        
        % Event Timing
        % ------------
        trialOnset              = [trialOnset; iTrialOnset];
        trialDuration           = [trialDuration; iTrialDuration];
        fixationSpotOnset       = [fixationSpotOnset; iFixationSpotOnset];
        fixationOnset           = [fixationOnset; iFixationOnset];
        fixationSpotDuration    = [fixationSpotDuration; iFixationSpotDuration];
        targetOnset             = [targetOnset; iTargetOnset];
        targetDuration          = [targetDuration; iTargetDuration];
        distractorOnset         = [distractorOnset; iDistractorOnset];
        distractorDuration      = [distractorDuration; iDistractorDuration];
        choiceStimulusOnset     = [choiceStimulusOnset; iChoiceStimulusOnset];
        choiceStimulusDuration  = [choiceStimulusDuration; iChoiceStimulusDuration];
        targetWindowEntered     = [targetWindowEntered; iTargetWindowEntered];
        SSD                     = [SSD; iSSD];
        realSSD                 = [realSSD; iRealSSD];
        stopOnset               = [stopOnset; iStopOnset];
        stopDuration            = [stopDuration; iStopDuration];
        responseCueOnset        = [responseCueOnset; iResponseCueOnset];
        responseOnset           = [responseOnset; iResponseOnset];
        feedbackOnset           = [feedbackOnset; iFeedbackOnset];
        feedbackDuration        = [feedbackDuration; iFeedbackDuration];
        abortOnset              = [abortOnset; iAbortOnset];
        trialOutcome            = [trialOutcome; iTrialOutcome];
        
        
        % Stimulus Properties
        % ------------
        choiceStimulusColor     = [choiceStimulusColor; iStimulusColorsArray];
        checkerboardArray       = [checkerboardArray; iCheckerboardArray];
        target1CheckerProportion = [target1CheckerProportion; iProportionRightTarg];
        targetAngle          = [targetAngle; iTargetAngle];
        distractorAngle      = [distractorAngle; iDistractorAngle];
    end
    %     SSD = [SSD; iSSD];
    %     realSSD = [realSSD; iRealSSD];
    
    
    
    if iTrial == totalTrial
        runningTask = 0;
    end
    
end % trial loop


% ********************************************************************
%                    End Trial Loop
% ********************************************************************


nTrial = length(trialOnset);









ssdTrack = SSD / flipTime;
ssdTrack(isnan(ssdTrack)) = [];
figure(67)
clf
plot(1:length(ssdTrack), ssdTrack)












% ********************************************************************
%                    Trial Data
% ********************************************************************


% Event Timing
% ------------




% Event Properties
% ---------------------------------------------------------------
stopTrialProportion  	= ones(nTrial, 1) .* stopTrialProportion;



% Location of Stimuli
% ---------------------------------------------------------------
fixationWindow = ones(nTrial, 1) * fixationWindow;
targetWindow = ones(nTrial, 1) * targetWindow;
distractorWindow = ones(nTrial, 1) * distractorWindow;
choiceStimulusWindow = ones(nTrial, 1) * choiceStimulusWindow;

fixationAngle           = ones(nTrial, 1) * fixationAngle;
fixationAmplitude     	= ones(nTrial, 1) * fixationAmplitude;
choiceStimulusAngle   	= ones(nTrial, 1) * choiceStimulusAngle;
choiceStimulusAmplitude	= ones(nTrial, 1) * choiceStimulusAmplitude;
targetAmplitude =  ones(nTrial, 1) * targetAmplitude;
distractorAmplitude	= ones(nTrial, 1) * targetAmplitude(1);  % for now target and distractor are equidistant from fixation

fixationSize              = ones(nTrial, 1) * fixationWidth;
targetSize                = ones(nTrial, 1) * targetWidth;
distractorSize            = ones(nTrial, 1) * distractorWidth;
choiceStimulusSize        = ones(nTrial, 1) * squareWidth * stimulusColumns;
fixationColor             = ones(nTrial, 1) * fixationColor;
targetColor               = ones(nTrial, 1) * targetColor;
distractorColor           = ones(nTrial, 1) * distractorColor;

% size(trialOutcome)
% size(trialOnset)
% size(trialDuration)
% size(abortOnset)
% size(fixationSpotOnset)
% size(fixationSpotDuration)
% size(targetOnset)
% size(targetDuration)
% size(distractorOnset)
% size(distractorDuration)
% size(choiceStimulusOnset)
% size(choiceStimulusDuration)
% size(stopOnset)
% size(stopDuration)
% size(SSD)
% size(realSSD)
% size(responseCueOnset)
% size(responseOnset)
% size(targetWindowEntered)
% size(feedbackOnset)
% size(feedbackDuration)
% size(stopTrialProportion)
% size(fixationAngle)
% size(targetAngle)
% size(distractorAngle)
% size(choiceStimulusAngle)
% size(fixationAmplitude)
% size(targetAmplitude)
% size(distractorAmplitude)
% size(choiceStimulusAmplitude)
% size(fixationSize)
% size(targetSize)
% size(distractorSize)
% size(choiceStimulusSize)
% size(fixationWindow)
% size(targetWindow)
% size(distractorWindow)
% size(choiceStimulusWindow)
% size(fixationColor)
% size(targetColor)
% size(distractorColor)
% size(target1CheckerProportion)
% size(choiceStimulusColor)
% Screen('CloseAll');

trialData = dataset(...
    {trialOutcome,              'trialOutcome'},...
    {round(trialOnset*1000),           'trialOnset'},...
    {round(trialDuration*1000),         'trialDuration'},...
    {round(abortOnset*1000),            'abortOnset'},...
    {round(fixationSpotOnset*1000),     'fixationSpotOnset'},...
    {round(fixationOnset*1000),         'fixationOnset'},...
    {round(fixationSpotDuration*1000),  'fixationSpotDuration'},...
    {round(targetOnset*1000),           'targetOnset'},...
    {round(targetDuration*1000),        'targetDuration'},...
    {round(distractorOnset*1000),       'distractorOnset'},...
    {round(distractorDuration*1000),    'distractorDuration'},...
    {round(choiceStimulusOnset*1000),   'choiceStimulusOnset'},...
    {round(choiceStimulusDuration*1000), 'choiceStimulusDuration'},...
    {round(stopOnset*1000),             'stopOnset'},...
    {round(stopDuration*1000),          'stopDuration'},...
    {round(SSD*1000),                   'SSD'},...
    {round(realSSD*1000),               'realSSD'},...
    {round(responseCueOnset*1000),      'responseCueOnset'},...
    {round(responseOnset*1000),         'responseOnset'},...
    {round(targetWindowEntered*1000),	'targetWindowEntered'},...
    {round(feedbackOnset*1000),         'feedbackOnset'},...
    {round(feedbackDuration*1000),      'feedbackDuration'},...
    {stopTrialProportion,       'stopTrialProportion'},...
    {fixationAngle,             'fixationAngle'},...
    {targetAngle,               'targetAngle'},...
    {distractorAngle,           'distractorAngle'},...
    {choiceStimulusAngle,       'choiceStimulusAngle'},...
    {fixationAmplitude,         'fixationAmplitude'},...
    {targetAmplitude,           'targetAmplitude'},...
    {distractorAmplitude,       'distractorAmplitude'},...
    {choiceStimulusAmplitude,   'choiceStimulusAmplitude'},...
    {fixationSize,              'fixationSize'},...
    {targetSize,                'targetSize'},...
    {distractorSize,            'distractorSize'},...
    {choiceStimulusSize,        'choiceStimulusSize'},...
    {fixationWindow,            'fixationWindow'},...
    {targetWindow,              'targetWindow'},...
    {distractorWindow,          'distractorWindow'},...
    {choiceStimulusWindow,      'choiceStimulusWindow'},...
    {fixationColor,             'fixationColor'},...
    {targetColor,               'targetColor'},...
    {distractorColor,           'distractorColor'},...
    {target1CheckerProportion,  'target1CheckerProportion'},...
    {checkerboardArray,         'checkerboardArray'},...
    {choiceStimulusColor,       'choiceStimulusColor'});


% ********************************************************************
% Session Data
% ********************************************************************

SessionData.taskID = 'ccm';
sessionData.task.effector = 'eyeMovement';

sessionData.timing.year = num2str(year(now));
sessionData.timing.month = num2str(month(now));
sessionData.timing.day = num2str(day(now));
sessionData.timing.hour = num2str(clockVector(4));
sessionData.timing.minute = num2str(clockVector(5));

sessionData.timing.totalDuration = trialOnset(end) + trialDuration(end) - trialOnset(1); % seconds


sessionData.subjectID = subjectID;
sessionData.sessionID = session;

sessionData.stimuli.stopSignalRGB       = stopSignalColor;
sessionData.stimuli.target1CheckerRGB   = [rightTargetCheckerColor(1), rightTargetCheckerColor(1), rightTargetCheckerColor(1)];
sessionData.stimuli.target2CheckerRGB   = [leftTargetCheckerColor(1), leftTargetCheckerColor(1), leftTargetCheckerColor(1)];
sessionData.stimuli.fixationRGB         = fixationColor;
sessionData.stimuli.targetRGB           = targetColor;








% ---------- Window Cleanup ----------
cleanup
Eyelink('CloseFile');
status = Eyelink('ReceiveFile',edffilename, 'data/');
status = Eyelink('ReceiveFile',edffilename);




if plotFlag
    post_session_psychometric(trialData)
end





if saveFlag
    save(saveFileName, 'trialData', 'sessionData');
    success = copy_human_files(edffilename, saveFileName);
end








end  % main function
















% **************************************************************************************
%               SUBFUNCTIONS
% **************************************************************************************




function cleanup
Eyelink('stoprecording');
Screen('CloseAll');
ShowCursor;
end




% *******************************************************************
function [window, timeout] = feedback(window, iTrialOutcome)


scrsz = get(0, 'ScreenSize');
screenWidth = scrsz(3);
screenHeight = scrsz(4);

incorrectTextColor = [250, 50, 50];
correctTextColor = [50, 220, 50];
Screen('TextFont', window, 'Times');
Screen('TextSize', window, 30);
Screen('TextStyle', window, 1);
if strcmp(iTrialOutcome, 'goCorrectTarget') || strcmp(iTrialOutcome, 'stopCorrect')
    [nx, ny, bbox] = DrawFormattedText(window, 'Good', 'center', 'center', correctTextColor);
    timeout = 0;
elseif strcmp(iTrialOutcome, 'goCorrectDistractor')
    DrawFormattedText(window, 'Wrong Target', 'center', 'center', incorrectTextColor);
    timeout = .5;
elseif strcmp(iTrialOutcome, 'goIncorrect')
    DrawFormattedText(window, 'You should have responded', 'center', 'center', incorrectTextColor);
    timeout = 1;
elseif strcmp(iTrialOutcome, 'stopIncorrectTarget') || strcmp(iTrialOutcome, 'stopIncorrectDistractor')
    DrawFormattedText(window, 'Missed Stop', 'center', 'center', incorrectTextColor);
    timeout = .5;
elseif strcmp(iTrialOutcome, 'stopIncorrectPreSSDTarget') || strcmp(iTrialOutcome, 'stopIncorrectPreSSDDistractor')
    DrawFormattedText(window, 'Missed Stop', 'center', 'center', [190, 190, 190]);
    timeout = .5;
elseif strcmp(iTrialOutcome, 'fixationAbort')
    DrawFormattedText(window, 'Please stay fixated on center spot', 'center', 'center', [190, 190, 190]);
    timeout = 1.5;
elseif strcmp(iTrialOutcome, 'saccadeAbort')
    DrawFormattedText(window, 'Your eyes wandered from center', 'center', 'center', [190, 190, 190]);
    timeout = 1;
elseif strcmp(iTrialOutcome, 'choiceStimulusAbort')
    DrawFormattedText(window, 'Don"t look up at the checkerboard', 'center', 'center', incorrectTextColor);
    timeout = 2;
elseif strcmp(iTrialType, 'stop') && (strcmp(iTrialOutcome, 'targetHoldAbort') || strcmp(iTrialOutcome, 'distractorHoldAbort'))
    DrawFormattedText(window, 'Missed Stop', 'center', 'center', incorrectTextColor);
    timeout = .2;
elseif strcmp(iTrialType, 'stop') && (strcmp(iTrialOutcome, 'targetHoldAbort') || strcmp(iTrialOutcome, 'distractorHoldAbort'))
    DrawFormattedText(window, 'Keep your eyes on the target until it disappears ', 'center', 'center', incorrectTextColor);
    timeout = .2;
elseif strcmp(iTrialOutcome, 'eyelinkError')
    DrawFormattedText(window, 'My fault- continue on', 'center', 'center', incorrectTextColor);
    timeout = 0;
end

% Screen('Flip', window);
% waitSecs(feedbackTime + timeout)
end


% *******************************************************************
function newSSDIndex = staircase(lastStopOutcome, lastSSDIndex, nSSD)

maxStepSize = 3;
iStepSize = randi(maxStepSize);

if strcmp(lastStopOutcome, 'stopCorrect')
    newSSDIndex = min(iStepSize + lastSSDIndex, nSSD);
elseif isempty(lastStopOutcome) ||  strcmp(lastStopOutcome, 'stopIncorrectTarget') || strcmp(lastStopOutcome, 'stopIncorrectDistractor') || ...
        strcmp(lastStopOutcome, 'stopIncorrectPreSSDTarget') || strcmp(lastStopOutcome, 'stopIncorrectPreSSDDistractor') || ...
        strcmp(lastStopOutcome, 'targetHoldAbort') || strcmp(lastStopOutcome, 'distractorHoldAbort')
    newSSDIndex = max(lastSSDIndex - iStepSize, 1);
end
end



% % *******************************************************************
% function [eyeX, eyeY] = locateEyePosition(el, eyeRecorded, dummymode)
% if dummymode == 0
%     error = Eyelink('CheckRecording');
%
%     if Eyelink( 'NewFloatSampleAvailable') > 0
%         % get the sample in the form of an event structure
%         event = Eyelink( 'NewestFloatSample');
%         event.gx
%         event.gy
%         if eyeRecorded ~= -1 % do we know which eye to use yet?
%             % if we do, get current gaze position from sample
%             x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
%             y = event.gy(eyeRecorded+1);
%             % do we have valid data and is the pupil visible?
%             if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
%                 eyeX = x;
%                 eyeY = y;
%             end
%         end
%     end
% else
%
%     % Query current mouse cursor position (our "pseudo-eyetracker") -
%     % (mx,my) is our gaze position.
%     [eyeX, eyeY, ~]=GetMouse(window); %#ok<*NASGU>
% end
%
% end

% *******************************************************************
function inWindow = inAcceptWindow(eyeX, eyeY, acceptWindow)
% determine if gx and gy are within fixation window
inWindow = eyeX > acceptWindow(1) &&  eyeX <  acceptWindow(3) && ...
    eyeY > acceptWindow(2) && eyeY < acceptWindow(4) ;
end










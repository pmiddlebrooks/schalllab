function choice_stop_task_images(SSDArrayScreenFlips, choiceStimulusAmplitude, choiceStimulusAngle, targetAmplitude, rightTargetAngle, plotFlag, saveFlag)


saveImages = 1;
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
    SSDArrayScreenFlips = 6 : 6 : 6 + 6*12;
    choiceStimulusAmplitude = 5;
    choiceStimulusAngle = 90;
    targetAmplitude = 12;
    rightTargetAngle = 0;
    plotFlag = 1;
    saveFlag = 1;
end

% rightTargProportionArray   = [0 1];
rightTargProportionArray     = [.48 .52];
% rightTargProportionArray     = [.41 .45 .48 .5 .52 .55 .59];
% rightTargProportionArray     = [.35 .42 .46 .5 .54 .58 .65];
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
pixelPerDegree     = screenWidthPixel / (2*theta);

targetAmplitudePixel = targetAmplitude * pixelPerDegree;
stimulusAmplitudePixel = choiceStimulusAmplitude * pixelPerDegree;



whichScreen     = 0;
% backGround      = [40 40 40];
backGround      = [20 20 20];
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
stopTrialProportion     = 1;
trialsPerBlock          = 5;
BLOCKS_TO_RUN           = 1;
totalTrial              = trialsPerBlock * BLOCKS_TO_RUN;
SSDArray                = SSDArrayScreenFlips * (flipTime); %subtracting 2 ms each cycle to enusre the SSD occurs before the next screen refresh
graceResponse           = 1.7; % seconds allowed to make a saccade
graceObtainFixation     = 2;  % seconds program will wait for user to obtain fixation
stopHoldDuration        = 1;
graceSaccadeDuration    = .1; % seconds allowed intra saccade time
postSaccadeHoldDuration = .4;  % duration to hold post-saccade fixation
stimulusScaleConversion = 1/ 10; % for now 10 pixels per every 100 pixels from fixation
feedbackTime            = .5;
dummymode               = 0;



% Keyboard assignments
KbName('UnifyKeyNames');
rightTargetKey          = KbName('m');
leftTargetKey           = Kbname('z');
escapeKey               = KbName('ESCAPE');


% Fixation Spot constants
fixationAngle       = 0;
fixationAmplitude   = 0;
% fixationWidth       = .5;
fixationWidth       = 1.5;
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
% goSignalColor       = [50 50 50];
goSignalColor       = backGround;

% Stop signal constants
stopSignalWidthPixel     = fixationWidthPixel - 2;
stopSignalColor     = [255 150 0];


% targetWidth         = .5 + (targetAmplitude * stimulusScaleConversion);
% targetWidth         = .7;
targetWidth         = 2;
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
targetColor         = [160 160 160];
distractorColor     = targetColor;

% Checkered Stimulus constants
stimulusColumns     = 10;
stimulusRows        = 10;
nSquare            = stimulusColumns * stimulusRows;
% squareWidthPixel         = 3; % pixels
squareWidthPixel         = 6; % pixels
squareWidth = squareWidthPixel / pixelPerDegree;














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
Screen('FillRect', window, backGround);
[vbl SOT] = Screen('Flip', window);
if saveImages
    image = Screen('GetImage', window);
    imwrite(image, 'blankScreen.tiff', 'tiff');
end
waitSecs(1)



runningTask = 1;
iTrial = 1;
taskStartTime = GetSecs;
newTrialVariables = 1; % A flag, set to zero if an abort occurs and we want all trial variables to stay the same
iBlock = 1;

% for iTrial = 0 : 2
while runningTask
    
    
    
    
    
    
    
    
    
    
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
    
    
    
    
    iTrialOnset             = GetSecs - taskStartTime;
    iTrialOnsetComputerTime = GetSecs;
    
    
    
    
    
    
    iFixationSpotOnset = fixationOnsetTime - iTrialOnsetComputerTime;
    
    if stagePreFixation
        
        if saveImages
            image = Screen('GetImage', window);
            imwrite(image, 'fixation.tiff', 'tiff');
        end
        
        stagePreFixation    = 0;
        stageFixation       = 1;
    end
    
    
    
    
    % ****************************************************************************************
    %            FIXATION STAGE
    % ****************************************************************************************
    % Subject has begun fixating, can either hold fixation until the
    % targets come on and advance to stageTargetOn, or can abort the trial
    % and start a new trial by not maintaining fixation
    if stageFixation
        
        % Prepare the screen for target onsets
        Screen('FillRect', window, fixationColor, fixationSquare);
        Screen('FillRect', window, targetColor, rightTargetSquare);
        Screen('FillRect', window, targetColor, leftTargetSquare);
        %                 Screen('FrameRect', window, targetColor, iTargetWindow);
        %                 Screen('FrameRect', window, targetColor, iDistractorWindow);
        %                 Screen('FrameRect', window, targetColor, fixationWindowPixel);
        waitSecs(iPreTargetFixDuration)
        [~, fixationOnsetTime, ~, ~, ~] = Screen('Flip', window);
        if saveImages
            image = Screen('GetImage', window);
            imwrite(image, 'targetsOn.tiff', 'tiff');
        end
        stageFixation       = 0;
        stageTargetOn       = 1;
    end
    
    
    
    % ****************************************************************************************
    %            TARGETS ON STAGE
    % ****************************************************************************************
    % Subject has obtained fixation, and the targets appeared,
    % can either conitue to fixate until choice stimulus comes on, or can abort the trial
    % and start a new trial by not maintaining fixation
    if stageTargetOn
        
        % Prepare the screen for choice stimulus onset
        Screen('FillRect', window, targetColor, rightTargetSquare);
        Screen('FillRect', window, targetColor, leftTargetSquare);
        Screen('FillRect', window, iStimulusColorsArray, stimulusSquaresArray);
        Screen('FillRect', window, fixationColor, fixationSquare);
        Screen('FillRect', window, goSignalColor, goSignalSquare);
        %        Screen('FrameRect', window, targetColor, iTargetWindow);
        %         Screen('FrameRect', window, targetColor, iDistractorWindow);
        %         Screen('FrameRect', window, targetColor, fixationWindowPixel);
        waitSecs(iPostTargetFixDuration)
        [~, fixationOnsetTime, ~, ~, ~] = Screen('Flip', window);
        if saveImages
            image = Screen('GetImage', window);
            imwrite(image, 'checkerboardOn.tiff', 'tiff');
        end
        stageChoiceOn = 1;
    end
    
    
    
    
    if stageChoiceOn
        switch iTrialType
            case 'go'
                % **********************************************************************************************
                %         TARGETS AND STIMULI ON STAGE
                % **********************************************************************************************
                % All stimuli come on at once (in this version at least). In the case
                % of a go trial, wait until subject makes a saccade to one of the
                % targets or aborts by timing out. In the case of a stop trial, wait
                % until subject makes a saccade (an error) or
                goInd = randi(2);
                switch goInd
                    case 1
                        iTrialOutcome = 'goCorrectTarget';
                   case 2
                        iTrialOutcome = 'goCorrectDistractor';
                end
                     waitSecs(stopHoldDuration)
                   stageFeedback       = 1;
                
                
                
                
            case 'stop'
                stopSignalOn = 0;
                iSSDIndex   = randi(length(iSSDIndexArray(iProportionRightTargIndex)));
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
                
                waitSecs(iSSD)
                [~, fixationOnsetTime, ~, ~, ~] = Screen('Flip', window);
                if saveImages
                    image = Screen('GetImage', window);
                    imwrite(image, 'stopSignalOnset.tiff', 'tiff');
                    imwrite(image, 'stopSignalOnset.tiff', 'tiff');
                end
                
                
                stageStopSignalOn = 1;
                
                if stageStopSignalOn
                    
                    waitSecs(stopHoldDuration)
                    stageFeedback       = 1;
                end
                
                stopInd = randi(2);
                switch stopInd
                    case 1
                        iTrialOutcome = 'stopCorrect';
                    case 2
                        targInd = randi(2);
                    switch targInd
                        case 1
                        iTrialOutcome = 'stopIncorrectTarget';    
                        case 2
                        iTrialOutcome = 'stopIncorrectDistractor';   
                    end
                end
                
        end
       iTrial = iTrial + 1;
    end
    
    
    
    
    
        
        
    
    

    
    % **********************************************************************************************
    %         FEEDBACK STAGE
    % **********************************************************************************************
    while stageFeedback
            feedbackTime = .5;
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
    
    
    
    % Allow option to quit task by pressing escape at the end of a trial
    tic
    while toc < 1
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        if keyIsDown && find(keyCode) == escapeKey
            runningTask = 0;
        end
    end
    
    
    
     
    if iTrial == totalTrial
        runningTask = 0;
    end
    
    
   
end % trial loop


% ********************************************************************
%                    End Trial Loop
% ********************************************************************


nTrial = length(trialOnset);




















Screen('CloseAll');
ShowCursor;







end  % main function
















% **************************************************************************************
%               SUBFUNCTIONS
% **************************************************************************************







% *******************************************************************
function [window, timeout] = feedback(window, iTrialOutcome, iTrialType)


scrsz = get(0, 'ScreenSize');
screenWidth = scrsz(3);
screenHeight = scrsz(4);

incorrectTextColor = [250, 50, 50];
correctTextColor = [50, 220, 50];
Screen('TextFont', window, 'Times');
Screen('TextSize', window, 30);
Screen('TextStyle', window, 1);
timeout = 0;
if strcmp(iTrialOutcome, 'goCorrectTarget') || strcmp(iTrialOutcome, 'stopCorrect')
    [nx, ny, bbox] = DrawFormattedText(window, 'Good', 'center', 'center', correctTextColor);
    timeout = 0;
elseif strcmp(iTrialOutcome, 'goCorrectDistractor')
    DrawFormattedText(window, 'Wrong Target', 'center', 'center', incorrectTextColor);
    timeout = .2;
elseif strcmp(iTrialOutcome, 'goIncorrect')
    DrawFormattedText(window, 'You should have responded', 'center', 'center', incorrectTextColor);
    timeout = 1;
elseif strcmp(iTrialOutcome, 'stopIncorrectTarget') || strcmp(iTrialOutcome, 'stopIncorrectDistractor')
    DrawFormattedText(window, 'Missed Stop', 'center', 'center', incorrectTextColor);
    timeout = .2;
elseif strcmp(iTrialOutcome, 'stopIncorrectPreSSDTarget') || strcmp(iTrialOutcome, 'stopIncorrectPreSSDDistractor')
    DrawFormattedText(window, 'Missed Stop', 'center', 'center', [190, 190, 190]);
    timeout = .2;
elseif strcmp(iTrialOutcome, 'fixationAbort')
    DrawFormattedText(window, 'Please stay fixated on center spot', 'center', 'center', [190, 190, 190]);
    timeout = 1.5;
elseif strcmp(iTrialOutcome, 'saccadeAbort')
    DrawFormattedText(window, 'Your eyes wandered from center', 'center', 'center', [190, 190, 190]);
    timeout = .5;
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












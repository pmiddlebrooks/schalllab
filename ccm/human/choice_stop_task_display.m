function choice_stop_task_display(stopMode, effector)



% stimulusAmplitude: in some units (e.g. pixels), the distance from the
%       center of the screen to the center of the stimulus
% stimulusAngle: in degrees, the angle from the center of the screen to the
%       center of the stimulus

% Format (numbering) of the squares in the stimulus (a 3 X 3 example):
%             0   1   2
%             3   4   5
%             6   7   8

% example SSDArrayScreenFlips = [22 26 30 34]

if nargin < 3
    magentaTargSide = 'right';
    SSDArrayScreenFlips = 6 : 6 : 6 + 6*12;
    checkerAmp = 3;
    checkerAngle = 90;
    targAmp = 10;
    rightTargAngle = 0;
end

rightTargPropArray     = [.35 .42 .46 .5 .54 .58 .65];
rightTargetRate = .5; % How often should right side be target?
fiftyPercentRate = .6; % how often should 50% signal strength be presented RELATIVE TO OTHER proportions?

nProportion             = length(rightTargPropArray);
% INITIAL_SSD_INDEX       = 6;
% iSSDIndexArray          = INITIAL_SSD_INDEX * ones(nProportion, 1);
INITIAL_SSD_INDEX       = [6 7 8 8 8 7 6];
iSSDIndexArray          = INITIAL_SSD_INDEX;


if strcmp(stopMode, 'auditory')
    stopHz = 500;
    % stopHz = 750;
    % stopHz = 1000;
    stopWavFile = [num2str(stopHz), 'Hz.wav'];
else
    stopHz = nan;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Set up Experiment Variables, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get distance values and conversions for stimuli presentation
metersFromScreen    = .64;
screenWidthMeters   = .4;
theta               = asind(screenWidthMeters / 2 / sqrt(((screenWidthMeters / 2)^2) + metersFromScreen^2));
screenSize          = get(0, 'ScreenSize');
screenWidthPixel  	= screenSize(3);
screenHeightPixel 	= screenSize(4);
matlabCenterX       = screenWidthPixel/2;
matlabCenterY       = screenHeightPixel/2;
pixelPerDegree      = screenWidthPixel / (2*theta);

targetAmplitudePixel = targAmp * pixelPerDegree;
stimulusAmplitudePixel = checkerAmp * pixelPerDegree;




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


% ***************************************************************
%     CONSTANTS
% ***************************************************************
pStopTrial              = .8;
nTrialPerBlock          = 10;
nBlock                  = 1;
nTrialTotal              = nTrialPerBlock * nBlock;
SSDArray                = SSDArrayScreenFlips * (flipTime); %subtracting 2 ms each cycle to enusre the ssd occurs before the next screen refresh
graceResponse           = 1.7; % seconds allowed to make a saccade
graceObtainFix          = 2;  % seconds program will wait for user to obtain fixation
stopHoldDuration        = 1;
graceSaccDuration       = .1; % seconds allowed intra saccade time
postSaccHoldDuration    = .4;  % duration to hold post-saccade fixation
stimulusScaleConversion = 1/ 10; % for now 10 pixels per every 100 pixels from fixation
feedbackTime            = .1;
dummymode               = 0;


% Keyboard assignments
KbName('UnifyKeyNames');
rightTargKey          = KbName('m');
leftTargKey           = Kbname('z');
escapeKey               = KbName('ESCAPE');


% Fixation Spot constants
fixAngle            = 0;
fixAmp              = 0;
fixationWidth       = .5;
fixationWidthPixel  = fixationWidth * pixelPerDegree;
fixWindowScale      = 3;  % fix window is 4 times size of fix point
fixWindow           = fixationWidth*fixWindowScale;
fixWindowPixel      = [-fixationWidthPixel*fixWindowScale, -fixationWidthPixel*fixWindowScale, fixationWidthPixel*fixWindowScale, fixationWidthPixel*fixWindowScale];
fixWindowPixel      = CenterRect(fixWindowPixel, centerPoint);
fixationColor       = [200 200 200];
fixationHoldBase    = .5;
fixationHoldAdd     = (0 : 10 : 500)  / 1000;



% targetWidth         = .5 + (targAmp * stimulusScaleConversion);
targetWidth         = 1;
targetWidthPixel   	= targetWidth * pixelPerDegree;
targWindowScale     = 3;
targWindow          = targWindowScale * targetWidth;
targWindowPixel 	= targWindow * pixelPerDegree;
distractorWidth     = targetWidth*.5;
distractorWidthPixel = distractorWidth * pixelPerDegree;
distWindow    = targWindow;
distWindowPixel = targWindowPixel;
% distAngle   = targAngle + 180;
leftTargAngle    	= rightTargAngle + 180;
targetColor         = [200 200 200];
distractorColor     = targetColor;

% Checkered Stimulus constants
nCheckerColumn     = 10;
nCheckerRow        = 10;
nSquare             = nCheckerColumn * nCheckerRow;
squareWidthPixel    = 3; % pixels
squareWidth = squareWidthPixel / pixelPerDegree;



% Visual vs auditory stop signal
switch stopMode
    case 'visual'
        % Go signal constants
        goSignalWidthPixel  = fixationWidthPixel - 2;
        goSignalColor       = [50 50 50];
        
        % Stop signal constants
        stopSignalWidthPixel = fixationWidthPixel - 2;
        stopSignalColor     = [255 150 0];
    case 'auditory'
        % Go signal constants
        goSignalWidthPixel  = fixationWidthPixel - 2;
        goSignalColor       = fixationColor;
        goSignalColor       = [50 50 50];
        
        % Load the audio into the buffer
        stopHz = 500;
        % stopHz = 750;
        % stopHz = 1000;
        stopWavFile = [num2str(stopHz), 'Hz.wav'];
        
        try
            setupAudio;
        catch err
            disp('setupAudio function elErr')
            rethrow(err)
        end
        
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             BEGIN TASK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% ***********************************************************************
%       INITIALIZE STUFF
% ***********************************************************************
% Variables that will get filled and appended each trial
trialOutcome            = {};
trialOnset              = [];
trialDuration           = [];
abortTime               = [];
fixOn                   = [];
fixWindowEntered        = [];
fixDuration             = [];
targOn                  = [];
targDuration            = [];
distOn                  = [];
distDuration            = [];
checkerOn               = [];
checkerDuration         = [];
stopSignalOn            = [];
stopDuration            = [];
ssd                     = [];
realSSD                 = []; % calculated from Screen Flips
responseCueOn           = [];
responseOnset           = [];
targWindowEntered       = [];
feedbackOnset           = [];
feedbackDuration        = [];

targAngle               = [];
distAngle               = [];

targ1CheckerProp        = [];
checkerColor            = {};
checkerArray            = {};

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

if strcmp(stopMode, 'visual')
    % ------------------------------
    %     STOP SIGNAL STIMULUS
    % ------------------------------
    stopSignalLeft      = matlabCenterX - stopSignalWidthPixel/2;
    stopSignalTop       = matlabCenterY - stopSignalWidthPixel/2;
    stopSignalRight     = matlabCenterX + stopSignalWidthPixel/2;
    stopSignalBottom    = matlabCenterY + stopSignalWidthPixel/2;
    stopSignalLocation  = [stopSignalLeft stopSignalTop stopSignalRight stopSignalBottom];
end
% ------------------------------
%     TARGET STIMULUI
% ------------------------------
rightTargetEyeLinkX      = targetAmplitudePixel * cosd(rightTargAngle);
rightTargetEyeLinkY      = targetAmplitudePixel * sind(rightTargAngle);
rightTargetMatlabX      = matlabCenterX + rightTargetEyeLinkX;
rightTargetMatlabY      = matlabCenterY - rightTargetEyeLinkY;

rightTargetLeft         = rightTargetMatlabX - targetWidthPixel/2;
rightTargetTop          = rightTargetMatlabY - targetWidthPixel/2;
rightTargetRight        = rightTargetMatlabX + targetWidthPixel/2;
rightTargetBottom       = rightTargetMatlabY + targetWidthPixel/2;
rightTargetSquare     = [rightTargetLeft, rightTargetTop, rightTargetRight, rightTargetBottom];

leftTargetEyeLinkX      = targetAmplitudePixel * cosd(leftTargAngle);
leftTargetEyeLinkY      = targetAmplitudePixel * sind(leftTargAngle);
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
stimulusEyeLinkX     = stimulusAmplitudePixel * cosd(checkerAngle);
stimulusEyeLinkY     = stimulusAmplitudePixel * sind(checkerAngle);
stimulusMatlabX     = matlabCenterX + stimulusEyeLinkX;
stimulusMatlabY     = matlabCenterY - stimulusEyeLinkY;
stimulusWindowPixel      = squareWidthPixel * nCheckerColumn; % For now, make it window the same size as the stimulus
checkerWindow = stimulusWindowPixel / pixelPerDegree;
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
% DrawFormattedText(window, 'Press the space to begin block 1.', 'center', 'center', [0 200 0]);
% Screen('Flip', window);
% junk=NaN;
% beforeStart = GetSecs;
% while isnan(junk) && GetSecs - beforeStart < 30;
%     [kd, sec, kc] = KbCheck;
%     if (kd==1) && (kc(32)==1)
%         junk = 1;
%     elseif (kd==1) && (kc(81)==1)
%         clear junk
%     end
% end
% Screen('FillRect', window, backGround);
% [vbl SOT] = Screen('Flip', window);



runningTask = 1;
iTrial = 1;
taskStartTime = GetSecs;
newTrialVariables = 1; % A flag, set to zero if an abort occurs and we want all trial variables to stay the same
iBlock = 1;

% for iTrial = 0 : 2
while runningTask
    
    
    
    
    
    
    
    % _________________________________________________________________
    %       SET UP TRIAL VARIABLES
    % _________________________________________________________________
    
    
    iTrialOutcome          	= nan;
    % Initialize variables that may or may not get filled
    iTrialDuration       	= nan;
    iabortTime              = nan;
    fixationAddIndex        = randperm(length(fixationHoldAdd));
    iPreTargetFixDuration   = fixationHoldAdd(fixationAddIndex(1));
    iPostTargetFixDuration  = fixationHoldBase + fixationHoldAdd(fixationAddIndex(2));
    iFixOn                  = nan;
    iFixWindowEntered       = nan;
    iFixDuration            = nan;
    iTargOn                 = nan;
    itargDuration           = nan;
    iDistOn                 = nan;
    iDistDuration           = nan;
    iCheckerOn              = nan;
    iCheckerDuration        = nan;
    iStopSignalOn          	= nan;
    iStopDuration           = nan;
    iSSD                    = nan;
    iRealSSD                = nan;
    iResponseCueOnset       = nan;
    iResponseOnset          = nan;
    iTargWindowEntered      = nan;
    iFeedbackOnset          = nan;
    iFeedbackDuration       = nan;
    
    
    
    
    
    
    
    
    % If the previous trial was aborted before it began (user did not
    % obtain fix), then keep all the variables the same. Only update
    % them if it's a new trial
    if newTrialVariables
        
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
        
        
        
        
%         cyanGun = 174;
        cyanGun = 210;
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
         
        
        if rightTargetTrial
            iTargetAngle        = rightTargAngle;
            iDistractorAngle   	= leftTargAngle;
            iTargetSquare       = rightTargetSquare;
            iDistractorSquare   = leftTargetSquare;
            iTargetLocation     = [rightTargetMatlabX - matlabCenterX, matlabCenterY - rightTargetMatlabY];
            iDistractorLocation = [leftTargetMatlabX - matlabCenterX, matlabCenterY - leftTargetMatlabY];
            iTargWindow       = [rightTargetMatlabX - targWindowPixel, rightTargetMatlabY - targWindowPixel, rightTargetMatlabX + targWindowPixel, rightTargetMatlabY + targWindowPixel];
            iDistWindow   = [leftTargetMatlabX - distWindowPixel, leftTargetMatlabY - distWindowPixel, leftTargetMatlabX + distWindowPixel, leftTargetMatlabY + distWindowPixel];
        else
            iTargetAngle        = leftTargAngle;
            iDistractorAngle   	= rightTargAngle;
            iTargetSquare       = leftTargetSquare;
            iDistractorSquare   = rightTargetSquare;
            iTargetLocation     = [leftTargetMatlabX - matlabCenterX, matlabCenterY - leftTargetMatlabY];
            iDistractorLocation = [rightTargetMatlabX - matlabCenterX, matlabCenterY - rightTargetMatlabY];
            iTargWindow       = [leftTargetMatlabX - targWindowPixel, leftTargetMatlabY - targWindowPixel, leftTargetMatlabX + targWindowPixel, leftTargetMatlabY + targWindowPixel];
            iDistWindow   = [rightTargetMatlabX - distWindowPixel, rightTargetMatlabY - distWindowPixel, rightTargetMatlabX + distWindowPixel, rightTargetMatlabY + distWindowPixel];
        end
    end
    
    
    %    GO OR STOP TRIAL?
    % -------------------------------
    randomProportion = rand;
    if randomProportion > pStopTrial
        iTrialType = 'go';
    else
        iTrialType = 'stop';
    end
    
    
    
    
    % Initialize stage logical variables each trial
    addTrialDataFlag    = 1;  % Gets set to false for pre-fixation aborts
    stageStartTrial    = 1;
    stageFixation       = 0;
    stageTargOn         = 0;
    stageChoiceOn       = 0;
    stageStopSignalOn   = 0;
    stageInFlight       = 0;
    stageKeyPress       = 0;
    stageOnTarget       = 0;
    stageOnDistractor   = 0;
    stageFeedback       = 0;
    
    preStopSignalResponse = 0;
    
    
    
    
    
    
    
    
    
    
    
    % ********************************************************************************************************************************
    %                                           TRIAL STAGES
    % ********************************************************************************************************************************
    % The various stages, specific for visual vs. auditory stop signals and
    % for saccadic vs. keypress responses
    
    
    
    
    
    
    
    if stageStartTrial
        % *********************************************
        %            START TRIAL STAGE
        % *********************************************
        % Turn on the fixation spot and wait forunit subject to start trial by
        % fixating and pressing space bar (which initiates a drift correction
        % and subsequently the trial
        
        
        Screen('FillRect', window, backGround);
        Screen('FillRect', window, fixationColor, fixationSquare);
        %     Screen('FrameRect', window, targetColor, fixWindowPixel);
        %         DrawFormattedText(window, 'Press the space to drift correct.', 'center', 'center', [0 200 0]);
        [~, fixWindowEnteredTime, ~, ~, ~] = Screen('Flip', window);
        
        
        
        
        iTrialOnset             = GetSecs - taskStartTime;
        iTrialOnsetComputerTime = GetSecs;
        
        
        
        iFixOn = fixWindowEnteredTime - iTrialOnsetComputerTime;
        tic
        while stageStartTrial
            
            
            
            if toc > graceObtainFix
                iFixWindowEntered      = GetSecs - iTrialOnsetComputerTime;
                stageStartTrial    = 0;
                stageFixation       = 1;
            end
            % We wait 1 ms each loop-iteration so that we
            % don't overload the system in realtime-priority:
            WaitSecs(0.001);
        end
    end
    
    
    
    
    
    
    
    
    if stageFixation
        % ****************************************************************************************
        %            FIXATION STAGE
        % ****************************************************************************************
        % Subject has begun fixating, can either hold fixation until the
        % targets come on and advance to stageTargOn, or can abort the trial
        % and start a new trial by not maintaining fixation
        
        % Prepare the screen for target onsets
        Screen('FillRect', window, fixationColor, fixationSquare);
        Screen('FillRect', window, targetColor, rightTargetSquare);
        Screen('FillRect', window, targetColor, leftTargetSquare);
        %                 Screen('FrameRect', window, targetColor, iTargWindow);
        %                 Screen('FrameRect', window, targetColor, idistWindow);
        %                 Screen('FrameRect', window, targetColor, fixWindowPixel);
        tic
        while stageFixation
            
            if toc >= iPreTargetFixDuration
                stageFixation       = 0;
                stageTargOn       = 1;
            end
        end
    end
    
    
    
    
    
    if stageTargOn
        % ****************************************************************************************
        %            PERIPHERAL TARGETS APPEAR
        % ****************************************************************************************
        % Subject has obtained fixation, and the targets appeared,
        % can either conitue to fixate until choice stimulus comes on, or can abort the trial
        % and start a new trial by not maintaining fixation
        [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
        iTargOn            = StimulusOnsetTime - iTrialOnsetComputerTime;
        iDistOn        = StimulusOnsetTime - iTrialOnsetComputerTime;
        
        % Prepare the screen for choice stimulus onset
        Screen('FillRect', window, targetColor, rightTargetSquare);
        Screen('FillRect', window, targetColor, leftTargetSquare);
        Screen('FillRect', window, iStimulusColorsArray, stimulusSquaresArray);
%         Screen('FillRect', window, fixationColor, fixationSquare);
%         Screen('FillRect', window, goSignalColor, goSignalSquare);
        incorrectTextColor = [250, 50, 50];
        tt = sprintf('%.2f \t %d \t %d \t %d \t %.2f', iRandomRightP, sum(randomRight), 100-sum(iCheckerboardArray), nMinority, iMinorityP);
            DrawFormattedText(window, tt, 'center', 'center', incorrectTextColor);
        %        Screen('FrameRect', window, targetColor, iTargWindow);
        %         Screen('FrameRect', window, targetColor, idistWindow);
        %         Screen('FrameRect', window, targetColor, fixWindowPixel);
        tic
        while stageTargOn
            
            
            if  toc >= iPostTargetFixDuration
                stageTargOn       = 0;
                stageChoiceOn       = 1;
            end
        end % while stageTargOn
    end
    
    
    
    
    
    
    
    if stageChoiceOn
        decideTime = .5 + rand(1);
        % ****************************************************************************************
        %            CHOICE STIMULUS APPEARS
        % ****************************************************************************************
        [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
        iFixDuration   = StimulusOnsetTime - iTrialOnsetComputerTime - iFixOn;
        iCheckerOn    = StimulusOnsetTime - iTrialOnsetComputerTime;
        iResponseCueOnset       = StimulusOnsetTime - iTrialOnsetComputerTime;
        
        preStopSignalSaccade = 0;
        tic
        switch iTrialType
            case 'go'
                while stageChoiceOn
                    
                    if toc > decideTime
                        
                        switch effector
                            case 'saccade'
                                iResponseOnset = GetSecs - iTrialOnsetComputerTime;
                                stageChoiceOn         = 0;
                                stageInFlight       = 1;
                            case 'keypress'
                                
                                iResponseOnset = GetSecs - iTrialOnsetComputerTime;
                                stageChoiceOn        = 0;
                                stageKeyPress       = 1;
                        end
                    end
                    
                end
            case 'stop'
                % Fill video buffer with stop signal
                if strcmp(stopMode, 'visual')
                    Screen('FillRect', window, targetColor, rightTargetSquare);
                    Screen('FillRect', window, targetColor, leftTargetSquare);
                    Screen('FillRect', window, iStimulusColorsArray, stimulusSquaresArray);
                    Screen('FillRect', window, fixationColor, fixationSquare);
                    Screen('FillRect', window, stopSignalColor, stopSignalLocation);
                end
                iStopSignalOnFlag = 0;
                iSSDIndex   = staircase(lastStopOutcomeArray{iPropRightTargIndex}, iSSDIndexArray(iPropRightTargIndex), nSSD);
                iSSD        = SSDArray(iSSDIndex)
                iSSDIndexArray(iPropRightTargIndex) = iSSDIndex;   % Update the array to staircase SSDs in each discriminability level
                
                
                while stageChoiceOn
                    
                    switch effector
                        case 'saccade'
                            
                            % turn on the Stop Signal when appropriate
                            if ~iStopSignalOnFlag && toc >= iSSD
                                switch stopMode
                                    case 'visual'
                                        % flash the stop signal
                                        [~, stopSignalOnTime, ~, ~, ~] = Screen('Flip', window);
                                    case 'auditory'
                                        % sound the stop signal
                                        % Start audio playback for 'repetitions' repetitions of the sound data,
                                        % start it immediately (0) and wait for the playback to start, return onset
                                        % timestamp.
                                        stopSignalOnTime = PsychPortAudio('Start', pahandle, repetitions, 0, 1);
                                        
                                end
                                stopSignalOnFlag = 1;
                                iRealSSD    = stopSignalOnTime - iTrialOnsetComputerTime - iCheckerOn;
                                iStopSignalOn  = stopSignalOnTime - iTrialOnsetComputerTime
                                stageChoiceOn = 0;
                                stageStopSignalOn = 1;
                                
                            end
                        case 'keypress'
                            
                            % turn on the Stop Signal when appropriate
                            if ~iStopSignalOnFlag && toc >= iSSD
                                switch stopMode
                                    case 'visual'
                                        % flash the stop signal
                                        [~, stopSignalOnTime, ~, ~, ~] = Screen('Flip', window);
                                    case 'auditory'
                                        % sound the stop signal
                                        % Start audio playback for 'repetitions' repetitions of the sound data,
                                        % start it immediately (0) and wait for the playback to start, return onset
                                        % timestamp.
                                        stopSignalOnTime = PsychPortAudio('Start', pahandle, repetitions, 0, 1);
                                end
                                iStopSignalOnFlag = 1;
                                iRealSSD    = stopSignalOnTime - iTrialOnsetComputerTime - iCheckerOn;
                                iStopSignalOn  = stopSignalOnTime - iTrialOnsetComputerTime
                                stageChoiceOn = 0;
                                stageStopSignalOn = 1;
                                
                            end
                    end %swithc effector
                end
                
        end % stage
    end
    
    
    
    
    
    
    
    
    
    
    
    if stageStopSignalOn
        while stageStopSignalOn
            % ****************************************************************************************
            %       IF IT'S A STOP TRIAL, THE STOP SIGNAL OCCURED
            % ****************************************************************************************
            
            WaitSecs(decideTime)
            
            saccadeOrStop = randi(2) - 1;
            switch effector
                case 'saccade'
                    % Made a saccade
                    if saccadeOrStop == 0
                        iResponseOnset = GetSecs - iTrialOnsetComputerTime;
                        stageStopSignalOn         = 0;
                        stageInFlight       = 1;
                    elseif saccadeOrStop == 1
                        iTrialOutcome       = 'stopCorrect';
                        lastStopOutcomeArray{iPropRightTargIndex} = iTrialOutcome;
                        stageStopSignalOn         = 0;
                        stageFeedback       = 1;
                    end
                case 'keypress'
                    
                    if saccadeOrStop == 1
                        iTrialOutcome       = 'stopCorrect';
                        lastStopOutcomeArray{iPropRightTargIndex} = iTrialOutcome;
                        stageStopSignalOn         = 0;
                        stageFeedback       = 1;
                    elseif saccadeOrStop == 0
                        iResponseOnset = GetSecs - iTrialOnsetComputerTime;
                        stageStopSignalOn         = 0;
                        stageKeyPress       = 1;
                    end
            end % swithc effector
        end % while stageStopSignalOn
        
        % Indicate that we've come far enough along to present new stimuli next
        % time
        newTrialVariables = 1;
    end
    
    
    
    
    
    if stageKeyPress
        % ****************************************************************************************
        %       IF IT'S A KEY PRESS VERSION OF TASK, PROCESS THE KEY PRESS
        % ****************************************************************************************
        correctTarg = randi(2) - 1;
        switch iTrialType
            case 'go'
                if correctTarg
                    % Flag it was correct right
                    iTrialOutcome = 'goCorrectTarget';
                else
                    % Flag as distractor left
                    iTrialOutcome = 'goCorrectDistractor';
                end
                stageFeedback = 1;
            case 'stop'
                if correctTarg
                    % Flag it incorrected stop to correct right
                    iTrialOutcome = 'stopIncorrectTarget';
                else
                    % Flag it incorrected stop to distractor left
                    iTrialOutcome = 'stopIncorrectDistractor';
                end
                lastStopOutcomeArray{iPropRightTargIndex} = iTrialOutcome;
                stageFeedback = 1;
        end
    end
    
    
    
    
    
    
    
    if stageInFlight
        % **********************************************************************************************
        %         IF EYES HAVE MOVED, SEE WHERE THEY GO
        % **********************************************************************************************
        correctTarg = randi(2) - 1;
        tic
        while stageInFlight
            
            
            switch effector
                case 'saccade'
                    if toc > graceSaccDuration && correctTarg
                        iTargWindowEntered = GetSecs - iTrialOnsetComputerTime;
                        stageInFlight       = 0;
                        stageOnTarget       = 1;
                    elseif toc > graceSaccDuration && ~correctTarg
                        iTargWindowEntered = GetSecs - iTrialOnsetComputerTime;
                        stageInFlight       = 0;
                        stageOnDistractor       = 1;
                    end
                    
            end % switch effector
        end
    end
    
    
    
    
    
    
    if stageOnTarget
        % **********************************************************************************************
        %        EYES HAVE LANDED ON TARGET
        % **********************************************************************************************
        tic
        while stageOnTarget
            
            
            
            
            switch iTrialType
                case 'go'
                    if toc >= postSaccHoldDuration
                        % Flag it was correct to the target
                        iTrialOutcome = 'goCorrectTarget';
                        stageOnTarget = 0;
                        stageFeedback = 1;
                    end
                case 'stop'
                    if toc >= postSaccHoldDuration
                        % Flag it was correct to the target
                        iTrialOutcome = 'stopIncorrectTarget';
                        stageOnTarget = 0;
                        stageFeedback = 1;
                        lastStopOutcomeArray{iPropRightTargIndex} = iTrialOutcome;
                    end
            end
        end
        
    end
    
    
    
    
    
    
    
    if stageOnDistractor
        % **********************************************************************************************
        %         EYES HAVE LANDED ON DISTRACTOR
        % **********************************************************************************************
        tic
        while stageOnDistractor
            
            
            
            switch iTrialType
                case 'go'
                    if toc >= postSaccHoldDuration
                        % Flag it was correct to the target
                        iTrialOutcome = 'goCorrectDistractor';
                        stageOnDistractor = 0;
                        stageFeedback = 1;
                    end
                case 'stop'
                    if toc >= postSaccHoldDuration
                        % Flag it was correct to the target
                        iTrialOutcome = 'stopIncorrectDistractor';
                        stageOnDistractor = 0;
                        stageFeedback = 1;
                        lastStopOutcomeArray{iPropRightTargIndex} = iTrialOutcome;
                    end
            end
        end
        
    end
    
    
    
    
    
    
    
    if stageFeedback
        % *****************************************************************
        %        DELIVER FEEDBACK
        % *****************************************************************
        while stageFeedback
            if iTrial > 20
                feedbackTime = .1;
            end
            [window, timeout] = feedback(window, iTrialOutcome, iTrialType);
            [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
            itargDuration         = StimulusOnsetTime - iTrialOnsetComputerTime + iTargOn;
            iDistDuration     = StimulusOnsetTime - iTrialOnsetComputerTime + iDistOn;
            iCheckerDuration = StimulusOnsetTime - iTrialOnsetComputerTime + iCheckerOn;
            iFeedbackOnset          = StimulusOnsetTime - iTrialOnsetComputerTime;
            iFeedbackDuration       = feedbackTime + timeout;
            if strcmp(iTrialType, 'stop')
                iStopDuration = StimulusOnsetTime - iTrialOnsetComputerTime + iStopSignalOn;
            end
            
            waitSecs(feedbackTime + timeout);
            stageFeedback = 0;
        end
        
        iTrialDuration = GetSecs - iTrialOnsetComputerTime;
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % Allow option to quit task by pressing escape at the end of a trial
    tic
    while toc < 1
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        if keyIsDown && find(keyCode) == escapeKey
            runningTask = 0;
        end
    end
    
    
    
end % trial loop




cleanup

    function cleanup
        Screen('CloseAll');
        ShowCursor;
    end






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
        elseif strcmp(iTrialType, 'go') && (strcmp(iTrialOutcome, 'targetHoldAbort') || strcmp(iTrialOutcome, 'distractorHoldAbort'))
            DrawFormattedText(window, 'Keep your eyes on the target until it disappears ', 'center', 'center', incorrectTextColor);
            timeout = .2;
        elseif strcmp(iTrialOutcome, 'eyelinkError')
            DrawFormattedText(window, 'My fault- continue on', 'center', 'center', incorrectTextColor);
            timeout = 0;
        end
        
        % Screen('Flip', window);
        % waitSecs(feedbackTime + timeout)
    end






% *******************************************************************************
    function setupAudio
        % Adapted from BasicSountOutputDemo.m and BasicSoundScheduleDemo.m in
        % psychtoolbox
        
        % Always init to 2 channels, for the sake of simplicity:
        nrchannels = 2;
        
        % Does a function for resampling exist?
        if exist('resample') %#ok<EXIST>
            % Yes: Select a target sampling rate of 44100 Hz, resample if
            % neccessary:
            freq = 44100;
            doresample = 1;
        else
            % No. We will choose the frequency of the wav file with the highest
            % frequency for actual playback. Wav files with deviating frequencies
            % will play too fast or too slow, b'cause we can't resample:
            % Init freq:
            freq = 0;
            doresample = 0;
        end
        
        % Perform basic initialization of the sound driver:
        InitializePsychSound(1);
        
        
        [audiodata, infreq] = wavread(char(stopWavFile));
        
        if doresample
            % Resampling supported. Check if needed:
            if infreq ~= freq
                % Need to resample this to target frequency 'freq':
                fprintf('Resampling from %i Hz to %i Hz... ', infreq, freq);
                audiodata = resample(audiodata, freq, infreq);
            end
        else
            % Resampling not supported by Matlab/Octave version:
            % Adapt final playout frequency to maximum frequency found, and
            % hope that all files match...
            freq = max(infreq, freq);
        end
        
        [samplecount, ninchannels] = size(audiodata);
        audiodata = repmat(transpose(audiodata), nrchannels / ninchannels, 1);
        
        buffer = PsychPortAudio('CreateBuffer', [], audiodata); %#ok<AGROW>
        [fpath, fname] = fileparts(char(stopWavFile));
        fprintf('Filling audiobuffer handle with soundfile %s ...\n', fname);
        
        
        % Open the default audio device [], with default mode [] (==Only playback),
        % and a required latencyclass of 1 == standard low-latency mode, as well as
        % a playback frequency of 'freq' and 'nrchannels' sound output channels.
        % This returns a handle 'pahandle' to the audio device:
        pahandle = PsychPortAudio('Open', [], [], 1, freq, nrchannels);
        
        % For the fun of demoing this as well, we switch PsychPortAudio to runMode
        % 1, instead of the default runMode 0. This will slightly increase the cpu
        % load and general system load, but provide better timing and even lower
        % sound onset latencies under certain conditions. It is not really needed
        % in this demo, just here to grab your attention for this feature. Type
        % PsychPortAudio RunMode? for more details...
        runMode = 1;
        PsychPortAudio('RunMode', pahandle, runMode);
        
        
        % Fill the audio playback buffer with the audio data 'wavedata':
        PsychPortAudio('FillBuffer', pahandle, audiodata);
        repetitions = 1;
    end %function setupAudio


end  % main function


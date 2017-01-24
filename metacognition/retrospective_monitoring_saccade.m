function retrospective_monitoring_saccade(soaArrayScreenFlips, nTarget, targetAmplitude, targetAngles, betAngle, targetDistribution, plotFlag, saveFlag)

% stimulusAmplitude: in some units (e.g. pixels), the distance from the
%       center of the screen to the center of the stimulus
% stimulusAngle: in degrees, the angle from the center of the screen to the
%       center of the stimulus


% example soaArrayScreenFlips = [22 26 30 34]




if nargin == 0
    soaArrayScreenFlips = [2 4 6 8];
    targetAmplitude = 10;
    maskAmplitude = targetAmplitude;
    targetAngles = 45;
    betAngle = 45;
    nTarget = 4;
    targetDistribution = 'maximize';
    plotFlag = 1;
    saveFlag = 1;
end

targetAngleArray    = get_mask_angle_array(targetAngles, nTarget, targetDistribution);
betAngleArray       = get_bet_angle_array(betAngle, targetDistribution);
% targetAngleArray = [45 135 225 315];
% betAngleArray = [45 135];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Set up Experiment Variables, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Keyboard assignments
KbName('UnifyKeyNames');
escapeKey               = KbName('ESCAPE');


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
maskAmplitudePixel     = targetAmplitudePixel;
betAmplitudePixel = targetAmplitude * pixelPerDegree;

% Get info for saving a data file.
subjectID     = input('Enter subject            ', 's');
clockVector     = clock;
session         = [datestr(now,'mmdd'),  'saccade'];
saveFileName    = ['data/', subjectID, session];

backGround      = [60 60 60];
whichScreen     = 0;
[window, centerPoint] = Screen('OpenWindow', whichScreen, backGround);
priorityLevel   = MaxPriority(window);
Priority(priorityLevel);

commandwindow;
flipTime        = Screen('GetFlipInterval',window);% get the flip rate of current monitor.
preFlipBuffer   = flipTime / 2;
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
dummymode               = 0;
trialsPerBlock          = 5;
BLOCKS_TO_RUN           = 20;
totalTrial              = trialsPerBlock * BLOCKS_TO_RUN;
soaArray                = soaArrayScreenFlips * (flipTime); %subtracting 2 ms each cycle to enusre the SSD occurs before the next screen refresh
graceResponse           = 1.7; % seconds allowed to make a saccade
graceObtainFixation     = 2;  % seconds program will wait for user to obtain fixation
graceSaccadeDuration    = .2; % seconds allowed intra saccade time
postSaccadeHoldDuration = .4;  % duration to hold post-saccade fixation
stimulusScaleConversion = 1 / 10; % for now 10 pixels per every 100 pixels from fixation
feedbackTime            = .4;


% Fixation Spot constants
% Fixation Spot constants
fixationAngle       = 0;
fixationAmplitude   = 0;
fixationWidth       = .5;
fixationWidthPixel  = fixationWidth * pixelPerDegree;
fixWindowScale      = 3;  % fix window is 4 times size of fix point
fixationWindow      = fixationWidth*fixWindowScale;
fixationWindowPixel	= [-fixationWidthPixel*fixWindowScale, -fixationWidthPixel*fixWindowScale, fixationWidthPixel*fixWindowScale, fixationWidthPixel*fixWindowScale];
fixationWindowPixel      = CenterRect(fixationWindowPixel, centerPoint);
decisionFixColor       = [200 200 200];
betFixColor       = [0 200 0];
fixationHoldBase    = .5;
fixationHoldAdd     = (0 : 10 : 500)  / 1000;



targetWidth         = 1;
targetWidthPixel   	= targetWidth * pixelPerDegree;
targetWindowScale   = 4;
targetWindow        = targetWindowScale * targetWidth;
targetWindowPixel 	= targetWindow * pixelPerDegree;
targetColor         = [50 0 0];

maskWidth         = 1.2;
maskWidthPixel   	= maskWidth * pixelPerDegree;
maskWindow          = targetWindow;
maskWindowPixel 	= maskWindow * pixelPerDegree;
maskColor         = [200 200 200];

betWidth         = 1;
betWidthPixel   	= betWidth * pixelPerDegree;
betWindow       = targetWindow;
betWindowPixel 	= betWindow * pixelPerDegree;
cyanGun = 175;
magentaGun = 255;
highBetColor = [0 cyanGun cyanGun];
lowBetColor = [magentaGun 0 magentaGun];








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
trialOnset              = [];
trialDuration           = [];

decisionFixSpotOnset  	= [];
decisionFixationOnset  	= [];
decisionFixSpotDuration = [];
decisionTargOnset     	= [];
soa                     = [];
decisionMaskOnset   	= [];
decisionResponseCueOnset = [];
decisionResponseOnset   = [];
targetWindowEntered     = [];
targetIndex   	= [];

betFixSpotOnset         = [];
betFixationOnset      	= [];
preBetFixDuration       = [];
betTargetOnset          = [];
betResponseOnset        = [];
betWindowEntered     = [];
betIndex   	= [];

feedbackOnset           = [];
feedbackDuration        = [];
abortOnset              = [];
trialOutcome            = {};

targetAngle             = [];
highBetAngle            = [];
lowBetAngle          = [];




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

maskSquare = nan(nTarget, 4);
betSquare = nan(nTarget, 4);
targetSquare = nan(nTarget, 4);

% ------------------------------
%     MASK STIMULUI
% ------------------------------
for iPosition = 1 : nTarget
    targetEyeLinkX(iPosition) = targetAmplitudePixel * cosd(targetAngleArray(iPosition));
    targetEyeLinkY(iPosition) = targetAmplitudePixel * sind(targetAngleArray(iPosition));
    targetMatlabX(iPosition)      = matlabCenterX + targetEyeLinkX(iPosition);
    targetMatlabY(iPosition)      = matlabCenterY - targetEyeLinkY(iPosition);
    
    targetLeft         = targetMatlabX(iPosition) - targetWidthPixel/2;
    targetTop          = targetMatlabY(iPosition) - targetWidthPixel/2;
    targetRight       = targetMatlabX(iPosition) + targetWidthPixel/2;
    targetBottom      = targetMatlabY(iPosition) + targetWidthPixel/2;
    targetSquare(iPosition,:)     = [targetLeft, targetTop, targetRight, targetBottom];
    
    maskLeft         = targetMatlabX(iPosition) - maskWidthPixel/2;
    maskTop          = targetMatlabY(iPosition) - maskWidthPixel/2;
    maskRight       = targetMatlabX(iPosition) + maskWidthPixel/2;
    maskBottom      = targetMatlabY(iPosition) + maskWidthPixel/2;
    maskSquare(iPosition,:)       = [targetLeft, targetTop, targetRight, targetBottom];
    
    maskWindowCoord(iPosition,:)       = [targetMatlabX(iPosition) - maskWindowPixel/2, targetMatlabY(iPosition) - maskWindowPixel/2, targetMatlabX(iPosition) + maskWindowPixel/2, targetMatlabY(iPosition) + maskWindowPixel/2];
end

% ------------------------------
%     BET STIMULUI
% ------------------------------
for iPosition = 1 : 2
    betEyeLinkX(iPosition) = targetAmplitudePixel * cosd(betAngleArray(iPosition));
    betEyeLinkY(iPosition) = targetAmplitudePixel * sind(betAngleArray(iPosition));
    betMatlabX(iPosition)      = matlabCenterX + betEyeLinkX(iPosition);
    betMatlabY(iPosition)      = matlabCenterY - betEyeLinkY(iPosition);
    
    betLeft         = betMatlabX(iPosition) - betWidthPixel/2;
    betTop          = betMatlabY(iPosition) - betWidthPixel/2;
    betRight        = betMatlabX(iPosition) + betWidthPixel/2;
    betBottom       = betMatlabY(iPosition) + betWidthPixel/2;
    betSquare(iPosition,:)     = [betLeft, betTop, betRight, betBottom];
    
    betWindowCoord(iPosition,:)       = [betMatlabX(iPosition) - betWindowPixel/2, betMatlabY(iPosition) - betWindowPixel/2, betMatlabX(iPosition) + betWindowPixel/2, betMatlabY(iPosition) + betWindowPixel/2];
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
    
    
    
    % Initialize variables that may or may not get filled
    % Decision Stage
    fixationAddIndex        = randperm(length(fixationHoldAdd));
    iPreTargetFixDuration   = fixationHoldBase + fixationHoldAdd(fixationAddIndex(1));
    iPostMaskFixDuration  = fixationHoldBase + fixationHoldAdd(fixationAddIndex(2))*2;
    iDecisionFixSpotOnset   = nan;
    iDecisionFixation       = nan;
    iDecisionFixSpotDuration   = nan;
    iDecisionTargOnset    	= nan;
    iRealSOA                = nan;
    iSOA                    = nan;
    iDecisionMaskOnset     	= nan;
    iDecisionResponseCueOnset       = nan;
    iDecisionResponseOnset          = nan;
    iTargetIndex          = nan;
    iDecisionIndex          = nan;
    
    % Bet Stage
    iBetFixSpotOnset            = nan;
    iBetFixation            = nan;
    iPreBetFixDuration      = fixationHoldBase + fixationHoldAdd(fixationAddIndex(3));
    iBetTargetOnset         = nan;
    iBetResponseOnset       = nan;
    
    % Outcome
    iTrialOutcome           = nan;
    iAbortOnset             = nan;
    iFeedbackOnset          = nan;
    iFeedbackDuration       = nan;
    
    
    
    
    
    
    
    
    % If the previous trial was aborted before it began (user did not
    % obtain fixation), then keep all the variables the same. Only update
    % them if it's a new trial
    if newTrialVariables
        
        %   WHERE IS THE DECISION TARGET?
        % -------------------------------
        iTarget = randi(nTarget);
        iTargetIndex = iTarget;
        
        iTargetSquare = targetSquare(iTarget, :);
        iTargetWindow     = [targetMatlabX(iTarget) - matlabCenterX, matlabCenterY - targetMatlabY(iTarget)];
        
        %   WHERE ARE THE HIGH AND LOW BET TARGETS?
        % -------------------------------
        iBetOrder = randperm(2);
        iHighBet = iBetOrder(1);
        iLowBet = iBetOrder(2);
        
        iHighBetSquare = betSquare(iHighBet,:);
        iLowBetSquare = betSquare(iLowBet,:);
        iHighBetLocation     = [betMatlabX(iHighBet) - matlabCenterX, matlabCenterY - betMatlabY(iHighBet)];
        iLowBetLocation     = [betMatlabX(iLowBet) - matlabCenterX, matlabCenterY - betMatlabY(iLowBet)];
    end
    
    
    
    
    % ****************************************************************************************
    %            BEGIN STAGES
    % ****************************************************************************************
    % Initialize stage logical variables each trial
    addTrialDataFlag    =  1;  % Gets set to true for pre-fixation aborts
    abortTrialFlag              = 0;
    
    stageDecisionPreFix         = 1;
    stageDecisionFix            = 0;
    stageDecisionTargetOn       = 0;
    stageDecisionMaskOn         = 0;
    stageCueToDecide            = 0;
    stageDecisionInFlight       = 0;
    stageDecisionOnTarget       = 0;
    
    stageBetPreFix              = 0;
    stageBetPostFix             = 0;
    stageBetTargetOn            = 0;
    stageBetInFlight            = 0;
    stageBetOnTarget            = 0;
    
    stageFeedback               = 0;
    
    
    
    % ****************************************************************************************
    %            PRE-FIXATION STAGE
    % ****************************************************************************************
    % Turn on the fixation spot and wait unit subject fixates the spot- or,
    % if subject does not fixate spot within a grace period, start a new
    % trial
    
    
    Screen('FillRect', window, backGround);
    Screen('FillRect', window, decisionFixColor, fixationSquare);
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
    
    
    
    
    
    
    iDecisionFixSpotOnset = fixationOnsetTime - iTrialOnsetComputerTime;
    tic
    while stageDecisionPreFix
        
        %--------------------------------------------------------------------------------------------
        % Unit code for checking eye position
        error = Eyelink('CheckRecording');
        if error
            fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
            iTrialOutcome = 'eyelinkError';
            stageFeedback = 1;
            break
        end
        if Eyelink( 'NewFloatSampleAvailable') > 0
            % get the sample in the form of an event structure
            event = Eyelink( 'NewestFloatSample');
            x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
            y = event.gy(eyeRecorded+1);
            % do we have valid data and is the pupil visible?
            if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                eyeX = x;
                eyeY = y;
            end
        end
        %----------------------------------------------------------------------------------------------
        
        
        if inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc <= graceObtainFixation
            Eyelink('Message', 'Fixation Start');
            iDecisionFixation      = GetSecs - iTrialOnsetComputerTime;
            stageDecisionPreFix = 0;
            stageDecisionFix = 1;
        elseif GetSecs > fixationOnsetTime + graceObtainFixation
            % If subject aborted, exit this stage and start a new trial
            % with all other parameters the same
            stageDecisionPreFix = 0;
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
    % target comes on and advance to stageTargetOn, or can abort the trial
    % and start a new trial by not maintaining fixation
    if stageDecisionFix
        
        % Prepare the screen for target onset
        Screen('FillRect', window, targetColor, iTargetSquare);
        Screen('FillRect', window, decisionFixColor, fixationSquare);
        tic
        while stageDecisionFix
            
            %--------------------------------------------------------------------------------------------
            % Unit code for checking eye position
            error = Eyelink('CheckRecording');
            if error
                fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
                iTrialOutcome = 'eyelinkError';
                stageFeedback = 1;
                break
            end
            if Eyelink( 'NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                event = Eyelink( 'NewestFloatSample');
                x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                y = event.gy(eyeRecorded+1);
                % do we have valid data and is the pupil visible?
                if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                    eyeX = x;
                    eyeY = y;
                end
            end
            %----------------------------------------------------------------------------------------------
            
            if ~inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc < iPreTargetFixDuration
                iTrialOutcome       = 'fixationAbort';
                iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
                abortTrialFlag = 1;
                Eyelink('Message', 'Fixation Abort');
                stageDecisionFix = 0;
                stageFeedback = 1;
                newTrialVariables            = 0;
            elseif inAcceptWindow(eyeX, eyeY, fixationWindowPixel) &&  toc >= iPreTargetFixDuration
                stageDecisionFix = 0;
                stageDecisionTargetOn = 1;
                [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
                iDecisionTargOnset            = StimulusOnsetTime - iTrialOnsetComputerTime;
            end
        end
    end
    
    
    
    % ****************************************************************************************
    %            TARGET ON STAGE
    % ****************************************************************************************
    % Subject has obtained fixation, and the targets appeared,
    % can either conitue to fixate until choice stimulus comes on, or can abort the trial
    % and start a new trial by not maintaining fixation
    if stageDecisionTargetOn
        %         Screen('FillRect', window, targetColor, iTargetSquare);
        Screen('FillRect', window, decisionFixColor, fixationSquare);
        % Load up each mask
        for iMask = 1 : nTarget
            Screen('FillRect', window, maskColor, maskSquare(iMask, :));
        end
        waitSecs(iSOA - preFlipBuffer);
        [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
        iDecisionMaskOnset = StimulusOnsetTime - iTrialOnsetComputerTime;
        stageDecisionTargetOn = 0;
        stageDecisionMaskOn = 1;
    end
    
    
    
    
    
    
    % ****************************************************************************************
    %            MASK ON STAGE
    % ****************************************************************************************
    if stageDecisionMaskOn
        % Load up each mask
        for iMask = 1 : nTarget
            Screen('FillRect', window, maskColor, maskSquare(iMask, :));
        end
        
        tic
        while stageDecisionMaskOn
            
            %--------------------------------------------------------------------------------------------
            % Unit code for checking eye position
            error = Eyelink('CheckRecording');
            if error
                fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
                iTrialOutcome = 'eyelinkError';
                stageFeedback = 1;
                break
            end
            if Eyelink( 'NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                event = Eyelink( 'NewestFloatSample');
                x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                y = event.gy(eyeRecorded+1);
                % do we have valid data and is the pupil visible?
                if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                    eyeX = x;
                    eyeY = y;
                end
            end
            %----------------------------------------------------------------------------------------------
            
            
            if ~inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc < iPostMaskFixDuration
                iTrialOutcome       = 'decisionFixationAbort';
                iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
                abortTrialFlag = 1;
                Eyelink('Message', 'Fixation Abort');
                stageDecisionMaskOn = 0;
                stageFeedback       = 1;
            elseif inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc >= iPostMaskFixDuration
                [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
                iDecisionFixSpotDuration   = StimulusOnsetTime - iTrialOnsetComputerTime - iDecisionFixSpotOnset;
                iDecisionResponseCueOnset       = StimulusOnsetTime - iTrialOnsetComputerTime;
                Eyelink('Message', 'Targets Start');
                stageDecisionMaskOn = 0;
                stageCueToDecide = 1;
            end
        end
    end
    
    
    
    
    
    if stageCueToDecide
        tic
        % **********************************************************************************************
        %         TARGETS AND STIMULI ON STAGE
        % **********************************************************************************************
        
        while stageCueToDecide
            
            %--------------------------------------------------------------------------------------------
            % Unit code for checking eye position
            error = Eyelink('CheckRecording');
            if error
                fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
                iTrialOutcome = 'eyelinkError';
                stageFeedback = 1;
                break
            end
            if Eyelink( 'NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                event = Eyelink( 'NewestFloatSample');
                x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                y = event.gy(eyeRecorded+1);
                % do we have valid data and is the pupil visible?
                if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                    eyeX = x;
                    eyeY = y;
                end
            end
            %----------------------------------------------------------------------------------------------
            
            
            % Made a saccade
            if ~inAcceptWindow(eyeX, eyeY, fixationWindowPixel)
                Eyelink('Message', 'Response Onset');
                iDecisionResponseOnset = GetSecs - iTrialOnsetComputerTime;
                stageCueToDecide = 0;
                stageDecisionInFlight = 1;
                % Waited too long to make a saccade
            elseif inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc > graceResponse
                iTrialOutcome       = 'decisionResponseAbort';
                iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
                abortTrialFlag = 1;
                Eyelink('Message', 'Response Timed Out');
                stageCueToDecide = 0;
                stageFeedback = 1;
            end
        end
    end
    
    
    
    
    
    % **********************************************************************************************
    %         IN FLIGHT STAGE
    % **********************************************************************************************
    tic
    while stageDecisionInFlight
        
        %--------------------------------------------------------------------------------------------
        % Unit code for checking eye position
        error = Eyelink('CheckRecording');
        if error
            fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
            iTrialOutcome = 'eyelinkError';
            stageFeedback = 1;
            break
        end
        if Eyelink( 'NewFloatSampleAvailable') > 0
            % get the sample in the form of an event structure
            event = Eyelink( 'NewestFloatSample');
            x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
            y = event.gy(eyeRecorded+1);
            % do we have valid data and is the pupil visible?
            if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                eyeX = x;
                eyeY = y;
            end
        end
        %----------------------------------------------------------------------------------------------
        
        
        if toc > graceSaccadeDuration
            iTrialOutcome       = 'saccadeAbort';
            iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
            abortTrialFlag = 1;
            stageDecisionInFlight = 0;
            stageFeedback = 1;
        else
            % Is the eye in any of the mask windows?
            for iPosition = 1 : nTarget
                if inAcceptWindow(eyeX, eyeY, maskWindowCoord(iPosition,:)) && toc <= graceSaccadeDuration
                    iDecisionIndex = iPosition;
                    stageDecisionInFlight = 0;
                    stageDecisionOnTarget = 1;
                    break
                end
            end
        end
    end
    
    
    
    % **********************************************************************************************
    %         ON TARGET STAGE
    % **********************************************************************************************
    Screen('FillRect', window, betFixColor, fixationSquare);
    tic
    while stageDecisionOnTarget
        
        %--------------------------------------------------------------------------------------------
        % Unit code for checking eye position
        error = Eyelink('CheckRecording');
        if error
            fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
            iTrialOutcome = 'eyelinkError';
            stageFeedback = 1;
            break
        end
        if Eyelink( 'NewFloatSampleAvailable') > 0
            % get the sample in the form of an event structure
            event = Eyelink( 'NewestFloatSample');
            x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
            y = event.gy(eyeRecorded+1);
            % do we have valid data and is the pupil visible?
            if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                eyeX = x;
                eyeY = y;
            end
        end
        %----------------------------------------------------------------------------------------------
        
        
        if inAcceptWindow(eyeX, eyeY, maskWindowCoord(iDecisionIndex,:)) && toc >= postSaccadeHoldDuration
            % Flag whether it was correct
            if iDecisionIndex == iTarget
                iDecision = 'Correct';
            else
                iDecision = 'Incorrect';
            end
            [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
            iBetFixSpotOnset = GetSecs - iTrialOnsetComputerTime;
            stageDecisionOnTarget = 0;
            stageBetPreFix = 1;
            iBetFixSpotOnset = StimulusOnsetTime - iTrialOnsetComputerTime;
        elseif ~inAcceptWindow(eyeX, eyeY, maskWindowCoord(iDecisionIndex,:)) && toc < postSaccadeHoldDuration
            % Flag it was target aborted
            iTrialOutcome = 'decisionTargetHoldAbort';
            iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
            abortTrialFlag = 1;
            stageDecisionOnTarget = 0;
            stageFeedback = 1;
        end
    end
    
    
    
    
    % **********************************************************************************************
    %         BET PRE FIXATION STAGE
    % **********************************************************************************************
    tic
    while stageBetPreFix
        
        %--------------------------------------------------------------------------------------------
        % Unit code for checking eye position
        error = Eyelink('CheckRecording');
        if error
            fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
            iTrialOutcome = 'eyelinkError';
            stageFeedback = 1;
            break
        end
        if Eyelink( 'NewFloatSampleAvailable') > 0
            % get the sample in the form of an event structure
            event = Eyelink( 'NewestFloatSample');
            x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
            y = event.gy(eyeRecorded+1);
            % do we have valid data and is the pupil visible?
            if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                eyeX = x;
                eyeY = y;
            end
        end
        %----------------------------------------------------------------------------------------------
        
        
        if inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc <= graceObtainFixation
            Eyelink('Message', 'Bet Fixation Start');
            iBetFixation      = GetSecs - iTrialOnsetComputerTime;
            stageBetPreFix = 0;
            stageBetPostFix = 1;
        elseif toc > iBetFixSpotOnset + graceObtainFixation
            % Flag it was target aborted
            iTrialOutcome = 'interstageAbort';
            iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
            abortTrialFlag = 1;
            stageBetPreFix = 0;
            stageFeedback = 1;
        end
        % We wait 1 ms each loop-iteration so that we
        % don't overload the system in realtime-priority:
        WaitSecs(0.001);
    end
    
    
    
    
    % ****************************************************************************************
    %            BET FIXATION STAGE
    % ****************************************************************************************
    % Subject has begun fixating, can either hold fixation until the
    % target comes on and advance to stageBetTargetOn, or can abort the trial
    % and start a new trial by not maintaining fixation
    if stageBetPostFix
        % Prepare the screen for bet targets onset
        Screen('FillRect', window, highBetColor, iHighBetSquare);
        Screen('FillRect', window, lowBetColor, iLowBetSquare);
        Screen('FillRect', window, betFixColor, fixationSquare);
        % Add a smidge of time for gaze to settle on fixation- else there
        % are lots of aborts
        WaitSecs(0.1);
        tic
        while stageBetPostFix
            
            %--------------------------------------------------------------------------------------------
            % Unit code for checking eye position
            error = Eyelink('CheckRecording');
            if error
                fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
                iTrialOutcome = 'eyelinkError';
                stageFeedback = 1;
                break
            end
            if Eyelink( 'NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                event = Eyelink( 'NewestFloatSample');
                x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                y = event.gy(eyeRecorded+1);
                % do we have valid data and is the pupil visible?
                if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                    eyeX = x;
                    eyeY = y;
                end
            end
            %----------------------------------------------------------------------------------------------
            
            if ~inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc < iPreBetFixDuration
                iTrialOutcome       = 'betFixationAbort';
                iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
                abortTrialFlag = 1;
                Eyelink('Message', 'Bet Fixation Abort');
                stageBetPostFix = 0;
                stageFeedback = 1;
            elseif inAcceptWindow(eyeX, eyeY, fixationWindowPixel) &&  toc >= iPreBetFixDuration
                stageBetPostFix = 0;
                stageBetTargetOn = 1;
                [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
                iBetTargetOnset            = StimulusOnsetTime - iTrialOnsetComputerTime;
            end
        end
    end
    
    
    % ****************************************************************************************
    %            BET TARGETS ON STAGE
    % ****************************************************************************************
    if stageBetTargetOn
        
        tic
        while stageBetTargetOn
            
            %--------------------------------------------------------------------------------------------
            % Unit code for checking eye position
            error = Eyelink('CheckRecording');
            if error
                fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
                iTrialOutcome = 'eyelinkError';
                stageFeedback = 1;
                break
            end
            if Eyelink( 'NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                event = Eyelink( 'NewestFloatSample');
                x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
                y = event.gy(eyeRecorded+1);
                % do we have valid data and is the pupil visible?
                if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                    eyeX = x;
                    eyeY = y;
                end
            end
            %----------------------------------------------------------------------------------------------
            
            
            if ~inAcceptWindow(eyeX, eyeY, fixationWindowPixel)
                Eyelink('Message', 'Response Onset');
                iBetResponseOnset = GetSecs - iTrialOnsetComputerTime;
                stageBetTargetOn = 0;
                stageBetInFlight = 1;
            elseif inAcceptWindow(eyeX, eyeY, fixationWindowPixel) && toc >= graceResponse
                trialOutcome = 'betResponseAbort';
                iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
                abortTrialFlag = 1;
                Eyelink('Message', 'Response Timed Out');
                stageBetTargetOn = 0;
                stageFeedback = 1;
            end
        end
    end
    
    % **********************************************************************************************
    %         IN FLIGHT STAGE
    % **********************************************************************************************
    tic
    while stageBetInFlight
        
        
        %--------------------------------------------------------------------------------------------
        % Unit code for checking eye position
        error = Eyelink('CheckRecording');
        if error
            fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
            iTrialOutcome = 'eyelinkError';
            stageFeedback = 1;
            break
        end
        if Eyelink( 'NewFloatSampleAvailable') > 0
            % get the sample in the form of an event structure
            event = Eyelink( 'NewestFloatSample');
            x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
            y = event.gy(eyeRecorded+1);
            % do we have valid data and is the pupil visible?
            if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                eyeX = x;
                eyeY = y;
            end
        end
        %----------------------------------------------------------------------------------------------
        
        
        if toc > graceSaccadeDuration
            iTrialOutcome       = 'saccadeAbort';
            iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
            abortTrialFlag = 1;
            stageBetInFlight = 0;
            stageFeedback = 1;
        else
            % Is the eye in any of the mask windows?
            for iPosition = 1 : 2
                if inAcceptWindow(eyeX, eyeY, betWindowCoord(iPosition,:)) && toc <= graceSaccadeDuration
                    iBetIndex = iPosition;
                    stageBetInFlight = 0;
                    stageBetOnTarget = 1;
                end
            end
        end
    end
    
    
    
    % **********************************************************************************************
    %         ON TARGET STAGE
    % **********************************************************************************************
    tic
    while stageBetOnTarget
        
        %--------------------------------------------------------------------------------------------
        % Unit code for checking eye position
        error = Eyelink('CheckRecording');
        if error
            fprintf(' *********** GOT AN ERROR IN LOCATEEYEPOSITION**********\n')
            iTrialOutcome = 'eyelinkError';
            stageFeedback = 1;
            break
        end
        if Eyelink( 'NewFloatSampleAvailable') > 0
            % get the sample in the form of an event structure
            event = Eyelink( 'NewestFloatSample');
            x = event.gx(eyeRecorded+1); % +1 as we're accessing MATLAB array
            y = event.gy(eyeRecorded+1);
            % do we have valid data and is the pupil visible?
            if x ~= el.MISSING_DATA && y ~= el.MISSING_DATA && event.pa(eyeRecorded+1) > 0
                eyeX = x;
                eyeY = y;
            end
        end
        %----------------------------------------------------------------------------------------------
        
        
        if inAcceptWindow(eyeX, eyeY, betWindowCoord(iBetIndex,:)) && toc >= postSaccadeHoldDuration
            % Flag whether it was correct
            if iBetIndex == iHighBet
                iBet = 'High';
            else
                iBet = 'Low';
            end
            [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
            switch iDecision
                case 'Correct'
                    if strcmp(iBet, 'High')
                        iTrialOutcome = 'correctHigh';
                    else
                        iTrialOutcome = 'correctLow';
                    end
                case 'Incorrect'
                    if strcmp(iBet, 'High')
                        iTrialOutcome = 'incorrectHigh';
                    else
                        iTrialOutcome = 'incorrectLow';
                    end
            end
            stageBetOnTarget = 0;
            stageFeedback = 1;
        elseif ~inAcceptWindow(eyeX, eyeY, betWindowCoord(iBetIndex,:)) && toc < postSaccadeHoldDuration
            % Flag it was bet aborted
            iTrialOutcome = 'betTargetHoldAbort';
            iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
            abortTrialFlag = 1;
            stageBetOnTarget = 0;
            stageFeedback = 1;
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % START HERE
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % **********************************************************************************************
    %         FEEDBACK STAGE
    % **********************************************************************************************
    while stageFeedback
        [window, timeout] = feedback(window, iTrialOutcome);
        [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
        iFeedbackOnset          = StimulusOnsetTime - iTrialOnsetComputerTime;
        iFeedbackDuration       = feedbackTime + timeout;
        
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
        decisionFixSpotOnset       = [decisionFixSpotOnset; iDecisionFixSpotOnset];
        decisionFixationOnset        = [decisionFixationOnset; iDecisionFixation];
        decisionFixSpotDuration    = [decisionFixSpotDuration; iDecisionFixSpotDuration];
        %         postTargetFixDuration   = [postTargetFixDuration; iPostMaskFixDuration];
        decisionTargOnset             = [decisionTargOnset; iDecisionTargOnset];
        targetIndex                 = [targetIndex; iTargetIndex];
        soa                     = [soa; iSOA];
        decisionMaskOnset               = [decisionMaskOnset; iDecisionMaskOnset];
        decisionResponseCueOnset        = [decisionResponseCueOnset; iDecisionResponseCueOnset];
        decisionResponseOnset        = [decisionResponseOnset; iDecisionResponseOnset];
        decisionIndex                 = [decisionIndex; iDecisionIndex];
        
        betFixSpotOnset           = [betFixSpotOnset; iBetFixSpotOnset];
        betFixationOnset           = [betFixationOnset; iBetFixation];
        preBetFixDuration           = [preBetFixDuration; iPreBetFixDuration];
        betTargetOnset           = [betTargetOnset; iBetTargetOnset];
        betResponseOnset           = [betResponseOnset; iBetResponseOnset];
        betIndex                 = [betIndex; iBetIndex];
        
        feedbackOnset           = [feedbackOnset; iFeedbackOnset];
        feedbackDuration        = [feedbackDuration; iFeedbackDuration];
        abortOnset              = [abortOnset; iAbortOnset];
        trialOutcome            = [trialOutcome; iTrialOutcome];
        
        targetAngle             = [targetAngle; targetAngleArray(iTarget)];
        highBetAngle             = [highBetAngle; betAngleArray(iHighBet)];
        lowBetAngle             = [lowBetAngle; betAngleArray(iLowBet)];
    end
    
    
    
    
    
    
    % Only increment trial number if it wasn't an unusable abort
    if ~abortTrialFlag && ~strcmp(iTrialOutcome, 'eyelinkError')
        iTrial                  = iTrial + 1;
    end
    
    if iTrial == totalTrial
        runningTask = 0;
    end
    
    
    %Take a break every "trialsPerBlock" trials
    if (mod(iTrial, trialsPerBlock) == 1) && (iTrial ~= 1) && ~abortTrialFlag && ~strcmp(iTrialOutcome, 'eyelinkError')
        iBlock = iBlock + 1;
        msg = sprintf('Press the space bar to begin block %d', iBlock);
        DrawFormattedText(window, msg, 'center', 'center', [0 200 0]);
        Screen('Flip', window);
        junky = NaN;
        WaitSecs(1);
        takeBreak = GetSecs;
        while isnan(junky) && GetSecs - takeBreak < 60;
            [kd, sec, kc] = KbCheck;
            if (kd==1) && (kc(32)==1)
                junky = 1;
            elseif (kd==1) && (kc(81)==1)
                clear junky
            end
            %             if kd && find(keyCode) == escapeKey
            %                 runningTask = 0;
            %                 junky = 1;
            %             end
        end
    end
    
end % trial loop


% ********************************************************************
%                    End Trial Loop
% ********************************************************************


nTrial = length(trialOnset);





















% ********************************************************************
%                    Trial Data
% ********************************************************************




% Location of Stimuli
% ---------------------------------------------------------------

fixationWindow = ones(nTrial, 1) * fixationWindow;
targetWindow = ones(nTrial, 1) * targetWindow;
maskWindow = ones(nTrial, 1) * maskWindow;
betWindow = ones(nTrial, 1) * betWindow;

fixationAngle           = ones(nTrial, 1) * fixationAngle;
fixationAmplitude     	= ones(nTrial, 1) * fixationAmplitude;
targetAmplitude =  ones(nTrial, 1) * targetAmplitude;
maskAmplitude =  ones(nTrial, 1) * maskAmplitude;


fixationSize              = ones(nTrial, 1) * fixationWidth;
targetSize       	= ones(nTrial, 1) * targetWidth;
maskSize          	= ones(nTrial, 1) * maskWidth;
betSize             = ones(nTrial, 1) * betWidth;

decisionFixColor 	= ones(nTrial, 1) * decisionFixColor;
targetColor       	= ones(nTrial, 1) * targetColor;
maskColor           = ones(nTrial, 1) * maskColor;
betFixColor        	= ones(nTrial, 1) * betFixColor;
highBetColor       	= ones(nTrial, 1) * highBetColor;
lowBetColor        	= ones(nTrial, 1) * lowBetColor;


% size(xxxxxxxxxxxxxxxxxx)


trialData = dataset(...
    {trialOutcome,          'trialOutcome'},...
    {round(trialOnset*1000),            'trialOnset'},...
    {round(trialDuration*1000),         'trialDuration'},...
    {round(abortOnset*1000),            'abortOnset'},...
    {round(decisionFixSpotOnset*1000), 	'decisionFixSpotOnset'},...
    {round(decisionFixationOnset*1000),  	'decisionFixationOnset'},...
    {round(decisionFixSpotDuration*1000),  	'decisionFixSpotDuration'},...
    {round(decisionTargOnset*1000),           'decisionTargOnset'},...
    {round(soa*1000),                   'soa'},...
    {round(decisionMaskOnset*1000),             'decisionMaskOnset'},...
    {round(decisionResponseCueOnset*1000), 'decisionResponseCueOnset'},...
    {round(decisionResponseOnset*1000),	'decisionResponseOnset'},...
    {round(betFixSpotOnset*1000),      	'betFixSpotOnset'},...
    {round(betFixationOnset*1000),           'betFixationOnset'},...
    {round(preBetFixDuration*1000),    	'preBetFixDuration'},...
    {round(betTargetOnset*1000),       	'betTargetOnset'},...
    {round(betResponseOnset*1000),     	'betResponseOnset'},...
    {round(feedbackOnset*1000),         'feedbackOnset'},...
    {round(feedbackDuration*1000),      'feedbackDuration'},...
    {targetIndex,      'targetIndex'},...
    {decisionIndex,      'decisionIndex'},...
    {betIndex,      'betIndex'},...
    {fixationAngle,      'fixationAngle'},...
    {targetAngle,        'targetAngle'},...
    {highBetAngle,       'highBetAngle'},...
    {lowBetAngle,        'lowBetAngle'},...
    {fixationSize,          'fixationSize'},...
    {targetSize,            'targetSize'},...
    {maskSize,              'maskSize'},...
    {betSize,               'betSize'},...
    {decisionFixColor,   	'decisionFixColor'},...
    {targetColor,           'targetColor'},...
    {maskColor,             'maskColor'},...
    {betFixColor,           'betFixColor'},...
    {highBetColor,          'highBetColor'},...
    {lowBetColor,           'lowBetColor'});


% ********************************************************************
% Session Data
% ********************************************************************


SessionData.taskID = 'retroMon';
SessionData.task.effector = 'eyeMovement';

SessionData.timing.year = num2str(year(now));
SessionData.timing.month = num2str(month(now));
SessionData.timing.day = num2str(day(now));
SessionData.timing.hour = num2str(clockVector(4));
SessionData.timing.minute = num2str(clockVector(5));


SessionData.timing.totalDuration = trialOnset(end) + trialDuration(end) - trialOnset(1); % seconds


SessionData.subjectID = subjectID;
SessionData.sessionID = session;

SessionData.stimuli.betFixColorRGB       = betFixColor;
SessionData.stimuli.highBetColorRGB   = highBetColor;
SessionData.stimuli.lowBetColorRGB   = lowBetColor;
SessionData.stimuli.maskColorRGB   = maskColor;
SessionData.stimuli.decisionFixColorRGB         = decisionFixColor;
SessionData.stimuli.targetRGB           = targetColor;






% ---------- Window Cleanup ----------
cleanup
Eyelink('CloseFile');
status = Eyelink('ReceiveFile',edffilename, 'data/');
status = Eyelink('ReceiveFile',edffilename);




if plotFlag
    %     post_session_psychometric(trialData)
end



if saveFlag
    save(saveFileName, 'trialData', 'SessionData');
    success = copy_human_files(edffilename, saveFileName);
end


end  % main function








function cleanup
Eyelink('stoprecording');
Screen('CloseAll');
ShowCursor;
end


% *******************************************************************
function [window, timeout] = feedback(window, iTrialOutcome)


backGround      = [60 60 60];
Screen('FillRect', window, backGround);
scrsz = get(0, 'ScreenSize');
screenWidth = scrsz(3);
screenHeight = scrsz(4);

incorrectTextColor = [250, 50, 50];
correctTextColor = [50, 220, 50];
Screen('TextFont', window, 'Times');
Screen('TextSize', window, 30);
Screen('TextStyle', window, 1);
if strcmp(iTrialOutcome, 'correctHigh')
    DrawFormattedText(window, '5', 'center', 'center', correctTextColor);
    timeout = 0;
elseif strcmp(iTrialOutcome, 'correctLow')
    DrawFormattedText(window, '3', 'center', 'center', correctTextColor);
    timeout = 0;
elseif strcmp(iTrialOutcome, 'incorrectHigh')
    DrawFormattedText(window, '-2', 'center', 'center', incorrectTextColor);
    timeout = 0;
elseif strcmp(iTrialOutcome, 'incorrectLow')
    DrawFormattedText(window, '2', 'center', 'center', incorrectTextColor);
    timeout = 0;
elseif strcmp(iTrialOutcome, 'decisionFixationAbort') || strcmp(iTrialOutcome, 'betFixationAbort')
    DrawFormattedText(window, 'Please stay fixated on center target', 'center', 'center', [190, 190, 190]);
    timeout = 2;
elseif strcmp(iTrialOutcome, 'decisionResponseAbort') || strcmp(iTrialOutcome, 'saccadeAbort') || ...
        strcmp(iTrialOutcome, 'betResponseAbort')
    DrawFormattedText(window, 'Hmmm, please look at one of the targets next time', 'center', 'center', [190, 190, 190]);
    timeout = 2;
elseif strcmp(iTrialOutcome, 'interstageAbort')
    DrawFormattedText(window, 'Please gaze at the fixation spot', 'center', 'center', [190, 190, 190]);
    timeout = 1;
elseif strcmp(iTrialOutcome, 'decisionTargetHoldAbort') || strcmp(iTrialOutcome, 'betTargetHoldAbort') || ...
        DrawFormattedText(window, 'Stay on the target until it disappears', 'center', 'center', incorrectTextColor);
    timeout = 3;
end

% Screen('Flip', window);
% waitSecs(feedbackTime + timeout)
end





% *******************************************************************
function inWindow = inAcceptWindow(eyeX, eyeY, acceptWindow)
% determine if gx and gy are within fixation window
inWindow = eyeX > acceptWindow(1) &&  eyeX <  acceptWindow(3) && ...
    eyeY > acceptWindow(2) && eyeY < acceptWindow(4) ;
end










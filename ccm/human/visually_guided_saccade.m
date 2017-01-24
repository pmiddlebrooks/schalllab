function visually_guided_saccade(targetAmplitude, nTrialPerLocation, saveFlag)

% stimulusAmplitude: in some units (e.g. pixels), the distance from the
%       center of the screen to the center of the stimulus
% stimulusAngle: in degrees, the angle from the center of the screen to the
%       center of the stimulus

% Format (numbering) of the squares in the stimulus (a 3 X 3 example):
%             0   1   2
%             3   4   5
%             6   7   8

% example SSDArrayScreenFlips = [22 26 30 34]






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Set up Experiment Variables, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stopKey=KbName('space');


% Get distance values and conversions for stimuli presentation
metersFromScreen    = .054;
screenWidthMeters   = .04;
theta               = asind(screenWidthMeters / 2 / sqrt(((screenWidthMeters / 2)^2) + metersFromScreen^2));
screenSize          = get(0, 'ScreenSize');
screenWidthPixels  	= screenSize(3);
screenHeightPixels 	= screenSize(4);
matlabCenterX       = screenWidthPixels/2;
matlabCenterY       = screenHeightPixels/2;
pixelsPerDegree     = screenWidthPixels / (2*theta);


% Get info for saving a data file.
% subjectName    = input('Enter initials of subject            ', 's');
subjectName =  'test';
clockVector = clock;
session     = [num2str(year(now)), '_', num2str(month(now)), '_', num2str(day(now)), '_', num2str(clockVector(4)), '_', num2str(clockVector(5))];
saveFileName = ['data/VGS_', subjectName, '_', session];


whichScreen     = 0;
backGround      = [0 0 0];
[window, centerPoint] = Screen('OpenWindow', whichScreen, backGround);

% commandwindow;
flipTime        = Screen('GetFlipInterval',window);% get the flip rate of current monitor.
preFlipBuffer   = flipTime / 2;
% a='Resolution & Refresh ok';
% if timing>.012 || timing<.011 %make sure we're running in 1024*768 at 85
% Hz, else stop
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
dummymode               = 0;
trialsPerBlock          = 100;
screenRefreshRate       = 1 / frameFrequency;
graceResponse           = 1.7; % seconds allowed to make a saccade
graceObtainFixation     = 2;  % seconds program will wait for user to obtain fixation
graceSaccadeDuration    = .1; % seconds allowed intra saccade time
postSaccadeHoldDuration = .4;  % duration to hold post-saccade fixation
stimulusScaleConversion = 10 / 100; % for now 10 pixels per every 100 pixels from fixation
feedbackTime            = 1.5;


% Fixation Spot constants
fixationWidth       = .7 * pixelsPerDegree;
fixationWindow      = [-fixationWidth * 1.5, -fixationWidth * 1.5, fixationWidth * 1.5, fixationWidth * 1.5];
fixationWindow      = CenterRect(fixationWindow, centerPoint);
fixationColor       = [200 200 200];
fixationHoldBase    = .5;
fixationHoldAdd     = (0 : 10 : 1000)  / 1000;



targetAmplitude     = targetAmplitude * pixelsPerDegree;
targetAngleArray    = repmat([0 45 90 135 180 225 270 315], 1, nTrialPerLocation);
targetWidth         = 1 + (targetAmplitude * stimulusScaleConversion);
targetWindowSize    = 4 * targetWidth;
targetColor         = [200 200 200];






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Set up Eyelink
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup Eyelink System- must have eyelink computer started and running
%eyelink
init = Eyelink('Initialize');
el = EyelinkInitDefaultsNC(window);
% % make sure that we get gaze data from the Eyelink
%     Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
edffilename = [subjectName,  datestr(now,'mmddyy'),  '.edf'];
openEDF = Eyelink('OpenFile', edffilename);
WaitSecs(0.1);
Eyelink('StartRecording');

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


if eyeRecorded == 0;
    Eyelink('Command', 'file_sample_filter=LEFT, GAZE, AREA, STATUS');
    Eyelink('Command', 'file_event_filter=LEFT,  FIXATION, SACCADE, BLINK, MESSAGE');
    Eyelink('Command', 'link_event_filter = LEFT, SACCADE,BLINK, MESSAGE');
else
    Eyelink('Command', 'file_sample_filter=RIGHT, GAZE, AREA, STATUS');
    Eyelink('Command', 'file_event_filter= RIGHT, FIXATION, SACCADE, BLINK, MESSAGE');
    Eyelink('Command', 'link_event_filter = RIGHT,SACCADE,BLINK, MESSAGE');
end

% screenpixstring = sprintf('screen_pixel_coords= %f,%f,%f,%f', centerPoint);      %you don't want eyelink to take those from physical.ini, because that means changing physical.ini for every screen resolution setting
% Eyelink('Command', 'clear_screen 0');

priorityLevel=MaxPriority(window);
Priority(priorityLevel);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           BEGIN TASK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% ***********************************************************************
%       INITIALIZE  STUFF
% ***********************************************************************
% Variables that will get filled each trial
trialOnset              = [];
trialDuration           = [];
fixationSpotOnset       = [];
fixationOnset           = [];
fixationSpotDuration    = [];
fixationHoldDuration    = [];
targetOnset             = [];
targetDuration          = [];
responseCueOnset        = [];
responseOnset           = [];
feedbackOnset           = [];
feedbackDuration        = [];
abortOnset              = [];

targetLocation          = [];






% ***********************************************************************
%     STIMULI POSITION INFORMATION
% ***********************************************************************

% ------------------------------
%     FIXATION STIMULUS
% ------------------------------
fixationEyeLinkX     = 0;
fixationEyeLinkY     = 0;
fixationLeft        = matlabCenterX - fixationWidth/2;
fixationTop         = matlabCenterY - fixationWidth/2;
fixationRight       = matlabCenterX + fixationWidth/2;
fixationBottom      = matlabCenterY + fixationWidth/2;
fixationSquare    = [fixationLeft fixationTop fixationRight fixationBottom];







% ***********************************************************************
%     START THE TRIAL LOOP
% ***********************************************************************
% Cue the user to begin the task by pressing the space bar
% DrawFormattedText(window, 'Press the space to begin.', 'center', 'center', [250 0 0]);
DrawFormattedText(window, 'About to start.', 'center', 'center', [250 0 0]);
Screen('Flip', window);
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
Waitsecs(2);
Screen('FillRect', window, backGround);
[vbl SOT] = Screen('Flip', window);



runningTask = 1;
iTrial = 0;
taskStartTime = GetSecs;
newTrialVariables = 1; % A flag, set to zero if an abort occurs and we want all trial variables to stay the same

while runningTask
    
    iTrialOnset             = GetSecs - taskStartTime;
    iTrialOnsetComputerTime = GetSecs;
    iTrial                  = iTrial + 1;
    % Initialize variables that may or may not get filled
    fixationAddIndex        = randperm(length(fixationHoldAdd));
    iFixationHoldDuration   = fixationHoldBase + fixationHoldAdd(fixationAddIndex(1));
    
    iTrialDuration = nan;
    iFixationSpotDuration = nan;
    iTargetOnset = nan;
    iTargetDuration = nan;
    iResponseCueOnset = nan;
    iResponseOnset          = nan;
    iAbortOnset             = nan;
    
    %     trialOnset              = [trialOnset; iTrialOnset];
    %     trialDuration           = [trialDuration; iTrialDuration];
    %     fixationSpotOnset       = [fixationSpotOnset; iFixationSpotOnset];
    %     fixationOnset           = [fixationOnset; iFixationOnset];
    %     fixationSpotDuration    = [fixationSpotDuration; iFixationSpotDuration];
    %     fixationHoldDuration    = [fixationHoldDuration; iFixationHoldDuration];
    %     targetOnset             = [targetOnset; iTargetOnset];
    %     targetDuration          = [targetDuration; iTargetDuration];
    %     responseCueOnset        = [responseCueOnset; iResponseCueOnset];
    %     responseOnset           = [responseOnset; iResponseOnset];
    %     feedbackOnset           = [feedbackOnset; iFeedbackOnset];
    %     feedbackDuration        = [feedbackDuration; iFeedbackDuration];
    %     abortOnset              = [abortOnset; iAbortOnset];
    
    
    %Take a break every "trialsPerBlock" trials
    if (mod(iTrial, trialsPerBlock) == 1) && (iTrial ~= 1)
        
        DrawFormattedText(w, 'Press the space to begin the next set.', 'center', centerY + 150, [250 0 0]);
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
        end
        Screen('FillRect', window, backGround);
        [vbl SOT] = Screen('Flip', window);
    end
    
    
    
    %     %Drift Correction
    %     drift = EyelinkDoDriftCorrect(el);
    
    %Start Recording Eye for current trial
    status = ['record_status_message "Trial ' num2str(iTrial) '"'];
    Eyelink('Command', status);
    Eyelink('Message', '%s%d', 'Trial=', iTrial);
    
    Screen('FillRect', window, backGround);
    [vbl SOT] = Screen('Flip', window);
    
    Eyelink('StartRecording');
    
    
    
    
    
    
    
    % If the previous trial was aborted before it began (user did not
    % obtain fixation), then keep all the variables the same. Only update
    % them if it's a new trial
    if newTrialVariables
        
        
        % ------------------------------
        %     TARGET STIMULUs
        % ------------------------------
        targetAngleIndexArray = randperm(length(targetAngleArray));
        targetAngleIndex = targetAngleIndexArray(1);
        targetAngle = targetAngleArray(targetAngleIndex);
        targetEyeLinkX      = targetAmplitude * cosd(targetAngle);
        targetEyeLinkY      = targetAmplitude  * sind(targetAngle);
        targetMatlabX      = matlabCenterX + targetEyeLinkX;
        targetMatlabY      = matlabCenterY - targetEyeLinkY;
        
        targetLeft         = targetMatlabX - targetWidth/2;
        targetTop          = targetMatlabY - targetWidth/2;
        targetRight        = targetMatlabX + targetWidth/2;
        targetBottom       = targetMatlabY + targetWidth/2;
        targetSquare     = [targetLeft, targetTop, targetRight, targetBottom];
        
    end
    
    
    %   WHICH IS THE TARGET?
    % -------------------------------
    iTrialTarget        = 1;
    iTargetLocation     = [targetMatlabX - matlabCenterX, matlabCenterY - targetMatlabY];
    iTargetWindow       = [targetMatlabX - targetWindowSize, targetMatlabY - targetWindowSize, targetMatlabX + targetWindowSize, targetMatlabY + targetWindowSize];
    
    
    
    
    
    
    
    
    
    
    
    % ********************************************
    %            BEGIN STAGES
    % ********************************************
    
    
    % Initialize stage logical variables each trial
    addTrialDataFlag    =  1;  % Gets set to false for pre-fixation aborts
    stagePreFixation    = 1;
    stageFixation       = 0;
    stageTargetOn       = 0;
    stageInFlight       = 0;
    stageOnTarget       = 0;
    stageFeedback       = 0;
    
    
    
    
    
    % ********************************************
    %            PRE-FIXATION STAGE
    % ********************************************
    % Turn on the fixation spot and wait unit subject fixates the spot- or,
    % if subject does not fixate spot within a grace period, start a new
    % trial
    Screen('FillRect', window, fixationColor, fixationSquare);
    Screen('FrameRect', window, fixationColor, fixationWindow);
    [~, fixationOnsetTime, ~, ~, ~] = Screen('Flip', window);
    iFixationSpotOnset = fixationOnsetTime - iTrialOnsetComputerTime;
    tic;
    while stagePreFixation
        
        % Unit code for checking eye position
        %-----------------------------------------------------------------------------------------------------
        if dummymode == 0
            error = Eyelink('CheckRecording');
            if(error~=0)
                disp('*********** GOT AN ERROR IN LOCATEEYEPOSITION**********')
                return;
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
                   % if data is valid, draw a circle on the screen at current gaze position
                    % using PsychToolbox's Screen function
                    gazeRect=[ x-9 y-9 x+10 y+10];
                    colour=round(rand(3,1)*255); % coloured dot
    Screen('FillRect', window, fixationColor, fixationSquare);
                    Screen('FillOval', window, colour, gazeRect);
                    Screen('Flip',  el.window, [], 1); % don't erase
                        eyeX = x;
                        eyeY = y;
                    end
                end
            end
        else
            
            % Query current mouse cursor position (our "pseudo-eyetracker") -
            % (mx,my) is our gaze position.
            [eyeX, eyeY, ~]=GetMouse(window); %#ok<*NASGU>
        end
        %-----------------------------------------------------------------------------------------------------
        
        
        if inAcceptWindow(eyeX, eyeY, fixationWindow) && toc <= graceObtainFixation
            Eyelink('Message', 'Fixation Start');
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
    
    
    
    
    
    
    
    
    
    % ********************************************
    %            FIXATION STAGE
    % ********************************************
    % Subject has begun fixating, can either hold fixation until the
    % targets come on and advance to stageTargetOn, or can abort the trial
    % and start a new trial by not maintaining fixation
    if stageFixation
        % Prepare the screen for stimuli onsets
        Screen('FillRect', window, targetColor, targetSquare);
        Screen('FrameRect', window, targetColor, iTargetWindow);
        Screen('FillRect', window, fixationColor, fixationSquare);
    Screen('FrameRect', window, fixationColor, fixationWindow);
tic
while stageFixation
            
            WaitSecs(.09)
            % Unit code for checking eye position
            %-----------------------------------------------------------------------------------------------------
            if dummymode == 0
                error = Eyelink('CheckRecording');
                if(error~=0)
                    disp('*********** GOT AN ERROR IN LOCATEEYEPOSITION**********')
                    return;
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
                [eyeX, eyeY, ~]=GetMouse(window); %#ok<*NASGU>
            end
            %-----------------------------------------------------------------------------------------------------
            
            
            if ~inAcceptWindow(eyeX, eyeY, fixationWindow) && toc < iFixationHoldDuration
                Eyelink('Message', 'Fixation Abort');
                iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
                iTrialOutcome       = 'fixationAbort';
                stageFixation       = 0;
                stageFeedback       = 1;
                newTrialVariables   = 0;
            elseif inAcceptWindow(eyeX, eyeY, fixationWindow) && toc >= iFixationHoldDuration
            status = Eyelink('Command', 'online_dcorr_trigger');
                stageFixation       = 0;
                stageTargetOn       = 1;
            end
            % We wait 1 ms each loop-iteration so that we
            % don't overload the system in realtime-priority:
            WaitSecs(0.001);
        end
    end
    
    
    
    
    if stageTargetOn
        [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
        iFixationSpotDuration   = StimulusOnsetTime - iTrialOnsetComputerTime - iFixationSpotOnset;
        iTargetOnset            = StimulusOnsetTime - iTrialOnsetComputerTime;
        iResponseCueOnset       = StimulusOnsetTime - iTrialOnsetComputerTime;
        Eyelink('Message', 'Targets Start');
        
        tic
        % **************************************************
        %         TARGET ON STAGE
        % **************************************************
        % All stimuli come on at once (in this version at least). In the case
        % of a go trial, wait until subject makes a saccade to one of the
        % targets or aborts by timing out. In the case of a stop trial, wait
        % until subject makes a saccade (an error) or
        while stageTargetOn
            
            
            % Unit code for checking eye position
            %-----------------------------------------------------------------------------------------------------
            if dummymode == 0
                error = Eyelink('CheckRecording');
                if(error~=0)
                    disp('*********** GOT AN ERROR IN LOCATEEYEPOSITION**********')
                    return;
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
                [eyeX, eyeY, ~]=GetMouse(window); %#ok<*NASGU>
            end
            %-----------------------------------------------------------------------------------------------------
            
            
            % Made a saccade
            if ~inAcceptWindow(eyeX, eyeY, fixationWindow)
                Eyelink('Message', 'Response Onset');
                iResponseOnset = GetSecs - iTrialOnsetComputerTime;
                stageTargetOn         = 0;
                stageInFlight       = 1;
                
                % Waited too long to make a saccade
            elseif inAcceptWindow(eyeX, eyeY, fixationWindow) && toc > graceResponse
                Eyelink('Message', 'Response Timed Out');
                iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
                iTrialOutcome       = 'timedOut';
                stageTargetOn         = 0;
                stageFeedback       = 1;
            end
        end
        
        
        newTrialVariables = 1;
    end
    
    
    
    
    
    % **************************************************
    %         IN FLIGHT STAGE
    % **************************************************
    tic
    while stageInFlight
        
        
        % Unit code for checking eye position
        %-----------------------------------------------------------------------------------------------------
        if dummymode == 0
            error = Eyelink('CheckRecording');
            if(error~=0)
                disp('*********** GOT AN ERROR IN LOCATEEYEPOSITION**********')
                return;
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
            [eyeX, eyeY, ~]=GetMouse(window); %#ok<*NASGU>
        end
        %-----------------------------------------------------------------------------------------------------
        
        
        if toc > graceSaccadeDuration
            iTrialOutcome       = 'saccadeAbort';
            iAbortOnset         = GetSecs - iTrialOnsetComputerTime;
            stageInFlight       = 0;
            stageFeedback       = 1;
        elseif inAcceptWindow(eyeX, eyeY, iTargetWindow) && toc <= graceSaccadeDuration
            stageInFlight       = 0;
            stageOnTarget       = 1;
        end
    end
    
    
    
    % **************************************************
    %         ON TARGET STAGE
    % **************************************************
    tic
    while stageOnTarget
        
        
        % Unit code for checking eye position
        %-----------------------------------------------------------------------------------------------------
        if dummymode == 0
            error = Eyelink('CheckRecording');
            if(error~=0)
                disp('*********** GOT AN ERROR IN LOCATEEYEPOSITION**********')
                return;
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
            [eyeX, eyeY, ~]=GetMouse(window); %#ok<*NASGU>
        end
        %-----------------------------------------------------------------------------------------------------
        
        
        
        if inAcceptWindow(eyeX, eyeY, iTargetWindow) && toc >= postSaccadeHoldDuration
            % Flag it was correct to the target
            iTrialOutcome = 'correct';
            % Update the remaining trial types
            targetAngleArray(targetAngleIndex) = [];
            stageOnTarget = 0;
            stageFeedback = 1;
        elseif ~inAcceptWindow(eyeX, eyeY, iTargetWindow) && toc < postSaccadeHoldDuration
            % Flag it was target aborted
            iTrialOutcome = 'targetHoldAbort';
            stageOnTarget = 0;
            stageFeedback = 1;
        end
    end
    
    
    
    
    
    
    % **************************************************
    %         FEEDBACK STAGE
    % **************************************************
    while stageFeedback
        [window, timeout] = feedback(window, iTrialOutcome);
        [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
        iTargetDuration         = StimulusOnsetTime - iTrialOnsetComputerTime + iTargetOnset;
        iFeedbackOnset          = StimulusOnsetTime - iTrialOnsetComputerTime;
        iFeedbackDuration       = feedbackTime + timeout;
        
        
        waitSecs(feedbackTime + timeout);
        stageFeedback = 0;
    end
    
    iTrialDuration = GetSecs - iTrialOnsetComputerTime;
    
    % Allow option to quit task by pressing escape at the end of a trial
    tic
    while toc < 1
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        if keyIsDown && find(keyCode) == stopKey
            runningTask = 0;
        end
    end
    if isempty(targetAngleArray)
        runningTask = 0;
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
        fixationHoldDuration    = [fixationHoldDuration; iFixationHoldDuration];
        targetOnset             = [targetOnset; iTargetOnset];
        targetDuration          = [targetDuration; iTargetDuration];
        responseCueOnset        = [responseCueOnset; iResponseCueOnset];
        responseOnset           = [responseOnset; iResponseOnset];
        feedbackOnset           = [feedbackOnset; iFeedbackOnset];
        feedbackDuration        = [feedbackDuration; iFeedbackDuration];
        abortOnset              = [abortOnset; iAbortOnset];
        
        
        % Stimulus Properties
        % ------------
        targetLocation          = [targetLocation; iTargetLocation];
    end
            Eyelink('StopRecording');

end % trial loop


nTrial = length(trialOnset);




% ********************************************************************
%                    Trial Data
% ********************************************************************


% Event Timing
% ------------




% Event Properties
% ---------------------------------------------------------------



% Location of Stimuli
% ---------------------------------------------------------------

fixationLocation          = [ones(nTrial, 1) * fixationEyeLinkX, ones(nTrial, 1) * fixationEyeLinkY];

fixationSize              = ones(nTrial, 1) * fixationWidth;
targetSize                = ones(nTrial, 1) * targetWidth;
fixationColor             = ones(nTrial, 1) * fixationColor;
targetColor               = ones(nTrial, 1) * targetColor;



trialData = dataset(...
    {trialOnset,            'trialOnset'},...
    {trialDuration,         'trialDuration'},...
    {abortOnset,            'abortOnset'},...
    {fixationSpotOnset,     'fixationSpotOnset'},...
    {fixationOnset,         'fixationOnset'},...
    {fixationHoldDuration   'fixationHoldDuration'},...
    {fixationSpotDuration,  'fixationSpotDuration'},...
    {targetOnset,           'targetOnset'},...
    {targetDuration,        'targetDuration'},...
    {responseCueOnset,      'responseCueOnset'},...
    {responseOnset,         'responseOnset'},...
    {feedbackOnset,         'feedbackOnset'},...
    {feedbackDuration,      'feedbackDuration'},...
    {fixationLocation,      'fixationLocation'},...
    {targetLocation,        'targetLocation'},...
    {fixationSize,          'fixationSize'},...
    {targetSize,            'targetSize'},...
    {fixationColor,         'fixationColor'},...
    {targetColor,           'targetColor'});


% ********************************************************************
% Session Data
% ********************************************************************


sessionData.task.effector = 'eyeMovement';

sessionData.timing.totalDuration = trialOnset(end) + trialDuration(end) - trialOnset(1); % seconds


sessionData.subjectID = subjectName;
sessionData.sessionID = session;







% ---------- Window Cleanup ----------

% Closes all windows.
Screen('CloseAll');

% Restores the mouse cursor.
ShowCursor;



if saveFlag
    save(saveFileName, 'trialData', 'sessionData');
end


end











% *******************************************************************
function [window, timeout] = feedback(window, iTrialOutcome)


scrsz = get(0, 'ScreenSize');
screenWidth = scrsz(3);
screenHeight = scrsz(4);

incorrectTextColor = [250, 50, 50];
correctTextColor = [50, 220, 50];
Screen('TextFont', window, 'Times');
Screen('TextSize', window, 50);
Screen('TextStyle', window, 1);
if strcmp(iTrialOutcome, 'correct')
    [nx, ny, bbox] = DrawFormattedText(window, 'Nice job', 'center', 'center', correctTextColor);
    timeout = 0;
elseif strcmp(iTrialOutcome, 'fixationAbort')
    DrawFormattedText(window, 'Need to stay fixated', 'center', 'center', incorrectTextColor);
    timeout = 1;
elseif strcmp(iTrialOutcome, 'timedOut')
    DrawFormattedText(window, 'You should have responded', 'center', 'center', incorrectTextColor);
    timeout = 2;
elseif strcmp(iTrialOutcome, 'saccadeAbort')
    DrawFormattedText(window, 'Hmmm, please look at the targets next time', 'center', 'center', [190, 190, 190]);
    timeout = 2;
elseif strcmp(iTrialOutcome, 'targetHoldAbort')
    DrawFormattedText(window, 'Stay on the target until it disappears', 'center', 'center', incorrectTextColor);
    timeout = 2;
end

% Screen('Flip', window);
% waitSecs(feedbackTime + timeout)
end





% *******************************************************************
function [eyeX, eyeY] = locateEyePosition(el, eyeRecorded, dummymode)
if dummymode == 0
    error = Eyelink('CheckRecording');
    if(error~=0)
        disp('*********** GOT AN ERROR IN LOCATEEYEPOSITION**********')
        return;
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
    [eyeX, eyeY, ~]=GetMouse(window); %#ok<*NASGU>
end

end

% *******************************************************************
function inWindow = inAcceptWindow(eyeX, eyeY, acceptWindow)
% determine if gx and gy are within fixation window
inWindow = eyeX > acceptWindow(1) &&  eyeX <  acceptWindow(3) && ...
    eyeY > acceptWindow(2) && eyeY < acceptWindow(4);
end










function task_images(soaArrayScreenFlips, nTarget, targetAmplitude, targetAngles, betAngle, targetDistribution, plotFlag, saveFlag)

% stimulusAmplitude: in some units (e.g. pixels), the distance from the
%       center of the screen to the center of the stimulus
% stimulusAngle: in degrees, the angle from the center of the screen to the
%       center of the stimulus


% example soaArrayScreenFlips = [22 26 30 34]


saveImages = 1;

if nargin == 0
    soaArrayScreenFlips = [2 4 6 8];
    targetAmplitude = 10;
    maskAmplitude = targetAmplitude;
    targetAngles = 45;
    betAngle = 45;
    nTarget = 4;
    targetDistribution = 'maximize';
    betTargetDistribution = 'mirror';
    plotFlag = 1;
    saveFlag = 1;
end

targetAngleArray    = get_mask_angle_array(targetAngles, nTarget, targetDistribution);
betAngleArray       = get_bet_angle_array(betAngle, betTargetDistribution);
% targetAngleArray = [45 135 225 315];
% betAngleArray = [45 135];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Set up Experiment Variables, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Keyboard assignments
KbName('UnifyKeyNames');
escapeKey               = KbName('ESCAPE');
backGround      = [60 60 60];
whichScreen     = 0;
[window, centerPoint] = Screen('OpenWindow', whichScreen, backGround);


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
maskAmplitudePixel     = targetAmplitudePixel;
betAmplitudePixel = targetAmplitude * pixelPerDegree;

% Get info for saving a data file.
clockVector     = clock;

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
BLOCKS_TO_RUN           = 1;
totalTrial              = trialsPerBlock * BLOCKS_TO_RUN;
soaArray                = soaArrayScreenFlips * (flipTime); %subtracting 2 ms each cycle to enusre the SSD occurs before the next screen refresh
graceResponse           = 1.7; % seconds allowed to make a saccade
graceObtainFixation     = 2;  % seconds program will wait for user to obtain fixation
graceSaccadeDuration    = .2; % seconds allowed intra saccade time
postSaccadeHoldDuration = .4;  % duration to hold post-saccade fixation
stimulusScaleConversion = 1 / 10; % for now 10 pixels per every 100 pixels from fixation
feedbackTime            = .4;
RETRO_FLAG              = 1;
PRO_FLAG                = 2;
PROP_RETRO             	= .5;

% Fixation Spot constants
fixationAngle       = 0;
fixationAmplitude   = 0;
fixationWidth       = .5;
fixationWidthPixel  = fixationWidth * pixelPerDegree;
fixWindowScale      = 3;  % fix window is 4 times size of fix point
fixationWindow      = fixationWidth*fixWindowScale;
fixationWindowPixel	= [-fixationWidthPixel*fixWindowScale, -fixationWidthPixel*fixWindowScale, fixationWidthPixel*fixWindowScale, fixationWidthPixel*fixWindowScale];
fixationWindowPixel = CenterRect(fixationWindowPixel, centerPoint);
decisionFixColor 	= [200 200 200];
perceptionFixColor 	= [200 200 200];
retroCueColor       = [0 120 170];
proCueColor         = [170 120 0];
betFixColor         = [0 200 0];
fixationHoldBase    = .5;
fixationHoldAdd     = (0 : 10 : 500)  / 1000;



targetWidth         = 1;
targetWidthPixel   	= targetWidth * pixelPerDegree;
targetWindowScale   = 4;
targetWindow        = targetWindowScale * targetWidth;
targetWindowPixel 	= targetWindow * pixelPerDegree;
targetColor         = [200 0 0];

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
%             BEGIN TASK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% ***********************************************************************
%       INITIALIZE STUFF
% ***********************************************************************
% Variables that will get filled each trial
trialOnset                  = [];
% trialDuration               = [];
% retroProFlag                = [];
% 
% perceptionFixSpotOnset  	= [];
% perceptionFixationOnset  	= [];
% perceptionFixSpotDuration   = [];
% perceptionTargOnset     	= [];
% soa                         = [];
% perceptionMaskOnset         = [];
% targetIndex                 = [];
% retroProCueOnset            = [];
% 
% betFixSpotOnset             = [];
% betFixationOnset            = [];
% preBetFixDuration           = [];
% betTargetOnset              = [];
% betResponseOnset            = [];
% betIndex                    = [];
% 
% decisionFixSpotOnset        = [];
% decisionFixationOnset       = [];
% decisionFixSpotDuration     = [];
% decisionMaskOnset           = [];
% % decisionResponseCueOnset = [];
% decisionResponseOnset       = [];
% 
% feedbackOnset               = [];
% feedbackDuration            = [];
% abortOnset                  = [];
% trialOutcome                = {};
% 
% targetAngle                 = [];
% highBetAngle                = [];
% lowBetAngle                 = [];



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
    
    
    
    
    % Initialize variables that may or may not get filled
    % Perception Stage
    iPreTargetFixSpotOnset 	= nan;
    iPreTargetFixation      = nan;
    fixationAddIndex        = randperm(length(fixationHoldAdd));
    iPreTargetFixDuration   = fixationHoldBase + fixationHoldAdd(fixationAddIndex(1));
    iPostMaskFixDuration    = fixationHoldBase + fixationHoldAdd(fixationAddIndex(2))*2;
    iPerceptionTargOnset  	= nan;
    iRealSOA                = nan;
    iSOA                    = nan;
    iPerceptionMaskOnset   	= nan;
    iTargetIndex            = nan;
    iRetroProCueOnset       = nan;
    
    % Bet Stage
    iBetFixSpotOnset     	= nan;
    iBetFixation            = nan;
    iPreBetFixDuration      = fixationHoldBase + fixationHoldAdd(fixationAddIndex(3));
    iBetTargetOnset         = nan;
    iBetResponseOnset       = nan;
    
    % Decision Stage
    iDecisionFixSpotOnset   = nan;
    iDecisionFixation       = nan;
    iPreMaskFixDuration     = fixationHoldBase + fixationHoldAdd(fixationAddIndex(3));
    iDecisionFixSpotDuration = nan;
    iDecisionMaskOnset     	= nan;
    iDecisionFixDuration    = fixationHoldBase + fixationHoldAdd(fixationAddIndex(1));
    iDecisionResponseCueOnset 	= nan;
    iDecisionResponseOnset    	= nan;
    iDecisionIndex          = nan;
    
    % Outcome
    iTrialOutcome           = nan;
    iAbortOnset             = nan;
    iFeedbackOnset          = nan;
    iFeedbackDuration       = nan;
    
    
    
    
    
    
    
    % If the previous trial was aborted before it began (user did not
    % obtain fixation), then keep all the variables the same. Only update
    % them if it's a new trial
    if newTrialVariables
        
        %   IS THIS A RETRO OR PRO TRIAL?
        % -------------------------------
        randomProportion = rand;
        if randomProportion > PROP_RETRO
            iRetroProFlag = RETRO_FLAG;
        else
            iRetroProFlag = PRO_FLAG;
        end
        
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
    addTrialDataFlag            = 1;  % Gets set to true for pre-fixation aborts
    abortTrialFlag              = 0;
    retroTrialFlag            	= 0;
    proTrialFlag              	= 0;
    
    stagePerceptionFix       	= 1;
    stagePerceptionTargetOn   	= 0;
    stagePerceptionMaskOn      	= 0;
    
    stageBetPreFix              = 0;
    stageBetPostFix             = 0;
    stageBetTargetOn            = 0;
    stageBetInFlight            = 0;
    stageBetOnTarget            = 0;
    
    stageDecisionPreFix         = 0;
    stageDecisionFix            = 0;
    stageDecisionMaskOn         = 0;
    stageCueToDecide            = 0;
    stageDecisionInFlight       = 0;
    stageDecisionOnTarget       = 0;
    
    stageFeedback               = 0;
    
    
    
    % ****************************************************************************************
    %            PRE-FIXATION STAGE
    % ****************************************************************************************
    % Turn on the fixation spot and wait unit subject fixates the spot- or,
    % if subject does not fixate spot within a grace period, start a new
    % trial
    

    Screen('FillRect', window, backGround);
    Screen('FillRect', window, perceptionFixColor, fixationSquare);
    [~, fixationOnsetTime, ~, ~, ~] = Screen('Flip', window);
    if saveImages
        image = Screen('GetImage', window);
        imwrite(image, 'perceptionFixation.tiff', 'tiff');
    end
    
    stagePerceptionFix = 1;
    
    
    iTrialOnset             = GetSecs - taskStartTime;
    iTrialOnsetComputerTime = GetSecs;
    
    
    
    
    
    iPerceptionFixSpotOnset = fixationOnsetTime - iTrialOnsetComputerTime;
    
    
    % ****************************************************************************************
    %            FIXATION STAGE
    % ****************************************************************************************
    % Subject has begun fixating, can either hold fixation until the
    % target comes on and advance to stageTargetOn, or can abort the trial
    % and start a new trial by not maintaining fixation

    if stagePerceptionFix
        
        % Prepare the screen for target onset
        Screen('FillRect', window, targetColor, iTargetSquare);
        Screen('FillRect', window, perceptionFixColor, fixationSquare);
        
        waitSecs(iPreTargetFixDuration)
        stagePerceptionTargetOn = 1;
        [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
        if saveImages
            image = Screen('GetImage', window);
            imwrite(image, 'perceptionTarget.tiff', 'tiff');
        end
        iPerceptionTargOnset            = StimulusOnsetTime - iTrialOnsetComputerTime;
    end
    
    
    % ****************************************************************************************
    %            TARGET ON STAGE
    % ****************************************************************************************
    % Subject has obtained fixation, and the targets appeared,
    % can either conitue to fixate until choice stimulus comes on, or can abort the trial
    % and start a new trial by not maintaining fixation

    if stagePerceptionTargetOn
        %         Screen('FillRect', window, targetColor, iTargetSquare);
        Screen('FillRect', window, perceptionFixColor, fixationSquare);
        % Load up each mask
        for iMask = 1 : nTarget
            Screen('FillRect', window, maskColor, maskSquare(iMask, :));
        end
        waitSecs(.1 - preFlipBuffer);
        [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
        if saveImages
            image = Screen('GetImage', window);
            imwrite(image, 'perceptionMask.tiff', 'tiff');
        end
        iPerceptionMaskOnset = StimulusOnsetTime - iTrialOnsetComputerTime;
        stagePerceptionTargetOn = 0;
        stagePerceptionMaskOn = 1;
    end
    
    
   
    
    
    
    % ****************************************************************************************
    %            MASK ON STAGE
    % ****************************************************************************************
    if stagePerceptionMaskOn
        % Load up each mask
        for iMask = 1 : nTarget
            Screen('FillRect', window, maskColor, maskSquare(iMask, :));
        end
        switch iRetroProFlag
            case RETRO_FLAG
                Screen('FillRect', window, retroCueColor, fixationSquare);
                retroTrialFlag = 1;
            case PRO_FLAG
                Screen('FillRect', window, proCueColor, fixationSquare);
                proTrialFlag = 1;
        end
        
            
            
        waitSecs(iPostMaskFixDuration)    
            [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
            switch iRetroProFlag
                case RETRO_FLAG
                    if saveImages
                        image = Screen('GetImage', window);
                        imwrite(image, 'retroCue.tiff', 'tiff');
                    end
                case PRO_FLAG
                    if saveImages
                        image = Screen('GetImage', window);
                        imwrite(image, 'proCue.tiff', 'tiff');
                    end
            end
            iPerceptionFixSpotDuration    = StimulusOnsetTime - iTrialOnsetComputerTime - iPerceptionFixSpotOnset;
            stagePerceptionMaskOn = 0;
    end






% **********************************************************************************************
% **********************************************************************************************
%                       RETROSPECTIVE MONITORING TRIAL
% **********************************************************************************************
% **********************************************************************************************
if retroTrialFlag
    stageCueToDecide = 1;
    iDecisionResponseCueOnset 	= StimulusOnsetTime - iTrialOnsetComputerTime;
    decInd = randi(2);
    switch decInd
        case 1
            iDecision = 'Correct';
        case 2
            iDecision = 'Incorrect';
    end
    
    if stageCueToDecide
        % **********************************************************************************************
        %         TARGETS AND STIMULI ON STAGE
        % **********************************************************************************************
        
            
            
                stageDecisionInFlight = 1;
    end
    
    
    
    
    
    % **********************************************************************************************
    %         IN FLIGHT STAGE
    % **********************************************************************************************
    if stageDecisionInFlight
        
                    stageDecisionOnTarget = 1;
    end
    
    
    
    % **********************************************************************************************
    %         ON TARGET STAGE
    % **********************************************************************************************
    Screen('FillRect', window, betFixColor, fixationSquare);

    if stageDecisionOnTarget
        
        waitSecs(graceObtainFixation) 
        
            [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
                    if saveImages
                        image = Screen('GetImage', window);
                        imwrite(image, 'retroBetFixation.tiff', 'tiff');
                    end
            stageBetPreFix = 1;
    end
    
    
    
    
    % **********************************************************************************************
    %         BET PRE FIXATION STAGE
    % **********************************************************************************************
    
    if stageBetPreFix
        
        
            stageBetPostFix = 1;
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
        
        
        if stageBetPostFix
            
                stageBetTargetOn = 1;
         WaitSecs(iPreBetFixDuration);
                [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
                   if saveImages
                        image = Screen('GetImage', window);
                        imwrite(image, 'retroBetTargets.tiff', 'tiff');
                    end
                iBetTargetOnset            = StimulusOnsetTime - iTrialOnsetComputerTime;
        end
    end
    
    
    % ****************************************************************************************
    %            BET TARGETS ON STAGE
    % ****************************************************************************************
    if stageBetTargetOn
        
                stageBetInFlight = 1;
    end
    
    % **********************************************************************************************
    %         IN FLIGHT STAGE
    % **********************************************************************************************
    if stageBetInFlight
        
        
                    stageBetOnTarget = 1;
    end
    
    
    
    % **********************************************************************************************
    %         ON TARGET STAGE
    % **********************************************************************************************
    
    if stageBetOnTarget
                           betInd = randi(2);
                   switch betInd
                       case 1
                           iBet = 'High';
                       case 2
                           iBet = 'Low';
                   end

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
         WaitSecs(postSaccadeHoldDuration);
        
            stageFeedback = 1;
    end
    
    
    
    
  
    
    
    
    % **********************************************************************************************
    % **********************************************************************************************
    %                       PROSPECTIVE MONITORING TRIAL
    % **********************************************************************************************
    % **********************************************************************************************
elseif proTrialFlag
        Screen('FillRect', window, betFixColor, fixationSquare);
         WaitSecs(graceObtainFixation);
                [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);

    iBetFixation      = GetSecs - iTrialOnsetComputerTime;
    stageBetPostFix = 1;
    
    
    
    
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

         WaitSecs(iPreBetFixDuration);
        [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
                   if saveImages
                        image = Screen('GetImage', window);
                        imwrite(image, 'proBetTargets.tiff', 'tiff');
                   end
                    stageBetTargetOn = 1;
    end
    
    
    % ****************************************************************************************
    %            BET TARGETS ON STAGE
    % ****************************************************************************************
    if stageBetTargetOn
        
            
                stageBetInFlight = 1;
    end
    
    % **********************************************************************************************
    %         IN FLIGHT STAGE
    % **********************************************************************************************

    if stageBetInFlight
        
        
                    stageBetOnTarget = 1;
    end
    
    
    
    % **********************************************************************************************
    %         ON TARGET STAGE
    % **********************************************************************************************
    Screen('FillRect', window, decisionFixColor, fixationSquare);
    
    
    if stageBetOnTarget
        
                      betInd = randi(2);
                   switch betInd
                       case 1
                           iBet = 'High';
                       case 2
                           iBet = 'Low';
                   end
     
          WaitSecs(graceObtainFixation);
       
            [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
                   if saveImages
                        image = Screen('GetImage', window);
                        imwrite(image, 'proDecisionFixation.tiff', 'tiff');
                   end
                    
            stageDecisionPreFix = 1;
    end
    
    
    
    
    
    if stageDecisionPreFix
            
            
                stageDecisionFix = 1;
     end
    
    
    
    % ****************************************************************************************
    %            FIXATION STAGE
    % ****************************************************************************************
    % Subject has begun fixating, can either hold fixation until the
    % masks come on and make decision, or can abort the trial
    % and start a new trial by not maintaining fixation
    if stageDecisionFix
        
        % Prepare the screen for target onset
        Screen('FillRect', window, decisionFixColor, fixationSquare);
        % Load up each mask
        for iMask = 1 : nTarget
            Screen('FillRect', window, maskColor, maskSquare(iMask, :));
        end
            
          WaitSecs(graceObtainFixation);
               [~, StimulusOnsetTime, ~, ~, ~] = Screen('Flip', window);
                   if saveImages
                        image = Screen('GetImage', window);
                        imwrite(image, 'proDecisionMasks.tiff', 'tiff');
                    end
       stageDecisionMaskOn = 1;     
        
    end
    
    
    
    
    
    
    
    % ****************************************************************************************
    %            MASK ON STAGE
    % ****************************************************************************************
    if stageDecisionMaskOn
               stageDecisionInFlight = 1;
    end
    
    
    
    
    
    
    % **********************************************************************************************
    %         IN FLIGHT STAGE
    % **********************************************************************************************
    
    if stageDecisionInFlight
        
        
        
        
                    stageDecisionOnTarget = 1;
     end
    
    
    
    % **********************************************************************************************
    %         ON TARGET STAGE
    % **********************************************************************************************
    
    
    if stageDecisionOnTarget
        
        
        
    decInd = randi(2);
    switch decInd
        case 1
            iDecision = 'Correct';
        case 2
            iDecision = 'Incorrect';
    end
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
          WaitSecs(graceObtainFixation);
            stageFeedback = 1;
    end
    
end













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



% Allow option to quit task by pressing escape at the end of a trial
tic
while toc < 1
    [ keyIsDown, seconds, keyCode ] = KbCheck;
    if keyIsDown && find(keyCode) == escapeKey
        runningTask = 0;
    end
end









% Only increment trial number if it wasn't an unusable abort
if ~abortTrialFlag && ~strcmp(iTrialOutcome, 'eyelinkError')
    iTrial  = iTrial + 1;
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
perceptionFixColor 	= ones(nTrial, 1) * perceptionFixColor;
targetColor       	= ones(nTrial, 1) * targetColor;
maskColor           = ones(nTrial, 1) * maskColor;
betFixColor        	= ones(nTrial, 1) * betFixColor;
highBetColor       	= ones(nTrial, 1) * highBetColor;
lowBetColor        	= ones(nTrial, 1) * lowBetColor;


% size(xxxxxxxxxxxxxxxxxx)

cleanup


function cleanup
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
elseif strcmp(iTrialOutcome, 'decisionFixationAbort') || strcmp(iTrialOutcome, 'betFixationAbort') || strcmp(iTrialOutcome, 'perceptionFixationAbort')
    DrawFormattedText(window, 'Please stay fixated on center target', 'center', 'center', [190, 190, 190]);
    timeout = 2;
elseif strcmp(iTrialOutcome, 'decisionResponseAbort') || strcmp(iTrialOutcome, 'saccadeAbort') || ...
        strcmp(iTrialOutcome, 'betResponseAbort')
    DrawFormattedText(window, 'Hmmm, please look at one of the targets next time', 'center', 'center', [190, 190, 190]);
    timeout = 2;
elseif strcmp(iTrialOutcome, 'proInterstageAbort') || strcmp(iTrialOutcome, 'retroInterstageAbort')
    DrawFormattedText(window, 'Please gaze at the fixation spot', 'center', 'center', [190, 190, 190]);
    timeout = 1;
elseif strcmp(iTrialOutcome, 'decisionTargetHoldAbort') || strcmp(iTrialOutcome, 'betTargetHoldAbort') || ...
        DrawFormattedText(window, 'Stay on the target until it disappears', 'center', 'center', incorrectTextColor);
    timeout = 3;
end

% Screen('Flip', window);
% waitSecs(feedbackTime + timeout)
end



end







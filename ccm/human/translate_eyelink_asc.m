function [Experiment] = translate_eyelink_asc(sessionASC, pixelPerDegreeX, pixelPerDegreeY)


fid = fopen(sessionASC);


PIXEL_PER_DEGREE_X = pixelPerDegreeX;
PIXEL_PER_DEGREE_Y = pixelPerDegreeY;

% data samples that get itteratively collected in each trial
trialEyeX = [];
trialEyeY = [];


% session data that gets itteratively collected accross trials, stuffed
% into a dataset for output

xEye = {};
yEye = {};
eyelinkTrialStart = [];
eyelinkTrialEnd = [];
eyelinkTargetOnset = [];
eyelinkChoiceStimulusOnset = [];
eyelinkStopSignalOnset = [];
eyelinkResponseOnset = [];
eyelinkFixationAbort = [];
eyelinkResponseTimedOut = [];

eyelinkSaccBegin = {};
eyelinkSaccadeEnd = {};
eyelinkFixationOnset = {};
eyelinkFixationEnd = {};

startedFlag = 0;

lastLine = 0;
sessionTrial = 0;
while ~lastLine
    iLine = fgetl(fid);
    if iLine == -1
        lastLine = 1;
    else
        [token, remain] = strtok(iLine);
        
        
        %  ********** SAMPLE DATA  **********
        % If the line starts with a number, it's sample data (eye position)
        if ~isempty(str2num(token))
            
            
            iSampleArray    = textscan(iLine, '%d %f %f %f');
            iSessionTime    = iSampleArray{1};
            iEyeXSample   	= iSampleArray{2};
            iEyeYSample    	= iSampleArray{3};
            iPupilSize      = iSampleArray{4};
            
            trialEyeX       = [trialEyeX; iEyeXSample];
            trialEyeY       = [trialEyeY; iEyeYSample];
            
            % Use mean of first 5 trial samples (20 ms) to center the eye traces
            if length(trialEyeX) == 5
                xEyeOffset = mean(trialEyeX);
                yEyeOffset = mean(trialEyeY);
            end
        else
            
            
            %  ********** EVENT DATA  **********
            % If the line starts with a string, it's an event
            switch token
                case 'START'
                    [timeStamp, remain] = strtok(remain);
                    %Initialize a bunch of stuff to get concatenated when a
                    %trial ends
                    trialEyeX = [];
                    trialEyeY = [];
                    jTrialStart = nan;
                    jTrialEnd = nan;
                    jTargetOnset = nan;
                    jChoiceStimulusOnset = nan;
                    jStopSignalOnset = nan;
                    jResponseOnset = nan;
                    jFixationAbort = nan;
                    jResponseTimedOut = nan;
                    jeyelinkSaccBegin = nan;
                    jEyelinkSaccadeEnd = nan;
                    jEyelinkFixationOnset = nan;
                    jEyelinkFixationEnd = nan;
                    
                    jTrialStart = str2num(timeStamp);
                    
                    % Reset some variables for new trial
                    if startedFlag
                        disp(sessionTrial)
                    else
                    startedFlag = 1;
                    end
                    
                    
                case 'END'
                    [timeStamp, remain] = strtok(remain);
                    
                    jTrialEnd = str2num(timeStamp) - jTrialStart;
                    
                    % On some trials eyelink has to calibrate again, and
                    % this causes matlab to skip the trial but Eyelink
                    % counts it as a very short trial (with a START and
                    % STOP). Don't include those trials, since they will
                    % misalign all the data from that trial forward.
                    if jTrialEnd > 300
                    sessionTrial = sessionTrial + 1;
                    % Concatenate the info from the trial:
                    
                    % Add trial data to the session data collection
                    trialEyeX = (trialEyeX - xEyeOffset) ./ PIXEL_PER_DEGREE_X;
                    trialEyeY = (trialEyeY - yEyeOffset) ./ PIXEL_PER_DEGREE_Y;
                    xEye = [xEye; trialEyeX];
                    yEye = [yEye; trialEyeY];
                    
                    eyelinkTrialStart = [eyelinkTrialStart; jTrialStart];
                    eyelinkTrialEnd = [eyelinkTrialEnd; jTrialEnd];
                    eyelinkTargetOnset = [eyelinkTargetOnset; jTargetOnset];
                    eyelinkChoiceStimulusOnset = [eyelinkChoiceStimulusOnset; jChoiceStimulusOnset];
                    eyelinkStopSignalOnset = [eyelinkStopSignalOnset; jStopSignalOnset];
                    eyelinkResponseOnset = [eyelinkResponseOnset; jResponseOnset];
                    eyelinkFixationAbort = [eyelinkFixationAbort; jFixationAbort];
                    eyelinkResponseTimedOut = [eyelinkResponseTimedOut; jResponseTimedOut];
                    eyelinkSaccBegin = [eyelinkSaccBegin; jeyelinkSaccBegin];
                    eyelinkSaccadeEnd = [eyelinkSaccadeEnd; jEyelinkSaccadeEnd];
                    eyelinkFixationOnset = [eyelinkFixationOnset; jEyelinkFixationOnset];
                    eyelinkFixationEnd = [eyelinkFixationEnd; jEyelinkFixationEnd];
                    
                    %                 if sessionTrial > 2
                    %                     plot(trialEyeX, trialEyeY)
                    %                     pause
                    %                     clf
                    %                 end
                     endedFlag = 1;
                     startedFlag = 0;
                    else
                        fprintf('Calibration at trial %d, discarding it\n', sessionTrial)
                    end
               case 'MSG'
                    [timeStamp, remain] = strtok(remain);
                    messageTime           = str2num(timeStamp);
                    % If it's a message, get the message- might want to use
                    % some of these
                    [jMSG, remain]      = strtok(remain);
                    if strcmp(jMSG, 'Fixation_Start')
                        % do something
                        % etc
                    elseif strcmp(jMSG, 'Targets_On')
                        % do something
                        jTargetOnset = messageTime - jTrialStart;
                    elseif strcmp(jMSG, 'Choice_Stimulus_On')
                        % do something
                        jChoiceStimulusOnset = messageTime - jTrialStart;
                    elseif strcmp(jMSG, 'Stop_Signal_On')
                        % do something
                        jStopSignalOnset = messageTime - jTrialStart;
                    elseif strcmp(jMSG, 'Response_Onset')
                        % do something
                        jResponseOnset = messageTime - jTrialStart;
                    elseif strcmp(jMSG, 'Fixation_Abort')
                        % do something
                        jFixationAbort = messageTime - jTrialStart;
                    elseif strcmp(jMSG, 'Response_Timed_Out')
                        % do something
                        jResponseTimedOut = messageTime - jTrialStart;
                    elseif strcmp(jMSG, 'xxxxxxx')
                        % do something
                        % etc
                    elseif strcmp(jMSG, 'xxxxxxx')
                        % do something
                        % etc
                    elseif strcmp(jMSG, 'xxxxxxx')
                        % do something
                        % etc
                    elseif strcmp(jMSG, 'xxxxxxx')
                        % do something
                        % etc
                    end
                    
                    
                case 'SFIX'
                    [eyeRecorded, remain] = strtok(remain);
                    [timeStamp, remain] = strtok(remain);
                    eventTime           = str2num(timeStamp);
                    % replace NaN with time if this is first fixation, else
                    % keep adding fixations for the trial
                    if isnan(jEyelinkFixationOnset)
                        jEyelinkFixationOnset = eventTime - jTrialStart;
                    else
                        jEyelinkFixationOnset = [jEyelinkFixationOnset, eventTime - jTrialStart];
                    end
                case 'EFIX'
                    [eyeRecorded, remain] = strtok(remain);
                    [timeStamp, remain] = strtok(remain);
                    eventTime           = str2num(timeStamp);
                    % xxxxxx
                    if isnan(jEyelinkFixationEnd)
                        jEyelinkFixationEnd = eventTime - jTrialStart;
                    else
                        jEyelinkFixationEnd = [jEyelinkFixationEnd, eventTime - jTrialStart];
                    end
                    
                case 'SSACC'
                    [eyeRecorded, remain] = strtok(remain);
                    [timeStamp, remain] = strtok(remain);
                    eventTime           = str2num(timeStamp);
                    % xxxxxx
                    if isnan(jeyelinkSaccBegin)
                        jeyelinkSaccBegin = eventTime - jTrialStart;
                    else
                        jeyelinkSaccBegin = [jeyelinkSaccBegin, eventTime - jTrialStart];
                    end
                    
                case 'ESACC'
                    [eyeRecorded, remain] = strtok(remain);
                    [timeStamp, remain] = strtok(remain);
                    eventTime           = str2num(timeStamp);
                    % xxxxxx
                    if isnan(jEyelinkSaccadeEnd)
                        jEyelinkSaccadeEnd = eventTime - jTrialStart;
                    else
                        jEyelinkSaccadeEnd = [jEyelinkSaccadeEnd, eventTime - jTrialStart];
                    end
                    
                case 'SBLINK'
                    % xxxxxx
                    
                case 'EBLINK'
                    % xxxxxx
                    
                case 'PRESCALER'
                    % xxxxxx
                    
                case 'VPRESCALER'
                    % xxxxxx
                    
                case 'EVENTS'
                    % xxxxxx
                    
                case 'SAMPLES'
                    % xxxxxx
                    
                case 'XXXXX'
                    % xxxxxx
                    
                case 'XXXXX'
                    % xxxxxx
                otherwise
                    %                 fprintf('Line %d: token: %s not recognized/n', i, token)
                    %                 return
            end % switch token
        end % if ~isempty(str2num(token))
    end % iLine == -1
    
end


eyeSampleHz = 500;

% Get saccade data.
nTrial = size(eyelinkTrialEnd, 1);
eyelinkTrialStartZero = 0;

saccBegin           = {};
saccToTargIndex   = [];
saccDuration   = {};
saccAmp   = {};
saccAngle   = {};
for jTrial = 1 : nTrial
    

    [jsaccBegin,...
        jSaccadeEnd,...
        jsaccAmp,...
        jsaccAngle,...
        SaccVelocity,...
        jsaccDuration,...
        SaccsNBlinks, ...
        jSaccadeToTargetIndex] = saccade_data(xEye{jTrial},...
        yEye{jTrial},...
        eyelinkTrialStartZero, ...
        eyelinkChoiceStimulusOnset(jTrial), ...
        eyelinkTrialEnd(jTrial), ...
        eyelinkTrialEnd(jTrial), ...
        3, ...
        eyeSampleHz);
    
    if isempty(jsaccBegin)
        jsaccBegin = {nan};
        jsaccAmp = {nan};
        jsaccAngle = {nan};
        jsaccDuration = {nan};
        jSaccadeToTargetIndex = [nan nan];
    end
       
    saccBegin           = [saccBegin; jsaccBegin];
    saccToTargIndex   = [saccToTargIndex; jSaccadeToTargetIndex(2)];
    saccDuration        = [saccDuration; jsaccDuration];
    saccAmp       = [saccAmp; jsaccAmp];
    saccAngle       = [saccAngle; jsaccAngle];

%     if ~isnan(jSaccadeToTargetIndex(2)) && strcmp(trialData.trialOutcome{jTrial}, 'stopCorrect')
% disp([jTrial, jSaccadeToTargetIndex(2)])
% disp(jsaccBegin)
% plot(xEye{jTrial}, yEye{jTrial})
% pause
%     end
end




Experiment = dataset(...
    {xEye,          'eyeX'},...
    {yEye,          'eyeY'},...
    {eyelinkTrialStart,          'eyelinkTrialStart'},...
    {eyelinkTrialEnd,          'eyelinkTrialEnd'},...
    {eyelinkTargetOnset,          'eyelinkTargetOnset'},...
    {eyelinkChoiceStimulusOnset,          'eyelinkChoiceStimulusOnset'},...
    {eyelinkStopSignalOnset,          'eyelinkStopSignalOnset'},...
    {saccToTargIndex,	'saccToTargIndex'},...
    {eyelinkFixationAbort,          'eyelinkFixationAbort'},...
    {eyelinkResponseOnset,         'eyelinkResponseOnset'},...
    {saccBegin,         'saccBegin'},...
    {saccDuration,      'saccDuration'},...
    {saccAmp,     'saccAmp'},...
    {saccAngle,     'saccAngle'},...
    {eyelinkSaccBegin,         'eyelinkSaccBegin'},...
    {eyelinkSaccadeEnd,         'eyelinkSaccadeEnd'},...
    {eyelinkFixationOnset,         'eyelinkFixationOnset'},...
    {eyelinkFixationEnd,         'eyelinkFixationEnd'},...
    {eyelinkResponseTimedOut,         'eyelinkResponseTimedOut'});

%     {xxxxxxx,         'xxxxxxxx'},...








function trialData = convert_data_file(sessionID)
%%
sessionID = 'np14n2';
sessionID = 'sp154n1';

% Determine subject and a place to save the converted data file
switch sessionID(1)
    case 'n'
        subjectID = 'nebby';
    case 's'
        subjectID = 'shuffles';
end
saveDir     = fullfile(local_data_path, subjectID);


% Alter sessionID to distinguis it from older versions of translated files
switch length(sessionID)
    case 6
        saveID = [sessionID(1:2),'0',sessionID(3:5),'0',sessionID(6)];
    case 7
        saveID = [sessionID(1:6),'0',sessionID(7)];
end





% load the original translated file
rexDataPath = 'rexroutines/data/';
rData = load(fullfile(rexDataPath,sessionID), 'allh', 'allv','allcodes','alltimes','allstart','allspk');

aData = load(fullfile(rexDataPath,[sessionID,'_analyses.mat']), 'target_angles', 'target_amp', 'bet_angles', 'brain_area','RF','hemisphere','targnum');
targAmp     = aData.target_amp/10;
maskAngles	= aData.target_angles;
betAngles   = aData.bet_angles;
nTarg       = aData.targnum;




betAngleLeftLog     = (betAngles < -89) & (betAngles > -270) | ...
    (betAngles > 90) & (betAngles < 269);
betAngleLeft        = betAngles(betAngleLeftLog);
betAngleRight       = betAngles(~betAngleLeftLog);
% betAngleLeftInd = find(betAngles == betAngleLeft);


% Mask target locations: this will (I think) need to be distinguished
% throuhout the code- 2 targets vs. 4 targets
leftTargInd = (maskAngles < -89) & (maskAngles > -270) | ...
    (maskAngles > 90) & (maskAngles < 269);
switch nTarg
    case 2
        decAngleLeft = maskAngles(leftTargInd);
        decAngleRight = maskAngles(~leftTargInd);
    case 4
        decAngleLeftUp = max(maskAngles(leftTargInd));
        decAngleLeftDown = min(maskAngles(leftTargInd));
        decAngleRightUp = max(maskAngles(~leftTargInd));
        decAngleRightDown = min(maskAngles(~leftTargInd));
end



% Info related to the session
SessionData.taskID          = 'maskbet';
SessionData.brainArea       = aData.brain_area;
SessionData.brainHemisphere = aData.hemisphere;
SessionData.receptiveField  = aData.RF;
SessionData.maskAngle       = maskAngles;
SessionData.betAngle        = betAngles;
switch subjectID
    case 'nebbey'
SessionData.recessed        = aData.recessed;
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEBUGSACC   = 0;
windowScale = .8; % allow 1 degree accept window per 1 degree target amplitude
msPerFrame 	= 1000/60; % The back-lit projector was 60 Hz
maxRTAllowed = 1500;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define time-stamped event codes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


eBkgOff        = 1026; % background light turned off
eDecTargOnPD    = 1028; % decision target turned on: in Retro.d, this occurs 1 ms after decTargOn (66xx), so it is essentially redundant
dBetTargOnPD    = 1029; % bet target turned on: in Retro.d, this occurs 1 ms after decTargOn (75xx), so it is essentially redundant
eDecTarg        = 1030; % a correct decision saccade (to target) -- also the first reward click (delivery)
eDecDist        = 1037; % an incorrect decision saccade (to a distractor)
eBetHigh        = 1050; % a high bet
eBetLow         = 1057; % a low bet
eRewardOn       = 1030; % the first click of reward
eRewardClick    = 1031; % each successive click
eTimeoutOn      = 1037;
eTimeoutOff     = 1038;
eAbortDec       = 1035; % Aborted during decision stage
eAbortBet       = 1036; % Aborted during bet stage



base                = [6020 6021 6022 6023 6024 6040 6041 6042 6043 6044];
eDecDownLeft        = base;
eDecDownRight       = eDecDownLeft + 5;
eDecUpLeft          = eDecDownLeft + 10;
eDecUpRight         = eDecDownLeft + 15;
eBetHighLeft       	= [6020 : 6039];
eBetHighRight     	= [6040 : 6059];

eDecFixOn            =  200;  % 62xx
eDecEyeInFixWindow   = 400;  % 64xx
eDecTargOn           = 600;  % 66xx
eDecFixOff           = 800;  % 68xx
eDecSaccOnset        = 1000;  % 70xx


eDecSaccIncorr       = 1100;  % 71xx
eDecSaccCorr         = 1200;  % 72xx

eBetFixOn            = 1300;  % 73xx
eBetEyeInFixWindow   = 1400;  % 74xx
eBetTargOn           = 1500;  % 75xx
eBetSaccOnset        = 1600;  % 76xx
eBetSaccMade         = 1700;  % 77xx
eBetSaccLow          = 1800;  % 78xx
eBetSaccHigh         = 1900;  % 79xx

% eBetFixOff           =  800;  % 68xx




nTrial  = size(rData.allh, 1);
td      = dataset();



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize matrices for events
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize matrices for task variables (non-timed "events")
td.decFixAngle     = zeros(nTrial,1);
td.decFixAmp       = zeros(nTrial,1);
td.decFixSize      = nan(nTrial,1);
td.decFixWindow    = ones(nTrial,1) * 3 + targAmp/5; % allow 5 degree fixation window for now
td.decTargAngle    = nan(nTrial,1);
td.decSaccAngle    = nan(nTrial,1);
td.decTargAmp      = ones(nTrial,1) * targAmp;
td.decTargSize     = nan(nTrial,1);
td.decMaskSize     = nan(nTrial,1);
td.decTargWindow   = ones(nTrial,1) * targAmp * windowScale; % target accept window is function of its amplitude
td.decMaskAngle    = repmat(maskAngles,nTrial,1); % All mask angles, including target and distractors
td.decMaskAmp      = ones(nTrial,1) * targAmp; %
td.betFixAngle     = zeros(nTrial,1);
td.betFixAmp       = zeros(nTrial,1);
td.betFixSize      = nan(nTrial,1);
td.betFixWindow    = ones(nTrial,1) * 3 + targAmp/5; % allow 5 degree fixation window for now
td.betHighAngle    = nan(nTrial,1);
td.betHighAmp      = ones(nTrial,1) * targAmp;
td.betLowAngle     = nan(nTrial,1);
td.betLowAmp       = ones(nTrial,1) * targAmp;
td.betSaccAngle  	= nan(nTrial,1);
td.betTargSize     = nan(nTrial,1);
td.betTargWindow   = ones(nTrial,1) * targAmp * windowScale; % target accept window is function of its amplitude

td.decOutcome       = cell(nTrial, 1);
td.betOutcome       = cell(nTrial, 1);
td.timeoutOn        = nan(nTrial,1);
td.timeoutOff       = nan(nTrial,1);

td.trialType        = repmat({'retro'}, nTrial, 1);
% Initialize matrices for timed events
td.decFixOn        = nan(nTrial,1);
td.decTargOn       = nan(nTrial,1);
td.decMaskOn       = nan(nTrial,1);
td.decFixOff       = nan(nTrial,1);
td.decResponseCueOn  = nan(nTrial,1);
td.decRT            = nan(nTrial, 1);
td.betFixOn        = nan(nTrial,1);
td.betTargOn       = nan(nTrial,1);
td.betResponseCueOn  = nan(nTrial,1);
td.betRT            = nan(nTrial, 1);
td.betFixOff       = nan(nTrial,1);
td.rewardOn        = num2cell(nan(nTrial,1));  % reward was delivered with multiple clicks


td.decFixWindowEntered     = nan(nTrial,1);
td.decTargWindowEntered    = nan(nTrial,1);
td.betFixWindowEntered     = nan(nTrial,1);
td.betTargWindowEntered    = nan(nTrial,1);


% Eye movement data
td.eyeX                = cell(nTrial,1);
td.eyeY                = cell(nTrial,1);
td.saccBegin           = cell(nTrial,1);
td.saccEnd             = cell(nTrial,1);
% td.saccDuration             = cell(nTrial,1);
td.saccAmp             = cell(nTrial,1);
td.saccAngle         	= cell(nTrial,1);
td.decSaccIndex         = nan(nTrial,1);
td.betSaccIndex         = nan(nTrial,1);
td.decResponseOnset     = nan(nTrial,1);
td.betResponseOnset     = nan(nTrial,1);

% Want to get rid of zeros at the end of all these trial vectors. First,
% flip the matrices to put zeros in front (will make it easier to remove)
allh        = fliplr(rData.allh);
allv        = fliplr(rData.allv);
allspk      = fliplr(rData.allspk);
td.spikeData        	= cell(nTrial,1);


baseRange = [20:24,40:44];












% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP THROUGH TRIALS AND FILL IN EVENTS, EYE MOVEMENTS, AND SPIKE DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : nTrial
    
    iCode = rData.allcodes(i,:)';
    iTime = rData.alltimes(i,:)';
    
    iBaseCode = iCode(iCode > 6019 & iCode < 6060);
    iDirCode = mod(iBaseCode,6000);
    
    
    
    
    
    
    
    
    % ==================================================
    % Non-timed Events data
    % ==================================================
    
    
    % Decision stage outcome
    if sum(iCode == eDecTarg)
        td.decOutcome{i} = 'target';
    elseif sum(iCode == eDecDist)
        td.decOutcome{i} = 'distractor';
    elseif sum(iCode == eAbortDec)
        td.decOutcome{i} = 'abort';
    end
    % Bet stage outcome
    if sum(iCode == eBetHigh)
        td.betOutcome{i} = 'high';
    elseif sum(iCode == eBetLow)
        td.betOutcome{i} = 'low';
    elseif sum(iCode == eAbortBet)
        td.betOutcome{i} = 'abort';
    end
    
    
    % High and low bet angles
    if iDirCode < 40
        td.betHighAngle(i) = betAngleLeft;
        td.betLowAngle(i) = betAngleRight;
    else
        td.betHighAngle(i) = betAngleRight;
        td.betLowAngle(i) = betAngleLeft;
    end
    
    
    % Target angle (all other mask angles are in td.decMaskAngle
    switch nTarg
        case 2
            if ismember(iDirCode, [baseRange baseRange + 10])
                td.decTargAngle(i) = decAngleLeft;
            elseif ismember(iDirCode, [baseRange + 5 baseRange + 15])
                td.decTargAngle(i) = decAngleRight;
            end
        case 4
            if ismember(iDirCode, baseRange)
                td.decTargAngle(i) = decAngleLeftDown;
            elseif ismember(iDirCode, baseRange + 5)
                td.decTargAngle(i) = decAngleRightDown;
            elseif ismember(iDirCode, baseRange + 10)
                td.decTargAngle(i) = decAngleLeftUp;
            elseif ismember(iDirCode, baseRange + 15)
                td.decTargAngle(i) = decAngleRightUp;
            end
    end
    
    
    
    
    % ==================================================
    % Timed Events data
    % ==================================================
    
    
    % Determine SOA for the trial: SOA 1-4 are 1-4 screen flips the masks
    % come on, and soa 5 is 5 screen flips and target goes of (masks do not
    % come on- essentially a memory-guided saccade trial)
    iSoaCode = 1 + mod(iDirCode, 10); % SOA screen frames = 1 + (Isolate the 1's digit of the basecode)
    if iSoaCode > 5
        iSoaCode = iSoaCode - 5;
    end

    
    % Decision Stage
    t 	= iTime((iCode >= iBaseCode + eDecFixOn) & (iCode < iBaseCode + eDecFixOn + 100));
    if ~isempty(t), td.decFixOn(i)   = t; end
    t  	= iTime((iCode >= iBaseCode + eDecEyeInFixWindow) & (iCode < iBaseCode + eDecEyeInFixWindow + 100));
    if ~isempty(t), td.decFixWindowEntered(i)   = t; end
    t  	= iTime((iCode >= iBaseCode + eDecTargOn) & (iCode < iBaseCode + eDecTargOn + 100));
    if ~isempty(t), td.decTargOn(i)   = t; end
    t  	= ceil(td.decTargOn(i) + iSoaCode * msPerFrame);
    if ~isempty(t), td.decMaskOn(i)   = t; end
    t  	= iTime((iCode >= iBaseCode + eDecFixOff) & (iCode < iBaseCode + eDecFixOff + 100));
    if ~isempty(t), td.decFixOff(i)   = t; end
    td.decResponseCueOn(i) = td.decFixOff(i);
    
    
    % Bet Stage
    t  	= iTime((iCode >= iBaseCode + eBetFixOn) & (iCode < iBaseCode + eBetFixOn + 100));
    if ~isempty(t), td.betFixOn(i)   = t; end
    t  	= iTime((iCode >= iBaseCode + eBetEyeInFixWindow) & (iCode < iBaseCode + eBetEyeInFixWindow + 100));
    if ~isempty(t), td.betFixWindowEntered(i)   = t; end
    t  	= iTime((iCode >= iBaseCode + eBetTargOn) & (iCode < iBaseCode + eBetTargOn + 100));
    if ~isempty(t), td.betTargOn(i)   = t; end
    td.betResponseCueOn(i) = td.betTargOn(i);
    
    % Reward
    t   = iTime((iCode == eRewardOn) | (iCode == eRewardClick));
    if ~isempty(t), td.rewardOn{i}   = t; end
    t   = iTime(iCode == eTimeoutOn);
    if ~isempty(t), td.timeoutOn(i)   = t; end
    t   = iTime(iCode == eTimeoutOff);
    if ~isempty(t), td.timeoutOff(i)   = t; end
    
    
    
    
    
    
    
    
    
    
    % ==================================================
    % Eye movement data
    % ==================================================
    td.eyeX{i} = flipud(allh(i,find(allh(i,:), 1) : end)');
    td.eyeY{i} = flipud(allv(i,find(allh(i,:), 1) : end)');
    
    
    
    
    if i == 50
        disp('huh')
    end
    % Use saccade algorithm to determine beginnings and ends of saccades,
    % etc
    [saccStart, saccEnd, ederiv] = find_saccades(td.eyeX{i}, td.eyeY{i});
    
    td.saccBegin{i} = saccStart;
    td.saccEnd{i} = saccEnd;
    
    yDiff = td.eyeY{i}(saccEnd) - td.eyeY{i}(saccStart);
    xDiff = td.eyeX{i}(saccEnd) - td.eyeX{i}(saccStart);
    td.saccAmp{i} = sqrt((yDiff .* yDiff) + (xDiff .* xDiff) );
    td.saccAngle{i} = atan2(yDiff, xDiff) .* 180/pi;
    
    
    
    % Find the decision saccade among the trial's saccades:
    decSaccInd = find((td.saccBegin{i} > td.decResponseCueOn(i)) & (td.saccBegin{i} < min(td.betFixOn(i), td.decResponseCueOn(i) + maxRTAllowed))); % all saccades after decision saccade cue
    
    if ~isempty(decSaccInd)
        foundDecSacc = false;
        dSacc = decSaccInd(1);
        lastDecSacc = decSaccInd(end);
        while ~foundDecSacc && dSacc <= lastDecSacc
            % Does the saccade start in the fixation window and end outside the
            % fixation window (in a possible target/distractor window)?
            startAtFix = in_window(td.eyeX{i}(td.saccBegin{i}(dSacc)), td.eyeY{i}(td.saccBegin{i}(dSacc)), 0, 0, td.decFixWindow(i));
            if ~startAtFix
                dSacc = dSacc + 1;
                continue
            end
            
            % Go through possible landing places of the saccade until we
            % find which location (if any) it landed
            endAtTarg = false;%zeros(length(maskAngles), 1);
            k = 0;
            while ~endAtTarg && k < length(maskAngles)
                k = k + 1;
                endAtTarg = in_window(td.eyeX{i}(td.saccEnd{i}(dSacc)), td.eyeY{i}(td.saccEnd{i}(dSacc)), targAmp, maskAngles(k), td.decTargWindow(i));
            end
            if ~endAtTarg
                dSacc = dSacc + 1;
                continue
            else
                foundDecSacc = true;
                td.decSaccIndex(i) = dSacc;
                td.decSaccAngle(i) = maskAngles(k);
                td.decResponseOnset(i)        = td.saccBegin{i}(dSacc);
                td.decRT(i)        = td.saccBegin{i}(dSacc) - td.decResponseCueOn(i);
            end
        end
    end % while ~foundDecSacc
    
    
    
    % Find the bet saccade among the trial's saccades:
    betSaccInd = find(td.saccBegin{i} > td.betResponseCueOn(i)); % all saccades after decision saccade cue
    
    if ~isempty(betSaccInd)
        foundBetSacc = false;
        bSacc = betSaccInd(1);
        lastBetSacc = betSaccInd(end);
        while ~foundBetSacc && bSacc <= lastBetSacc
            
            % Does the saccade start in the fixation window and end outside the
            % fixation window (in a possible target/distractor window)?
            startAtFix = in_window(td.eyeX{i}(td.saccBegin{i}(bSacc)), td.eyeY{i}(td.saccBegin{i}(bSacc)), 0, 0, td.decFixWindow(i));
            if ~startAtFix
                bSacc = bSacc + 1;
                continue
            end
            
            % Go through possible landing places of the saccade until we
            % find which location (if any) it landed
            endAtTarg = false;%zeros(length(betAngles), 1);
            k = 0;
            while ~endAtTarg && k < length(betAngles)
                k = k + 1;
                endAtTarg = in_window(td.eyeX{i}(td.saccEnd{i}(bSacc)), td.eyeY{i}(td.saccEnd{i}(bSacc)), targAmp, betAngles(k), td.betTargWindow(i));
            end
            if ~endAtTarg
                bSacc = bSacc + 1;
                continue
            else
                foundBetSacc = true;
                td.betSaccIndex(i) = bSacc;
                td.betSaccAngle(i) = betAngles(k);
                td.betResponseOnset(i)        = td.saccBegin{i}(bSacc);
                td.betRT(i)        = td.saccBegin{i}(bSacc) - td.betResponseCueOn(i);
            end
            
        end % while ~foundBetSacc
    end
    
    
    
    
    
    % ==================================================
    % Use to visual saccade data for debugging
    % ==================================================
    startDebugTrial = 1;
    if (DEBUGSACC && i >= startDebugTrial) %|| td.decRT(i) > 1000
        fprintf('\n \n Trial %d\n ', i)
        fprintf('target: %.1f \t highbet: %.1f \n', td.decTargAngle(i), td.betHighAngle(i))
        fprintf('targsac: %.1f \t betsac: %.1f \n', td.decSaccAngle(i), td.betSaccAngle(i))
        fprintf('decision: %s \t bet: %s \n', td.decOutcome{i}, td.betOutcome{i})
        
        figure(4)
        clf
        hold on
        set(gca,'ylim',[-20 20],'xlim',[-20 20])
        plot(td.eyeX{i},td.eyeY{i})
        [td.saccBegin{i},td.saccEnd{i},td.eyeX{i}(td.saccBegin{i}),td.eyeX{i}(td.saccEnd{i}),td.eyeY{i}(td.saccBegin{i}),td.eyeY{i}(td.saccEnd{i})]
        %             plot(td.eyeX{i}(td.saccBegin{i}),td.eyeY{i}(td.saccBegin{i}),'og','markersize',10,'linewidth',2)
        %             plot(td.eyeX{i}(td.saccEnd{i}),td.eyeY{i}(td.saccEnd{i}),'xr','markersize',10,'linewidth',2)
        if ~isnan(td.decSaccIndex(i))
            plot(td.eyeX{i}(td.saccBegin{i}(td.decSaccIndex(i))),td.eyeY{i}(td.saccBegin{i}(td.decSaccIndex(i))),'.g','markersize',50)
            plot(td.eyeX{i}(td.saccEnd{i}(td.decSaccIndex(i))),td.eyeY{i}(td.saccEnd{i}(td.decSaccIndex(i))),'.r','markersize',50)
        end
        if ~isnan(td.betSaccIndex(i))
            plot(td.eyeX{i}(td.saccBegin{i}(td.betSaccIndex(i))),td.eyeY{i}(td.saccBegin{i}(td.betSaccIndex(i))),'.b','markersize',30)
            plot(td.eyeX{i}(td.saccEnd{i}(td.betSaccIndex(i))),td.eyeY{i}(td.saccEnd{i}(td.betSaccIndex(i))),'.k','markersize',30)
        end
        
        % Draw the target and fixation accept boxes
        rectangle('position',[-td.decFixWindow(i)/2, -td.decFixWindow(i)/2, td.decFixWindow(i), td.decFixWindow(i)])
        for r = 1 : length(td.decMaskAngle(i,:));
            %             xMax = td.decTargAmp(i) * cosd(td.decMaskAngle(i,r)) + td.decTargWindow(i)/2;
            xMin = td.decTargAmp(i) * cosd(td.decMaskAngle(i,r)) - td.decTargWindow(i)/2;
            %             yMax = td.decTargAmp(i) * sind(td.decMaskAngle(i,r)) + td.decTargWindow(i)/2;
            yMin = td.decTargAmp(i) * sind(td.decMaskAngle(i,r)) - td.decTargWindow(i)/2;
            rectangle('position',[xMin, yMin, td.decTargWindow(i), td.decTargWindow(i)])
        end
        
        for s = 1 : length(td.saccBegin{i})
            plot(td.eyeX{i}(td.saccBegin{i}(s)),td.eyeY{i}(td.saccBegin{i}(s)),'og','markersize',10,'linewidth',2)
            plot(td.eyeX{i}(td.saccEnd{i}(s)),td.eyeY{i}(td.saccEnd{i}(s)),'xr','markersize',10,'linewidth',2)
            pause
        end
        %         pause
    end
    
    
    % ==================================================
    % SPIKE data
    % ==================================================
    td.spikeData{i} = fliplr(allspk(i,find(allspk(i,:), 1) : end))';
end
% td.saccDuration = cellfun(@(x,y) x - y, td.saccEnd, td.saccBegin, 'uni', false);




% Calculate variables based on other variables
td.trialOnset       = rData.allstart' - rData.allstart(1);
td.trialDuration    = [diff(rData.allstart)'; rData.alltimes(end,end)];
td.soa              = td.decMaskOn - td.decTargOn;

SessionData.soaArray = unique(td.soa(~isnan(td.soa)));




% oldTrialData = dataset(...
%     {trialOutcome,                      'trialOutcome'},...
%     {retroProFlag,                      'retroProFlag'},...
%     {round(trialOnset*1000),            'trialOnset'},...
%     {round(trialDuration*1000),         'trialDuration'},...
%     {round(abortOnset*1000),            'abortOnset'},...
%     {round(perceptionFixSpotOnset*1000), 'perceptionFixSpotOnset'},...
%     {round(perceptionFixationOnset*1000),'perceptionFixationOnset'},...
%     {round(perceptionFixSpotDuration*1000),'perceptionFixSpotDuration'},...
%     {round(perceptionTargOnset*1000),  	'perceptionTargOnset'},...
%     {round(soa*1000),                   'soa'},...
%     {round(perceptionMaskOnset*1000), 	'perceptionMaskOnset'},...
%     {round(retroProCueOnset*1000),     	'retroProCueOnset'},...
%     {round(betFixSpotOnset*1000),      	'betFixSpotOnset'},...
%     {round(betFixationOnset*1000),    	'betFixationOnset'},...
%     {round(preBetFixDuration*1000),    	'preBetFixDuration'},...
%     {round(betTargetOnset*1000),       	'betTargetOnset'},...
%     {round(betResponseOnset*1000),     	'betResponseOnset'},...
%     {round(decisionFixSpotOnset*1000), 	'decisionFixSpotOnset'},...
%     {round(decisionFixationOnset*1000),	'decisionFixationOnset'},...
%     {round(decisionFixSpotDuration*1000),'decisionFixSpotDuration'},...
%     {round(decisionMaskOnset*1000),   	'decisionMaskOnset'},...
%     {round(decisionResponseOnset*1000),	'decisionResponseOnset'},...
%     {round(feedbackOnset*1000),         'feedbackOnset'},...
%     {round(feedbackDuration*1000),      'feedbackDuration'},...
%     {targetIndex,                       'targetIndex'},...
%     {decisionIndex,                     'decisionIndex'},...
%     {betIndex,                          'betIndex'},...
%     {fixationAngle,                     'fixationAngle'},...
%     {targetAngle,                       'targetAngle'},...
%     {highBetAngle,                      'highBetAngle'},...
%     {lowBetAngle,                       'lowBetAngle'},...
%     {fixationSize,          'fixationSize'},...
%     {targetSize,            'targetSize'},...
%     {maskSize,              'maskSize'},...
%     {betSize,               'betSize'},...
%     {perceptionFixColor,   	'perceptionFixColor'},...
%     {decisionFixColor,   	'decisionFixColor'},...
%     {targetColor,           'targetColor'},...
%     {maskColor,             'maskColor'},...
%     {betFixColor,           'betFixColor'},...
%     {highBetColor,          'highBetColor'},...
%     {lowBetColor,           'lowBetColor'});

%%

% eyeSampleHz = 1000;
% tda = saccade_data(td, eyeSampleHz);



trialData = td;

save(fullfile(saveDir,saveID), 'trialData', 'SessionData')

function inWindow = in_window(eyeX, eyeY, amp, angle, windowWidth)
xMax = amp * cosd(angle) + windowWidth/2;
xMin = amp * cosd(angle) - windowWidth/2;
yMax = amp * sind(angle) + windowWidth/2;
yMin = amp * sind(angle) - windowWidth/2;
if eyeX <= xMax &&...
        eyeX >= xMin &&...
        eyeY <= yMax &&...
        eyeY >= yMin
    inWindow = true;
else
    inWindow = false;
end
return
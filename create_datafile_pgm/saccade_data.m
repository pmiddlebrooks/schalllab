    function trialData = saccade_data(trialData, taskID, eyeSampleHz)

if nargin < 3
    eyeSampleHz = 1000; % Default eyetracker rate of 1000 hz
end

% Need to multiply ms time stamps by this number to match the sampling
% rate, and divide all rates by it.
MS_PER_SAMPLE = 1000 / eyeSampleHz;

% Set velocity threshold (degrees per second) (may want to make input)
velocityThreshold = 30;
%physiologically impossible saccade cutoff
velocityCutoff = 1000;

if strcmp(taskID, 'maskbet')
    maskAngles = unique(trialData.decTargAngle(~isnan(trialData.decTargAngle)));
    betAngles = unique(trialData.betHighAngle(~isnan(trialData.betHighAngle)));
end




% Get maximum values for preallocating matrices and other operations
% NOTE; if trial lengths are too long here maxSamples will be too large.
% This is a problem b/c then we will preallocate all of our output to be
% huge matrices and it will lead to out of memory issues.  But trial length
% should be pretty short and constant right?  So who cares?  Well if you
% end up pausing the task during data collection you will have one really
% long trial which screws everything up.  To get around that, we check for
% outliers to the trial length and remove them.  Then we limit how long
% Saccade detection will proceed for in the ITI later in the code if we
% find outliers.
nTrial       = size(trialData.trialOutcome, 1);
trialOnset = trialData.trialOnset;
if nTrial > 1
    outlierLimit    = mean(diff(trialOnset)) + (std(diff(trialOnset)) * 5);  % Use 5 standard deviations as a cutoff
    samples         = (diff(trialOnset));
    samples(samples > outlierLimit) = [];
    maxSamples      = max(samples);
else
    maxSamples      = Eot_;
    trialOnset     = 1;  % If inputing
end









%       Step through each trial and process saccades:
% ----------------------------------------------------------------

% Preallocate matrices with all values set to NaN
saccBegin       = num2cell(nan(nTrial, maxSamples), 2);
saccEnd         = num2cell(nan(nTrial, maxSamples), 2);
saccAmp         = num2cell(nan(nTrial, maxSamples), 2);
saccAngle       = num2cell(nan(nTrial, maxSamples), 2);
saccSpeed       = num2cell(nan(nTrial, maxSamples), 2);
saccDuration  	= num2cell(nan(nTrial, maxSamples), 2);
saccBlinks  	= num2cell(nan(nTrial, maxSamples), 2);




for i = 1 : nTrial
    
    % Calculate velocity of eye traces in each dimension (horizontal and
    % vertical components)
    deltaEyeX = [0 ; diff(trialData.eyeX{i})] ./ MS_PER_SAMPLE;
    deltaEyeY = [0 ; diff(trialData.eyeY{i})] ./ MS_PER_SAMPLE;
    
    % Calculate velocity using pythagorean formula on the horizontal and
    % vertical velocity vectors
    
    velocityXY = sqrt((deltaEyeX.^2) + (deltaEyeY.^2));
    clear deltaEyeX deltaEyeY %conserve memory for low powered machines
    
    % Convert velocity from degrees/ms to degrees/second
    velocityXY = velocityXY * 1000;
    
    % Determine where velocity is greater than or equal to threshold
    saccadeTime = velocityXY >= velocityThreshold;
    
    
    % Get threshold test
    saccStartsEnds = [0 ; diff(saccadeTime)];
    
    
    
    % Identify beginning of saccade
    iSaccBegins = find(saccStartsEnds == 1)' * MS_PER_SAMPLE;
    % Identify end of saccade
    iSaccEnds = (find(saccStartsEnds == -1) - 1)' * MS_PER_SAMPLE;
    
    
    
    
    
    % If no saccades store a NaN
    if isempty(iSaccBegins) || isempty(iSaccEnds)
        iSaccBegins = nan;
        iSaccEnds = nan;
        
        % If we have saccades...
    elseif ~isempty(iSaccBegins) && ~isempty(iSaccEnds)
        % ...do the starts and ends of the saccades match each other?
        if iSaccEnds(1) < iSaccBegins(1) % then we matched a saccade end from the last trial with a saccade begin from this trial
            iSaccEnds(1) = []; % delete it b/c saccades are classified based on when they start not end (see below)
        end
        if length(iSaccBegins) > length(iSaccEnds) % then a saccade was in flight during the end of the current trial
            %         % If a multiple trials were sent to saccade_data.m (assuming a single trial
            %         % takes less than 20 seconds, we can look ahead to the next trial for the
            %         % saccade end....
            %         if length(AllStarts) > 20000
            %             missedEndIndex = find(AllEnds >= nextTrialStart,1,'first'); % get the index of the end from the next trial
            %             iMissedEnd = AllEnds(missedEndIndex); % and then get the saccade end time
            %             iSaccEnds(end + 1) = iMissedEnd; % NOTE: the next trial should fall into the conditional statement above this one.
            %             % ...but if a single trial was sent, to saccade_data.m, we can't look into
            %             % the next trial and thus can't figure out when the saccade
            %             % ended. Just cutoff the last saccade start then
            %         else
            iSaccBegins(end) = [];
            
            %         end
        end
    end
    
    % Preallocate matrices with all values set to NaN
    iSaccAmps = nan(1, length(iSaccBegins));
    iSaccDirs = nan(1, length(iSaccBegins));
    iSaccDurs = nan(1, length(iSaccBegins));
    iSaccVels = nan(1, length(iSaccBegins));
    
    %     iTrial
    %     iSaccBegins
    %     iSaccEnds
    
    % Compute saccade dynamics
    for jSaccade = 1 : sum(~isnan(iSaccBegins))
        jBegin = iSaccBegins(jSaccade) / MS_PER_SAMPLE;
        jEnd   = iSaccEnds(jSaccade) / MS_PER_SAMPLE;
%         deltaX    = trialData.eyeX{i}(jBegin) - trialData.eyeX{i}(jEnd);
%         deltaY    = trialData.eyeY{i}(jBegin) - trialData.eyeY{i}(jEnd);
        deltaX    = trialData.eyeX{i}(jEnd) - trialData.eyeX{i}(jBegin);
        deltaY    = trialData.eyeY{i}(jEnd) - trialData.eyeY{i}(jBegin);
        %         deltaX    = eyeX(jBegin) - eyeX(jEnd);
        %         deltaY    = eyeY(jBegin) - eyeY(jEnd);
        
        % Compute saccade amplitude
        iSaccAmps(1, jSaccade) = sqrt(deltaX^2 + deltaY^2);
        
        % Compute saccade direction
        iSaccDirs(1, jSaccade) = mod((180/pi * atan2(deltaY, deltaX)), 360);
        
        % Compute saccade duration
        iSaccDurs(1, jSaccade) = jEnd - jBegin;
        
        % Compute saccade velocity
        iSaccVels(1, jSaccade) = max(velocityXY(jBegin : jEnd));
        
    end
    
    % Use a few common sense checks to remove physiologically impossible saccades
    iSaccsNBlinks = iSaccBegins;
    iSaccBegins(iSaccVels > velocityCutoff) = [];
    iSaccEnds  (iSaccVels > velocityCutoff) = [];
    iSaccAmps  (iSaccVels > velocityCutoff) = [];
    iSaccDirs  (iSaccVels > velocityCutoff) = [];
    iSaccDurs  (iSaccVels > velocityCutoff) = [];
    iSaccVels  (iSaccVels > velocityCutoff) = [];
    
    iSaccBegins(iSaccAmps == 0) = [];
    iSaccEnds  (iSaccAmps == 0) = [];
    iSaccDirs  (iSaccAmps == 0) = [];
    iSaccDurs  (iSaccAmps == 0) = [];
    iSaccVels  (iSaccAmps == 0) = [];
    iSaccAmps  (iSaccAmps == 0) = [];
    
    % Make sure we don't have too many saccades (see comment about calculating
    % maxSamples above)
    if length(iSaccBegins) > maxSamples
        
        iSaccBegins(maxSamples+1 : end) = [];
        iSaccEnds  (maxSamples+1 : end) = [];
        iSaccDirs  (maxSamples+1 : end) = [];
        iSaccDurs  (maxSamples+1 : end) = [];
        iSaccVels  (maxSamples+1 : end) = [];
        iSaccAmps  (maxSamples+1 : end) = [];
        
    end
    
    
    
    
    
    
    
    
    % Insert into matrices for whole recording session
    saccBegin{i}      = iSaccBegins;
    saccEnd{i}     = iSaccEnds;
    saccAmp{i}    = iSaccAmps;
    saccAngle{i}    = iSaccDirs;
    saccSpeed{i}     = iSaccVels;
    saccDuration{i}     = iSaccDurs;
    saccBlinks{i} = iSaccsNBlinks;
    
end


trialData.saccBegin         = saccBegin;
trialData.saccAmp         = saccAmp;
trialData.saccAngle         = saccAngle;
trialData.saccDuration         = saccDuration;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         &
% FIND THE SACCADE OF INTEREST ON ALL TRIALS: This is task-dependent(see below)   &
%                                                                         &
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saccToTargIndex is the first saccade which begins after Target_ and
% before Correct_ or EOT_.  Note that I have imposed the additional
% criteria that saccades must be at least 5 degrees in amplitude for them to
% count, so we may have saccade events which do not actually determine
% trial outcome and are therefore not classified as saccToTargIndex. This
% is because the dist of error sacccade amps is bimodal with some small
% saccades happening around fixation.  These fixational saccades tend to
% have longer amplitudes and probably reflect the fact that it is hard to
% keep the eye immobilized for 2 X the allowed time for a saccade.  If they
% are left in error RTs have a long tail which violates assumption of
% independance in race model.

if ~strcmp(taskID, 'maskbet')
    saccToTargIndex         = nan(nTrial, 1);
    Sacc_of_interest_end    = nan(nTrial, 1);
    %figure out saccade starting and ending positions
    % sacc_pos                = nan(nTrial, 4);
    % X_starts = 1;
    % X_ends   = 2;
    % Y_starts = 3;
    % Y_ends   = 4;
    
    for i = 1 : nTrial
        if ~strcmp(trialData.trialOutcome{i}, 'stopCorrect') && ...
                ~strcmp(trialData.trialOutcome{i}, 'fixationAbort') && ...
                ~strcmp(trialData.trialOutcome{i}, 'noFixation')
            % nSacc = min(length(saccBegin{iTrial}), length(saccEnd{iTrial}));
            nSacc = length(saccBegin{i});
            
            afterTarg   = repmat(trialData.responseCueOn(i), 1, nSacc);
            beforeReward     = repmat(trialData.trialDuration(i), 1, nSacc);
            % size(saccBegin{iTrial})
            % size(afterTarg)
            % size(saccAmp{iTrial})
            
            possibleSacc  = (saccBegin{i} - afterTarg) > 0   &...
                (saccBegin{i} - beforeReward) <= 0 &...
                saccAmp{i} > 2;
            
            
            foundIt = 0;
            jSacc = find(possibleSacc, 1, 'first');
            if ~isempty(jSacc) && ~isnan(jSacc)
                while ~foundIt && jSacc <= length(find(~isnan(saccEnd{i})))
                    % If started in fixation window...And ended outside the fixation window...
                                        
                    if (trialData.eyeX{i}(saccBegin{i}(jSacc)) <= trialData.fixWindow(i)/2 && ...
                            trialData.eyeX{i}(saccBegin{i}(jSacc)) >= -trialData.fixWindow(i)/2) && ...
                            (trialData.eyeY{i}(saccBegin{i}(jSacc)) <= trialData.fixWindow(i)/2 && ...
                            trialData.eyeY{i}(saccBegin{i}(jSacc)) >= -trialData.fixWindow(i)/2) && ...
                            (trialData.eyeX{i}(saccEnd{i}(jSacc)) < -trialData.fixWindow(i)/2 || ...
                            trialData.eyeX{i}(saccEnd{i}(jSacc)) > trialData.fixWindow(i)/2 || ...
                            trialData.eyeY{i}(saccEnd{i}(jSacc)) < -trialData.fixWindow(i)/2 || ...
                            trialData.eyeY{i}(saccEnd{i}(jSacc)) > trialData.fixWindow(i)/2)% && ...
                        
                        foundIt = 1;
                        %     if ~isempty(curr_Sacc)
                        saccToTargIndex(i)  = jSacc;
                        Sacc_of_interest_end(i) = saccEnd{i}(jSacc);
                        
                        %                     sacc_pos(iTrial,X_starts) = trialData.eyeX{iTrial}(saccBegin{iTrial}(jSacc) / MS_PER_SAMPLE);
                        %                     sacc_pos(iTrial,X_ends)   = trialData.eyeX{iTrial}(Sacc_of_interest_end(iTrial) / MS_PER_SAMPLE);
                        %                     sacc_pos(iTrial,Y_starts) = trialData.eyeY{iTrial}(saccBegin{iTrial}(jSacc) / MS_PER_SAMPLE);
                        %                     sacc_pos(iTrial,Y_ends)   = trialData.eyeY{iTrial}(Sacc_of_interest_end(iTrial) / MS_PER_SAMPLE);
                        %     end
                    else
                        jSacc = jSacc + 1;
                    end
                end
            end
        end
    end
    
    
    
    trialData.saccToTargIndex         = saccToTargIndex;
    
    trialWithResponse       = ~isnan(saccToTargIndex);
    responseOnset = nan(nTrial, 1);
    responseOnset(trialWithResponse) = cellfun(@(x, y) x(y), saccBegin(trialWithResponse), num2cell(saccToTargIndex(trialWithResponse), 2));
    trialData.responseOnset         = responseOnset;
    
    
    
    
elseif strcmp(taskID, 'maskbet')
    trialData.decSaccIndex         = nan(nTrial, 1);
    trialData.decSaccAngle         = nan(nTrial, 1);
    trialData.decResponseOnset   	= nan(nTrial, 1);
    trialData.decRT                 = nan(nTrial, 1);
    %     decSacc_of_interest_end    = nan(nTrial, 1);
    trialData.betSaccIndex         = nan(nTrial, 1);
    trialData.betSaccAngle         = nan(nTrial, 1);
    trialData.betResponseOnset      = nan(nTrial, 1);
    trialData.betRT                 = nan(nTrial, 1);
    %     betSacc_of_interest_end    = nan(nTrial, 1);
    
    for i = 1 : nTrial
        
        % nSacc = min(length(saccBegin{iTrial}), length(saccEnd{iTrial}));
        nSacc = length(saccBegin{i});
        
        
        
        
        % ******************************************
        % DECISION STAGE saccade
        afterTarg   = repmat(trialData.decResponseCueOn(i), 1, nSacc);
        beforeBet     = repmat(trialData.decResponseCueOn(i)+2000, 1, nSacc);  % Allowing 2000 ms to make a decsion saccade (more than what is allowed in task, so should be enough)
        
        possibleSacc  = (saccBegin{i} - afterTarg) > 0   &...
            (saccBegin{i} - beforeBet) <= 0 &...
            saccAmp{i} > 2;
        
        
        foundDecSacc = false;
        dSacc = find(possibleSacc, 1, 'first');
        if ~isempty(dSacc) && ~isnan(dSacc)
            while ~foundDecSacc && dSacc <= length(possibleSacc)
                % If started in fixation window...And ended outside the fixation window...
                
                startAtFix = in_window(...
                    trialData.eyeX{i}(saccBegin{i}(dSacc)), ...
                    trialData.eyeY{i}(saccBegin{i}(dSacc)), ...
                    trialData.decFixAmp(i), ...
                    trialData.decFixAngle(i), ...
                    trialData.decFixWindow(i));
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
                    endAtTarg = in_window(...
                        trialData.eyeX{i}(saccEnd{i}(dSacc)), ...
                        trialData.eyeY{i}(saccEnd{i}(dSacc)), ...
                        trialData.decTargAmp(i), ...
                        maskAngles(k), ...
                        trialData.decTargWindow(i));
                end
                if ~endAtTarg
                    dSacc = dSacc + 1;
                    continue
                else
                    foundDecSacc = true;
                    trialData.decSaccIndex(i) = dSacc;
                    trialData.decSaccAngle(i) = maskAngles(k);
                    trialData.decResponseOnset(i) = trialData.saccBegin{i}(dSacc);
                    trialData.decRT(i)        = trialData.saccBegin{i}(dSacc) - trialData.decResponseCueOn(i);
                end
            end
        end
        
        
        
        
        % ******************************************
        % BET STAGE SACCADE
        afterTarg   = repmat(trialData.betResponseCueOn(i), 1, nSacc);
        beforeReward     = repmat(trialData.trialDuration(i), 1, nSacc);
        
        possibleSacc  = (saccBegin{i} - afterTarg) > 0   &...
            (saccBegin{i} - beforeReward) <= 0 &...
            saccAmp{i} > 2;
        
        
        foundBetSacc = false;
        bSacc = find(possibleSacc, 1, 'first');
        if ~isempty(bSacc) && ~isnan(bSacc)
            while ~foundBetSacc && bSacc <= length(possibleSacc)
                % If started in fixation window...And ended outside the fixation window...
                
                startAtFix = in_window(...
                    trialData.eyeX{i}(saccBegin{i}(bSacc)), ...
                    trialData.eyeY{i}(saccBegin{i}(bSacc)), ...
                    trialData.decFixAmp(i), ...
                    trialData.decFixAngle(i), ...
                    trialData.decFixWindow(i));
                if ~startAtFix
                    bSacc = bSacc + 1;
                    continue
                end
                
                % Go through possible landing places of the saccade until we
                % find which location (if any) it landed
                endAtTarg = false;%zeros(length(maskAngles), 1);
                k = 0;
                while ~endAtTarg && k < length(betAngles)
                    k = k + 1;
                    endAtTarg = in_window(...
                        trialData.eyeX{i}(saccEnd{i}(bSacc)), ...
                        trialData.eyeY{i}(saccEnd{i}(bSacc)), ...
                        trialData.decTargAmp(i), ...
                        betAngles(k), ...
                        trialData.betTargWindow(i));
                end
                if ~endAtTarg
                    bSacc = bSacc + 1;
                    continue
                else
                    foundBetSacc = true;
                    trialData.betSaccIndex(i) = bSacc;
                    trialData.betSaccAngle(i) = betAngles(k);
                    trialData.betResponseOnset(i) = trialData.saccBegin{i}(bSacc);
                    trialData.betRT(i)        = trialData.saccBegin{i}(bSacc) - trialData.betResponseCueOn(i);
                end
            end
        end
        
        
        trialData.saccBegin         = saccBegin;
        trialData.saccAmp         = saccAmp;
        trialData.saccAngle         = saccAngle;
        trialData.saccDuration         = saccDuration;
        %         trialData.decSaccIndex         = decSaccIndex;
        %         trialData.betSaccIndex         = betSaccIndex;
        
        %         trialWithDecResponse       = ~isnan(decSaccIndex);
        %         decResponseOnset = nan(nTrial, 1);
        %         decResponseOnset(trialWithDecResponse) = cellfun(@(x, y) x(y), saccBegin(trialWithDecResponse), num2cell(decSaccIndex(trialWithDecResponse), 2));
        %         trialData.decResponseOnset         = decResponseOnset;
        %
        %         trialWithBetResponse       = ~isnan(betSaccIndex);
        %         betResponseOnset = nan(nTrial, 1);
        %         betResponseOnset(trialWithBetResponse) = cellfun(@(x, y) x(y), saccBegin(trialWithBetResponse), num2cell(betSaccIndex(trialWithBetResponse), 2));
        %         trialData.betResponseOnset         = betResponseOnset;
        
    end
    
    
    
    
    
    
    % %use same algorithm to find the SecondSacc for saving to file
    % [rows cols] = size(saccBlinks);
    % after_Sacc_of_interest = repmat(saccToTargIndex(:,1),1,cols);
    % possibleSacc  = (saccBlinks - after_Sacc_of_interest) > 0;
    % SecondSacc(1:rows,1:2)   = nan;
    % for iTrial = 1 : nTrial
    %     curr_Sacc = find(possibleSacc(iTrial,:),1,'first');
    %     if ~isempty(curr_Sacc)
    %         SecondSacc(iTrial,2) = curr_Sacc;
    %         SecondSacc(iTrial,1) = saccBlinks(iTrial,curr_Sacc);
    %     end
    % end
    
    
    
    
    
    
end










%
% % subfunction
%
% function [iSaccBegins,...
%     iSaccEnds,...
%     iSaccAmps,...
%     iSaccDirs,...
%     iSaccVels,...
%     iSaccDurs,...
%     iSaccsNBlinks] = processTrial(iTrialStart,...
%     nextTrialStart,...
%     EyeX_,...
%     EyeY_,...
%     velocityXY,...
%     AllStarts,...
%     AllEnds,...
%     maxSamples,...
%     MS_PER_SAMPLE)
%
% %physiologically impossible saccade cutoff
% velocityCutoff = 1000;
%
% % Which saccades happened on this trial?
% iSaccBegins = AllStarts(AllStarts >= iTrialStart &...
%     AllStarts < nextTrialStart);
% iSaccEnds = AllEnds(AllEnds >= iTrialStart &...
%     AllEnds < nextTrialStart);
%
% % If we have saccades...
% if ~isempty(iSaccBegins) && ~isempty(iSaccEnds)
%     % ...do the starts and ends of the saccades match each other?
%     if iSaccEnds(1) < iSaccBegins(1) % then we matched a saccade end from the last trial with a saccade begin from this trial
%         iSaccEnds(1) = []; % delete it b/c saccades are classified based on when they start not end (see below)
%     end
%
%     if length(iSaccBegins) > length(iSaccEnds) % then a saccade was in flight during the end of the current trial
%         % If a multiple trials were sent to saccade_data.m (assuming a single trial
%         % takes less than 20 seconds, we can look ahead to the next trial for the
%         % saccade end....
%         if length(AllStarts) > 20000
%             missedEndIndex = find(AllEnds >= nextTrialStart,1,'first'); % get the index of the end from the next trial
%             iMissedEnd = AllEnds(missedEndIndex); % and then get the saccade end time
%             iSaccEnds(end + 1) = iMissedEnd; % NOTE: the next trial should fall into the conditional statement above this one.
%             % ...but if a single trial was sent, to saccade_data.m, we can't look into
%             % the next trial and thus can't figure out when the saccade
%             % ended. Just cutoff the last saccade start then
%         else
%             iSaccBegins(end) = [];
%         end
%     end
%
%     % if this ever happens my logic is flawed.
%     if length(~isnan(iSaccBegins)) ~= length(~isnan(iSaccEnds))
%         warning(sprintf('David screwed up!\nSaccade begin & end mismatch!\nPANIC! PANIC!'));
%     end
%
%     % If no saccades store a NaN
% elseif isempty(iSaccBegins) || isempty(iSaccEnds)
%     iSaccBegins = NaN;
%     iSaccEnds = NaN;
% end
%
%
%
% % Preallocate matrices with all values set to NaN
% clear iSaccAmps iSaccDirs iSaccDurs iSaccVels
% iSaccAmps(1, 1:length(iSaccBegins)) = NaN;
% iSaccDirs(1, 1:length(iSaccBegins)) = NaN;
% iSaccDurs(1, 1:length(iSaccBegins)) = NaN;
% iSaccVels(1, 1:length(iSaccBegins)) = NaN;
%
% % Compute saccade dynamics
% for jSaccade = 1 : sum(~isnan(iSaccBegins))
%     jBegin = iSaccBegins(jSaccade) / MS_PER_SAMPLE;
%     jEnd   = iSaccEnds(jSaccade) / MS_PER_SAMPLE;
%     deltaX    = EyeX_(jBegin) - EyeX_(jEnd);
%     deltaY    = EyeY_(jBegin) - EyeY_(jEnd);
%
%     % Compute saccade amplitude
%     iSaccAmps(1, jSaccade) = sqrt(deltaX^2 + deltaY^2);
%
%     % Compute saccade direction
%     iSaccDirs(1, jSaccade) = mod((180/pi * atan2(deltaY, deltaX)), 360);
%
%     % Compute saccade duration
%     iSaccDurs(1, jSaccade) = jEnd - jBegin;
%
%     % Compute saccade velocity
%     iSaccVels(1, jSaccade) = max(velocityXY(jBegin:jEnd));
%
% end
%
% % Use a few common sense checks to remove physiologically impossible saccades
% iSaccsNBlinks = iSaccBegins;
% iSaccBegins(iSaccVels > velocityCutoff) = [];
% iSaccEnds  (iSaccVels > velocityCutoff) = [];
% iSaccAmps  (iSaccVels > velocityCutoff) = [];
% iSaccDirs  (iSaccVels > velocityCutoff) = [];
% iSaccDurs  (iSaccVels > velocityCutoff) = [];
% iSaccVels  (iSaccVels > velocityCutoff) = [];
%
% iSaccBegins(iSaccAmps == 0) = [];
% iSaccEnds  (iSaccAmps == 0) = [];
% iSaccDirs  (iSaccAmps == 0) = [];
% iSaccDurs  (iSaccAmps == 0) = [];
% iSaccVels  (iSaccAmps == 0) = [];
% iSaccAmps  (iSaccAmps == 0) = [];
%
% % Make sure we don't have too many saccades (see comment about calculating
% % maxSamples above)
% if length(iSaccBegins) > maxSamples
%
%     iSaccBegins(maxSamples+1 : end) = [];
%     iSaccEnds  (maxSamples+1 : end) = [];
%     iSaccDirs  (maxSamples+1 : end) = [];
%     iSaccDurs  (maxSamples+1 : end) = [];
%     iSaccVels  (maxSamples+1 : end) = [];
%     iSaccAmps  (maxSamples+1 : end) = [];
%
% end
%
%
%
% end
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
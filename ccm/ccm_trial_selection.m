function  [trialList] = ccm_trial_selection(trialData, selectOpt)
% function  [trialList] = ccm_trial_selection(trialData, outcome, rightCheckerPct, ssd, targDir)
%
%
%
%
% [trialList] = ccm_trial_selection(trialData, selctOpt)
%
% Returns a list of the trial numbers with conditions specified in options
% structure. If called without any arguments, returns a default options structure.
% If options are input but one is not specified, it assumes default.
%
% Requires the following dataset VarNames: {'trialOutcome', 'targ1CheckerProp', 'ssd', 'targAngle', 'saccToTargIndex', 'saccAngle'};
%
% Possible conditions are (default listed first):
%     selectOpt.outcome  = array of strings indicating the outcomes to include:
%           {'collapse',
%           'valid'
%           'goCorrectTarget', 'goCorrectDistractor',
%           'goIncorrect',
%           'stopCorrect',
%           'stopIncorrectTarget', 'stopIncorrectDistractor',
%           'stopIncorrectTarget', 'stopIncorrectDistractor',
%           'targetHoldAbort', 'distractorHoldAbort',
%           'fixationAbort', 'saccadeAbort', 'checkerStimulusAbort'}
%     selectOpt.choiceAccuracy  = default is collapse across all choices. Segments
%           with respect to weather a choice was correct or error (does not count
%           if no choice was made.
%           options: 'collapse' (default), 'correct', 'error'.
%     selectOpt.rightCheckerPct  = range of checkerboard percentage of right target checkers:
%           {'collapse', 'right', 'left'
%           a double array containing the values, e.g. [40 50 60]
%     selectOpt.allowRtPreSsd    = whether to allow noncanceled stop trials with RTs before SSDs
%           true (default) or false
%     selectOpt.ssd    = range of SSDs to include in the trial list:
%           {'collapse', 'any', 'none', or
%           a double array containing the values, e.g. [43 86 129]
%     selectOpt.targDir  = the angle of the CORRECT TARGET
%           {'collapse','left','right','up','down','upLeft','upRight',downLeft','downRight'<list of possible angels, e.g. [45 -45]
%     selectOpt.responseDir  = the angle of target to which a response was made
%           {'collapse','left','right','up','down','upLeft','upRight',downLeft','downRight'<list of possible angels, e.g. [45 -45]





% If not input, return the default options structure
if nargin < 2
    selectOpt.outcome           = 'valid';
    selectOpt.choiceAccuracy   	= 'collapse';
    selectOpt.rightCheckerPct   = 'collapse';
    selectOpt.ssd               = 'any';
    selectOpt.allowRtPreSsd     = true;
    selectOpt.targDir           = 'collapse';
    selectOpt.responseDir       = 'collapse';
    selectOpt.doStops           = true;
    if nargin < 1
        trialList                   = selectOpt;
        return
    end
end


if strcmp(selectOpt.outcome, 'stopStop')
    error('stopStop is not valid. Use stopCorrect')
end
% if ~iscell(outcome)
%     error('ccm_trial_selection: need to enter the outcomArray as a cell, e.g. {''valid''}')
% end


nTrial = size(trialData.trialOutcome, 1);
trialLogical = ones(nTrial, 1);

trialData = cell_to_mat(trialData);


if isfield(selectOpt, 'choiceAccuracy')
    choiceLogical = zeros(nTrial, 1);
    switch selectOpt.choiceAccuracy
        case 'collapse'
            choiceArray = {};
            choiceLogical = ones(nTrial, 1);
            % do nothing
        case 'correct'
            choiceArray = {...
                'goCorrectTarget',...
                'stopIncorrectTarget',...
                'targetHoldAbort'};
        case 'error'
            choiceArray = {...
                'goCorrectDistractor',...
                'stopIncorrectDistractor',...
                'distractorHoldAbort'};
    end
    if ~isempty(choiceArray)
        for iChoiceIndex = 1 : length(choiceArray)
            iChoice = choiceArray{iChoiceIndex};
            
            choiceLogical = choiceLogical + strcmp(trialData.trialOutcome, iChoice);
        end
    end
    trialLogical = trialLogical & choiceLogical;
end





% Trials w.r.t. the outcome
if strcmp(selectOpt.outcome, 'valid')
    selectOpt.outcome = {...
        'goCorrectTarget', 'goCorrectDistractor', ...
        'targHoldAbort', 'distractorHoldAbort', ...
        'stopCorrect', ...
        'stopIncorrectTarget', 'stopIncorrectDistractor'};
    
end

if strcmp(selectOpt.outcome, 'collapse')
    outcomeLogical = ones(nTrial, 1);
else
    outcomeLogical = zeros(nTrial, 1);
    for iOutcomeIndex = 1 : length(selectOpt.outcome)
        iOutcome = selectOpt.outcome{iOutcomeIndex};
        
        outcomeLogical = outcomeLogical + strcmp(trialData.trialOutcome, iOutcome);
    end
end
trialLogical = trialLogical & outcomeLogical;




% Trials w.r.t. the checkerboard proportion of target1
% checkers
if strcmp(selectOpt.rightCheckerPct, 'collapse')
    % do nothing
else
    if strcmp(selectOpt.rightCheckerPct, 'right')
        target1Proportion = selectOpt.rightCheckerPct/100 > 50;
    elseif strcmp(selectOpt.rightCheckerPct, 'left')
        target1Proportion = selectOpt.rightCheckerPct/100 < 50;
    else
        target1Proportion = selectOpt.rightCheckerPct / 100;
    end
    trialLogical = trialLogical & ismember(trialData.targ1CheckerProp, target1Proportion);
end






% Trials w.r.t. response direction (or no response)
if isfield(selectOpt, 'responseDir') && ~strcmp(selectOpt.responseDir, 'collapse')
    
    if ismember('responseDirection', fieldnames(trialData))
        saccTrial = strcmp(trialData.responseDirection, selectOpt.responseDir);
        responseTrial = ~isnan(trialData.rt);
    elseif ~strcmp(selectOpt.responseDir ,'none')
        saccAngle = nan(nTrial, 1);
        responseTrial = ~isnan(trialData.saccToTargIndex);
        
        %    trialData.saccAngle(nanResp) = cellfun(@(x) [x 0], trialData.saccAngle(nanResp), 'uni', false);  % For trials without responses, need to pretend there is one to sort data (these trials won't be included in trialList)
        saccAngle(responseTrial) = cellfun(@(x,y) x(y), trialData.saccAngle(responseTrial), num2cell(trialData.saccToTargIndex(responseTrial)));
        %    saccAngle = cell2num(trialData.saccAngle);
        %     switch selectOpt.responseDir
        if strcmp(selectOpt.responseDir, 'collapse')
            % Do nothing
            saccTrial = ones(nTrial, 1);
        elseif strcmp(selectOpt.responseDir, 'right')
            saccTrial = ((saccAngle > 270) & (saccAngle <= 360)) | ...
                ((saccAngle >= 0) & (saccAngle < 90)) | ...
                ((saccAngle > -90) & (saccAngle < 0)) | ...
                ((saccAngle >= -360) & (saccAngle < -270));
        elseif strcmp(selectOpt.responseDir, 'left')
            saccTrial = ((saccAngle > 90) & (saccAngle <= 270)) | ...
                ((saccAngle < -90) & (saccAngle > -270));
        elseif strcmp(selectOpt.responseDir, 'leftUp')
            saccTrial = saccAngle == decAngleLeftUp;
        elseif strcmp(selectOpt.responseDir, 'leftDown')
            saccTrial = saccAngle == decAngleLeftDown;
        elseif strcmp(selectOpt.responseDir, 'rightUp')
            saccTrial = saccAngle == decAngleRightUp;
        elseif strcmp(selectOpt.responseDir, 'rightDown')
            saccTrial = saccAngle == decAngleRightDown;
        else
            saccTrial = ismember(saccAngle,selectOpt.responseDir);
        end
    elseif strcmp(selectOpt.responseDir ,'none')
        saccTrial = isnan(trialData.saccToTargIndex);
    end
    trialLogical = trialLogical & saccTrial;
    % For go trials and noncanceled stop trials, get rid of trials without
    % responses (for canceled stop trials all of them should be trials
    % without responses, so leave them in)
    if strcmp(selectOpt.ssd, 'none')
        trialLogical = trialLogical & responseTrial;
    end
end




% Trials w.r.t. the SSDs and whether to allow RTs that preceded SSD on noncanceled stop trials
if selectOpt.doStops
    if strcmp(selectOpt.ssd, 'none')
        % take any trials without a stop signal (nan values for selectOpt.ssd)
        trialLogical = trialLogical & isnan(trialData.ssd);
    else
        
    % Trials w.r.t. whether to allow RTs that preceded SSD on noncanceled stop
    % trials
    if ~selectOpt.allowRtPreSsd
        trialLogical = trialLogical & ~isnan(trialData.ssd)  & ~isnan(trialData.rt) & trialData.rt > trialData.ssd;
    end
        
    if strcmp(selectOpt.ssd, 'any')
        % Do nothing- might want trials without regard to stop/go
    elseif strcmp(selectOpt.ssd, 'collapse')
        ssd = min(trialData.ssd) : max(trialData.ssd);
        trialLogical = trialLogical & ismember(trialData.ssd, ssd);
    else
        %     ssd = [ssd - 1, ssd, ssd + 1];
        % For now, use range to within an extra frame refresh
        ssd = [selectOpt.ssd - 5 : selectOpt.ssd + 5];
        trialLogical = trialLogical & ismember(trialData.ssd, ssd);
    end
    end
end


% Trials w.r.t. target angle (useful for 50%
% checkerboard displays wherein a random target is assigned
% Get list(s) of trials w.r.t. the SSDs
targAngle = trialData.targAngle;
if strcmp(selectOpt.targDir, 'collapse')
    % Do nothing
    targTrial = ones(nTrial, 1);
elseif strcmp(selectOpt.targDir, 'right')
    targTrial = ((targAngle > 270) & (targAngle <= 360)) | ...
        ((targAngle >= 0) & (targAngle < 90)) | ...
        ((targAngle > -90) & (targAngle < 0)) | ...
        ((targAngle >= -360) & (targAngle < -270));
elseif strcmp(selectOpt.targDir, 'left')
    targTrial = ((targAngle > 90) & (targAngle <= 270)) | ...
        ((targAngle < -90) & (targAngle > -270));
elseif strcmp(selectOpt.targDir, 'leftUp')
    targTrial = targAngle == decAngleLeftUp;
elseif strcmp(selectOpt.targDir, 'leftDown')
    targTrial = targAngle == decAngleLeftDown;
elseif strcmp(selectOpt.targDir, 'rightUp')
    targTrial = targAngle == decAngleRightUp;
elseif strcmp(selectOpt.targDir, 'rightDown')
    targTrial = targAngle == decAngleRightDown;
else
    targTrial = ismember(targAngle,selectOpt.targDir);
end
trialLogical = trialLogical & targTrial;





% If there was supposed to be an RT but there wasn't (if there is a NaN for
% RT, get rid of those trials

% % RTs occur on the following trial types:
% if any(ismember(outcome, {...
%       'goCorrectTarget', ...
%       'goCorrectDistractor', ...
%       'goIncorrectTarget', ...
%       'goIncorrectDistractor', ...
%       'targetHoldAbort', ...
%       'distractorHoldAbort', ...
%       'stopIncorrectTarget', ...
%       'stopIncorrectDistractor'}))
%    rtTrial = strcmp(trialData.trialOutcome, 'goCorrectTarget') | ...
%       strcmp(trialData.trialOutcome, 'goCorrectDistractor') | ...
%       strcmp(trialData.trialOutcome, 'goIncorrectTarget') | ...
%       strcmp(trialData.trialOutcome, 'goIncorrectDistractor') | ...
%       strcmp(trialData.trialOutcome, 'targetHoldAbort') | ...
%       strcmp(trialData.trialOutcome, 'distractorHoldAbort') | ...
%       strcmp(trialData.trialOutcome, 'stopIncorrectTarget') | ...
%       strcmp(trialData.trialOutcome, 'stopIncorrectDistractor');
%    validRT = ~isnan(trialData.rt) & rtTrial;
%    trialLogical = trialLogical & validRT;
% end







trialList = find(trialLogical);


% toc
% disp('angle')



















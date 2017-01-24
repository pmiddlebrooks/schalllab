function  [trialList] = maskbet_trial_selection(trialData, selectOpt)%, decOutcome, betOutcome, soa, decSaccDir, decTargDir, betSaccDir, betHighDir)
% function  [trialList] = maskbet_trial_selection(trialData, selectOpt)

% [trialList] = maskbet_trial_selection(trialData, decOutcome, betOutcome, soa, decSaccDir, decTargDir, betSaccDir, betHighDir)
%
% Returns a list of the trial numbers with conditions specified in options
% structure. If called without any arguments, returns a default options structure.
% If options are input but one is not specified, it assumes default.
%
% Possible conditions are (default listed first):
%     options.decOutcome  = 'collapse','correct','incorrect','abort';
%     options.betOutcome  = 'collapse''high','low','abort';
%     options.soa    = 'collapse',<list of soa's: e.g. [17 34]>;
%     options.decSaccDir  = 'collapse','left','right','up','down','upLeft','upRight',downLeft','downRight'<list of possible angels, e.g. [45 -45];
%     options.decTargDir  = 'collapse','left','right','up','down','upLeft','upRight',downLeft','downRight'<list of possible angels, e.g. [45 -45];
%     options.betSaccDir  = 'collapse','left','right','up','down','upLeft','upRight',downLeft','downRight'<list of possible angels, e.g. [45 -45];
%     options.betHighDir  = 'collapse','left','right','up','down','upLeft','upRight',downLeft','downRight'<list of possible angels, e.g. [45 -45];





if nargin < 1
    selectOpt.decOutcome  = 'collapse';
    selectOpt.betOutcome  = 'collapse';
    selectOpt.soa         = 'collapse';
    selectOpt.decSaccDir  = 'collapse';
    selectOpt.decTargDir  = 'collapse';
    selectOpt.betSaccDir  = 'collapse';
    selectOpt.betHighDir  = 'collapse';
    trialList           = selectOpt;
    return
end




nTrial = size(trialData.trialOutcome, 1);
trialLogical = ones(nTrial, 1);





% Mask target locations: this will (I think) need to be distinguished
% throuhout the code- 2 targets vs. 4 targets
maskAngle  = trialData.decMaskAngle(1,:);
nDecTarg       = length(maskAngle);
leftTargInd = (maskAngle < -89) & (maskAngle > -270) | ...
    (maskAngle > 90) & (maskAngle < 269);
% switch nTarg
%     case 2
%         decAngleLeft = maskAngles(leftTargInd);
%         decAngleRight = maskAngles(~leftTargInd);
%     case 4
decAngleLeftUp = max(maskAngle(leftTargInd));
decAngleLeftDown = min(maskAngle(leftTargInd));
decAngleRightUp = max(maskAngle(~leftTargInd));
decAngleRightDown = min(maskAngle(~leftTargInd));
% end

% Bet target locations: this will (I think) need to be distinguished
% throuhout the code- 2 targets vs. 4 targets
betAngle  = unique(trialData.betHighAngle(~isnan(trialData.betHighAngle)));;
nBetTarg       = length(betAngle);
leftTargInd = (betAngle < -89) & (betAngle > -270) | ...
    (betAngle > 90) & (betAngle < 269);
% switch nTarg
%     case 2
%         decAngleLeft = maskAngles(leftTargInd);
%         decAngleRight = maskAngles(~leftTargInd);
%     case 4
decAngleLeftUp = max(maskAngle(leftTargInd));
decAngleLeftDown = min(maskAngle(leftTargInd));
decAngleRightUp = max(maskAngle(~leftTargInd));
decAngleRightDown = min(maskAngle(~leftTargInd));
% end




% ================================================
% Decision stage outcome
% ================================================
if isfield(selectOpt, 'decOutcome')
    switch selectOpt.decOutcome
        case 'collapse'
            trialLogical = trialLogical & ...
                (strcmp(trialData.decOutcome, 'target') | strcmp(trialData.decOutcome, 'distractor'));
        case {'target','distractor','abort'}
            trialLogical = trialLogical & strcmp(trialData.decOutcome, selectOpt.decOutcome);
        otherwise
            %do nothing
    end
end


% ================================================
% Bet stage outcome
% ================================================
if isfield(selectOpt, 'betOutcome')
    switch selectOpt.betOutcome
        case 'collapse'
            trialLogical = trialLogical & ...
                (strcmp(trialData.betOutcome, 'high') | strcmp(trialData.betOutcome, 'low'));
        case {'high','low','abort'}
            trialLogical = trialLogical & strcmp(trialData.betOutcome, selectOpt.betOutcome);
        otherwise
            %do nothing
    end
end






% ================================================
% Get list(s) of trials w.r.t. the SOA
% ================================================
if isfield(selectOpt, 'soa')
    if strcmp(selectOpt.soa, 'collapse')
        trialLogical = trialLogical & ~isnan(trialData.soa);
    else
        trialLogical = trialLogical & ismember(trialData.soa, selectOpt.soa);
    end
end








% ================================================
% Decision stage target Angle Range
% ================================================
if isfield(selectOpt, 'decTargDir')
    targAngle = trialData.decTargAngle;
    targTrial = ones(nTrial, 1);
    switch selectOpt.decTargDir
        case 'collapse'
            % Do nothing
        case 'right'
            targTrial = ((targAngle > 270) & (targAngle <= 360)) | ...
                ((targAngle >= 0) & (targAngle < 90)) | ...
                ((targAngle > -90) & (targAngle < 0)) | ...
                ((targAngle >= -360) & (targAngle < -270));
        case 'left'
            targTrial = ((targAngle > 90) & (targAngle <= 270)) | ...
                ((targAngle < -90) & (targAngle > -270));
        case 'leftUp'
            targTrial = targAngle == decAngleLeftUp;
        case 'leftDown'
            targTrial = targAngle == decAngleLeftDown;
        case 'rightUp'
            targTrial = targAngle == decAngleRightUp;
        case 'rightDown'
            targTrial = targAngle == decAngleRightDown;
        otherwise
            targTrial = ismember(targAngle,selectOpt.decTargDir);
    end
    trialLogical = trialLogical & targTrial;
end



% ================================================
% Bet stage High bet target Angle Range
% ================================================
if isfield(selectOpt, 'betHighDir')
    targAngle = trialData.betHighAngle;
    targTrial = ones(nTrial, 1);
    switch selectOpt.betHighDir
        case 'collapse'
            % Do nothing
        case 'right'
            targTrial = ((targAngle > 270) & (targAngle <= 360)) | ...
                ((targAngle >= 0) & (targAngle < 90)) | ...
                ((targAngle > -90) & (targAngle < 0)) | ...
                ((targAngle >= -360) & (targAngle < -270));
        case 'left'
            targTrial = ((targAngle > 90) & (targAngle <= 270)) | ...
                ((targAngle < -90) & (targAngle > -270));
        case 'leftUp'
            targTrial = targAngle == betHighAngleLeftUp;
        case 'leftDown'
            targTrial = targAngle == betHighAngleLeftDown;
        case 'rightUp'
            targTrial = targAngle == betHighAngleRightUp;
        case 'rightDown'
            targTrial = targAngle == betHighAngleRightDown;
        otherwise
            targTrial = ismember(targAngle,selectOpt.betHighDir);
    end
    trialLogical = trialLogical & targTrial;
end



% ================================================
% Decision stage saccade (response) Angle Range
% ================================================
if isfield(selectOpt, 'decSaccDir')
    saccAngle = trialData.decSaccAngle;
    decSaccTrial = ones(nTrial, 1);
    %     switch selectOpt.decSaccDir
    if strcmp(selectOpt.decSaccDir, 'collapse')
        % Do nothing
    elseif strcmp(selectOpt.decSaccDir, 'right')
        decSaccTrial = ((saccAngle > 270) & (saccAngle <= 360)) | ...
            ((saccAngle >= 0) & (saccAngle < 90)) | ...
            ((saccAngle > -90) & (saccAngle < 0)) | ...
            ((saccAngle >= -360) & (saccAngle < -270));
    elseif strcmp(selectOpt.decSaccDir, 'left')
        decSaccTrial = ((saccAngle > 90) & (saccAngle <= 270)) | ...
            ((saccAngle < -90) & (saccAngle > -270));
    elseif strcmp(selectOpt.decSaccDir, 'leftUp')
        decSaccTrial = saccAngle == decAngleLeftUp;
    elseif strcmp(selectOpt.decSaccDir, 'leftDown')
        decSaccTrial = saccAngle == decAngleLeftDown;
    elseif strcmp(selectOpt.decSaccDir, 'rightUp')
        decSaccTrial = saccAngle == decAngleRightUp;
    elseif strcmp(selectOpt.decSaccDir, 'rightDown')
        decSaccTrial = saccAngle == decAngleRightDown;
    else
        decSaccTrial = ismember(saccAngle,selectOpt.decSaccDir);
    end
    trialLogical = trialLogical & decSaccTrial;
end




% ================================================
% Bet stage saccade (response) Angle Range
% ================================================
if isfield(selectOpt, 'betSaccDir')
    saccAngle = trialData.betSaccAngle;
    betSaccTrial = ones(nTrial, 1);
    switch selectOpt.betSaccDir
        case 'collapse'
            % Do nothing
        case 'right'
            betSaccTrial = ((saccAngle > 270) & (saccAngle <= 360)) | ...
                ((saccAngle >= 0) & (saccAngle < 90)) | ...
                ((saccAngle > -90) & (saccAngle < 0)) | ...
                ((saccAngle >= -360) & (saccAngle < -270));
        case 'left'
            betSaccTrial = ((saccAngle > 90) & (saccAngle <= 270)) | ...
                ((saccAngle < -90) & (saccAngle > -270));
        case 'leftUp'
            betSaccTrial = saccAngle == betAngleLeftUp;
        case 'leftDown'
            betSaccTrial = saccAngle == betAngleLeftDown;
        case 'rightUp'
            betSaccTrial = saccAngle == betAngleRightUp;
        case 'rightDown'
            betSaccTrial = saccAngle == betAngleRightDown;
        otherwise
            betSaccTrial = ismember(saccAngle,selectOpt.betSaccDir);
    end
    trialLogical = trialLogical & betSaccTrial;
end




% If there was supposed to be an RT but there wasn't (if there is a NaN for
% RT, get rid of those trials






trialList = find(trialLogical);








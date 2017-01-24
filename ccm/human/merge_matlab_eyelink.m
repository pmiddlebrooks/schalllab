function merge_matlab_eyelink(sessionID, eyeOrKey)


% dataFolder = 'data/';
humanDataPath = '/Volumes/middlepg/HumanData/ChoiceStopTask/';
% humanDataPath = '~/schalllab/';

switch eyeOrKey
    case 'eye'
        effectorTag = 'em';
        effector = 'saccade';
    case 'key'
        effectorTag = 'kp';
        effector = 'keypress';
end
[humanDataPath, sessionID, effector, '.mat'];
MatlabData = dir([humanDataPath, sessionID, effector, '.mat']);
[humanDataPath, sessionID, effectorTag, '.asc'];
EyelinkData = dir([humanDataPath, sessionID, effectorTag, '.asc']);


% Translate the eyelink file associated with that session
sessionASC = [humanDataPath, EyelinkData.name];
pixelPerDegreeX = 31.5;
pixelPerDegreeY = 31.5;


[Experiment] = translate_eyelink_asc(sessionASC, pixelPerDegreeX, pixelPerDegreeY);
% Experiment(147,:) = [];
matlabFile = [humanDataPath, MatlabData.name];
load(matlabFile, 'trialData')

size(trialData.trialOutcome)
size(Experiment)

% If necessary, take off last trial from Experiment (means eyelink started
% a trial and so has data, but matlab cut the trial off because the task
% was ended mid-trial
trialData(strcmp(trialData.trialOutcome, 'eyelinkError'), :) = [];
if size(Experiment, 1) + 1 == size(trialData.trialOutcome, 1) && isnan(trialData.trialOutcome(end))
    %     Experiment(end,:) = [];
    trialData(end,:) = [];
end



% Some sessions, there are either eyelink errors (so there are extra
% early-aborted eyelink Experiment trials) or task errors (so there are extra
% early-aborted trialData trials). Find thos discrepant trials by look for
% cases when checker onset is mismatched between the two data sets (finding
% NaN's that aren't aligned).
t = find(isnan(trialData.checkerOn));
e = find(isnan(Experiment.eyelinkChoiceStimulusOnset));
if length(t) ~= length(intersect(t,e)) || length(e) ~= length(intersect(t,e))
    % Get rid of all mismatches
    mismatch = 1;
    while mismatch
        % Iteratively find each mismatch to do so
        findOne = 1;
    i = 1;
        while findOne
            if t(i) < e(i)
                trialData(t(i),:) = [];
                findOne = 0;
            elseif e(i) < t(i)
                Experiment(e(i),:) = [];
                findOne = 0;
            else
                i = i + 1;
            end
        end
        % If that fixes it (if there was only one discrepant trial, exit the
        % loop. Otherwise
        t = find(isnan(trialData.checkerOn));
        e = find(isnan(Experiment.eyelinkChoiceStimulusOnset));
        intersect(t,e)
        if length(t) == length(intersect(t,e)) && length(e) == length(intersect(t,e))
            mismatch = 0;
        end
    end
end


% Troubleshooting code:
% __________________________________________________________________
% [trialData.checkerOn, Experiment.eyelinkChoiceStimulusOnset(1:end-3)]
% dt = find(isnan(Experiment.eyelinkChoiceStimulusOnset));
% Experiment(dt(8),:) = []

% size(trialData.trialOutcome)
% size(Experiment)
% [t, e, intersect(t,e)]
% [trialData.checkerOn, Experiment.eyelinkChoiceStimulusOnset]
% __________________________________________________________________

trialData = [trialData, Experiment];

% Older code did not add stopSignalOn times when a response was made before
% the stop signal appeared. Fix that here.
if sum(~isnan(trialData.ssd) & isnan(trialData.stopSignalOn))
    r = find(~isnan(trialData.ssd) & isnan(trialData.stopSignalOn));
    trialData.stopSignalOn(r) = trialData.checkerOn(r) + trialData.ssd(r) + 7;
end
save(matlabFile, 'trialData', '-append');


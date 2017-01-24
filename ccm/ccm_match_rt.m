function data = ccm_match_rt(goRT, stopRT, nStopCorrect)

debug = 1;
testCode = 0;
if testCode
    %%
    Opt = ccm_options;
    Opt.plotFlag = false;
    
    Unit = ccm_session_data('Broca','bp093n02', Opt);
%     [trialData, SessionData] = load_data('Broca', 'bp093n02');
    %%
    sigInd          = 5;
    opt             = ccm_concat_neural_conditions;
    opt.epochName   = 'responseOnset';
    opt.eventMarkName = 'checkerOn';
    opt.colorCohArray = Unit(1).(opt.epochName).pSignalArray(sigInd);
    
    opt.conditionArray = {'goTarg'};
    dataGoTarg      = ccm_concat_neural_conditions(Unit(1), opt);
    
    opt.conditionArray = {'stopTarg'};
    opt.ssdArray = Unit(1).ssdArray;
    dataStopTarg    = ccm_concat_neural_conditions(Unit(1), opt);
    
    opt.conditionArray = {'stopStop'};
    opt.epochName   = 'stopSignalOn';
    opt.eventMarkName = 'checkerOn';
    dataStopStop    = ccm_concat_neural_conditions(Unit(1), opt);
    
    %%
%     nStopCorrect = size(dataStopStop.raster, 1);
    
    goRT    = dataGoTarg.eventLatency;
    stopRT  = dataStopTarg.eventLatency;
end

nGo       = length(goRT);
nStopTarg = length(stopRT);

% Get a random sample of the go trials (and RTs) equal in length to the stop trials
% goTrial   = randperm(nGo);
% goTrialSample   = goTrial(1 : nStopTarg + nStopCorrect);
% goRTSample      = goRT(goTrialSample);
% goTrialSample   = goTrial;
goTrialSample   = 1:length(goRT);
goRTSample      = goRT;

%%
% Create grid of all possible go/stop RT pairs, and a matrix of the go
% trial number available
[goMeshRT, stopMeshRT]      = meshgrid(goRTSample, stopRT);
[goMeshTrial, stopMeshTrial] = meshgrid(goTrialSample, 1:nStopTarg);
%%
% Take the difference between each possible pair to find the lowest
% difference (nearest-neighbor matched RTs)
[deltaRT, ind] = sort(abs(goMeshRT(:) - stopMeshRT(:)));
%%
% Sort trials in in same order as sorted deltaRT
goMeshTrial     = goMeshTrial(ind);
stopMeshTrial   = stopMeshTrial(ind);
goMeshRT        = goMeshRT(ind);
stopMeshRT      = stopMeshRT(ind);

%%
% Loop through the trials with lowest matched RT differences, and build a
% list of those trials without repeats
goFastTrial    = nan(nStopTarg, 1);
stopMatchTrial  = nan(nStopTarg, 1);
goFastRT       = nan(nStopTarg, 1);
stopMatchRT     = nan(nStopTarg, 1);
i = 1;
matchInd = 1;
while sum(isnan(goFastTrial))
    if ~ismember(goMeshTrial(i), goFastTrial) && ~ismember(stopMeshTrial(i), stopMatchTrial)
        goFastTrial(matchInd)      = goMeshTrial(i);
        stopMatchTrial(matchInd)    = stopMeshTrial(i);
        goFastRT(matchInd)         = goMeshRT(i);
        stopMatchRT(matchInd)       = stopMeshRT(i);
        matchInd = matchInd + 1;
    end
    i = i + 1;
end
goSlowTrial = setdiff(unique(goMeshTrial), goFastTrial);
goSlowRT = goRT(goSlowTrial);

data.goFastTrial = goFastTrial;
data.goSlowTrial = goSlowTrial;
data.stopMatchTrial = stopMatchTrial;




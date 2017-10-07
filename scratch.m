%%
figure(10);
hold all
cellfun(@(x,y) plot(x,y, 'k'), goEyeX{1,1}, goEyeY{1,1})

% plot(goEyeX{1,1}, goEyeY{1,1})

degX = cellfun(@(x) sqrt(x), cellfun(@(x,y) x^2 + y^2, goEyeX{1,1}, goEyeY{1,1}))


velX{d, i} = cellfun(@(x) [0; diff(x(:))], degX, 'uni', false);


scatter(trialData.rt(goTrial{d, i}), cellfun(@max, goVel{d, i}))


%% Add hemisphere to translated data file
subject = 'joule';
session = {'jp061n02'};
tebaPath = '/Volumes/SchallLab/data/';


for i = 1 : length(session)
    [~, SessionData] = load_data(subject, session{i});
    
    SessionData.hemisphere = 'left';
    
    save(fullfile(local_data_path, subject, [session{i}, '.mat']), 'SessionData', '-append')
    save(fullfile(tebaPath, subject, [session{i}, '.mat']), 'SessionData', '-append')
    
    
end


%%
sdf53 = nanmean(spike_density_function(cell2mat(alignedRasters(trialData.targ1CheckerProp == .53,:)), Kernel), 1);
sdf58 = nanmean(spike_density_function(cell2mat(alignedRasters(trialData.targ1CheckerProp == .58,:)), Kernel), 1);
figure(1)
eopchLBegin = epochBegin(trialData.targ1CheckerProp == .58);
epochLEnd = epochEnd(trialData.targ1CheckerProp == .58);
for i = 1 : length(epochBegin)
    
    clf
    hold on
    plot(sdf53);
    plot(sdf58, 'k')
    plot([eopchLBegin(i) eopchLBegin(i)], [0 50])
    plot([epochLEnd(i) epochLEnd(i)], [0 50])
    pause
end
%%
opt = ccm_options;
opt.trialOutcome = 'valid';
opt.ssd = 'collapse';
trial = ccm_trial_selection(td, opt);
td = td(trial,:);
nSsd = arrayfun( @(x)(length(find(td.ssd==x))), E.ssdArray);
[E.ssdArray, nSsd]
criteriaDiff = 20;

belowCritera = find(diff(E.ssdArray) < criteriaDiff);

% Are there any that have runs of more than 2 SSDs that are less than
% criteria? If so, we need to give up and keep one.
remove = 1+find(diff(belowCritera) < 2);
belowCritera(remove) = [];

ssdIndAltered = [belowCritera; belowCritera+1];
ssdKeep = setxor(ssdIndAltered, 1:length(E.ssdArray));

ssdWeighted = nan(length(belowCritera), 1);
for i = 1 : length(belowCritera)
    
    ssdInd = [belowCritera(i) belowCritera(i)+1];
    ssdWeighted(i) = round(sum(E.ssdArray(ssdInd) .* nSsd(ssdInd) / sum(nSsd(ssdInd))));
    td.ssd(td.ssd == E.ssdArray(ssdInd(1)) | td.ssd == E.ssdArray(ssdInd(2))) = ssdWeighted(i);
end

newSSD = sort([E.ssdArray(ssdKeep); ssdWeighted]);
newSSD2 = unique(td.ssd);

[i,ia] = ismember(newSSD, td.ssd)
nSsdNew = arrayfun( @(x)(length(find(td.ssd==x))), newSSD);
[newSSD, nSsdNew]

%%
[td1, S, E] = load_data('broca','bp256n01');
[td2, S, E] = load_data('broca','bp256n02');
[td3, S, E] = load_data('broca','bp256n03');
% [td1, S, E] = load_data('broca','bp255n01');
% [td2, S, E] = load_data('broca','bp255n02');
% [td3, S, E] = load_data('broca','bp255n03');

opt = ccm_options;

opt.trialOutcome = 'goCorrectTarget';
% opt.targAngle = 0;
opt.targHemifield = 'right';

rCorr1 = cmd_trial_selection(td1, opt);
rCorr2 = cmd_trial_selection(td2, opt);
rCorr3 = cmd_trial_selection(td3, opt);

rtRCorr1 = nanmean(td1.rt(rCorr1))
rtRCorr2 = nanmean(td2.rt(rCorr2))
rtRCorr3 = nanmean(td3.rt(rCorr3))


% opt.targAngle = 180;
opt.targHemifield = 'left';

lCorr1 = cmd_trial_selection(td1, opt);
lCorr2 = cmd_trial_selection(td2, opt);
lCorr3 = cmd_trial_selection(td3, opt);

rtLCorr1 = nanmean(td1.rt(lCorr1))
rtLCorr2 = nanmean(td2.rt(lCorr2))
rtLCorr3 = nanmean(td3.rt(lCorr3))

leftRT = [td1.rt(lCorr1); td2.rt(lCorr2); td3.rt(lCorr3)];
leftGroup = [ones(length(td1.rt(lCorr1)), 1); 2*ones(length(td2.rt(lCorr2)), 1); 3* ones(length(td3.rt(lCorr3)), 1)];
rightRT = [td1.rt(rCorr1); td2.rt(rCorr2); td3.rt(rCorr3)];
rightGroup = [ones(length(td1.rt(rCorr1)), 1); 2*ones(length(td2.rt(rCorr2)), 1); 3* ones(length(td3.rt(rCorr3)), 1)];
% [pL, tableL, statsL] = anova1(leftRT, leftGroup, 'display', 'off');
% [pR, tableR, statsR] = anova1(rightRT, rightGroup, 'display', 'off');
save(fullfile(local_data_path, 'anodal.mat'), 'leftRT', 'leftGroup', 'rightRT', 'rightGroup')

%%
figure(1)
hold all
plot([rtLCorr1, rtLCorr2, rtLCorr3], '--k')
plot([rtRCorr1, rtRCorr2, rtRCorr3], '--b')


%%
ssd = new.ssd;

ssdList = unique(ssd(~isnan(ssd)))
nSSD = nan(length(ssdList), 1);
for i = 1 : length(ssdList)
    nSSD(i) = sum(trialData.ssd == ssdList(i));
end

%%

% plot Average unit dynamics
dynTime = mean(cell2mat(prd.dyn{iTrialCatGo}.stopICorr.stopStim.targetGO.sX));
dynAct = mean(cell2mat(prd.dyn{iTrialCatGo}.stopICorr.stopStim.targetGO.sY));
plot(dynTime, dynAct, 'Color','k','LineStyle',unitLnStyle{unitGoCorr},'LineWidth',unitLnWidth(unitGoCorr))

dynTime = prd.dyn{iTrialCatGo}.stopICorr.goStim.targetGO.sX;
dynAct = prd.dyn{iTrialCatGo}.stopICorr.goStim.targetGO.sY;
cellfun(@(x,y) plot(x,y, 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',1), dynTime,dynAct, 'uni', false)


plot(dynTime(iSsd+1:end), dynAct(iSsd+1:end), 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',unitLnWidth(unitStop))
plot(dynTime(iSsd+1:end), dynAct(iSsd+1:end), 'Color','b','LineStyle',unitLnStyle{unitStop},'LineWidth',unitLnWidth(unitStop))

normAct = cell(size(dynActGoCorr));
lastNonNan = cellfun(@(x) find(isnan(x), 1), dynActGoCorr);
lastNonNan = max(lastNonNan) - 1;
for i = 1 : size(dynActGoCorr, 1)
    iLast = find(isnan(dynActGoCorr{i}), 1) - 1;
    iScale = lastNonNan / iLast;
    normAct{i} = dynActGoCorr{i} * iScale;
end
cellfun(@(x,y) plot(x,y, 'Color','k','LineStyle',unitLnStyle{unitGoCorr},'LineWidth',1), dynTimeGoCorr,normAct, 'uni', false)

meanFn = nanmean(cell2mat(normAct));

%%

% [trialData, S, E] = load_data('joule','jp098n02');
[trialData, S, E] = load_data('broca','bp244n02');
%%

pSignalArray    = E.pSignalArray;
mEpochName = 'checkerOn';
mEpochName = 'responseOnset';
Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;

[unitIndex, unitArrayNew] = neuronexus_plexon_mapping(S.spikeUnitArray, 32);


selectOpt = ccm_options;
selectOpt.rightCheckerPct = pSignalArray(end) * 100;
selectOpt.ssd = 'none';
selectOpt.outcome     = {'goCorrectTarget', 'targetHoldAbort', 'toneAbort'};


iGoTrial = ccm_trial_selection(trialData, selectOpt);
alignListGo = trialData.(mEpochName)(iGoTrial);

sdf = [];
epochWindow = [-199:400];
normWindow = [-299:0];
for kUnit = unitIndex
    alignListTarg = trialData.targOn(iGoTrial);
    [alignedRasters, ~] = spike_to_raster(trialData.spikeData(iGoTrial, kUnit), alignListTarg);
    iSdfTarg = nanmean(spike_density_function(alignedRasters, Kernel));
    iNormTarg = max(iSdfTarg);
    
    alignListSacc = trialData.responseOnset(iGoTrial);
    [alignedRasters, ~] = spike_to_raster(trialData.spikeData(iGoTrial, kUnit), alignListSacc);
    iSdfSacc = nanmean(spike_density_function(alignedRasters, Kernel));
    iNormSacc = max(iSdfSacc);
    
    alignListChecker = trialData.targOn(iGoTrial);
    [alignedRasters, ~] = spike_to_raster(trialData.spikeData(iGoTrial, kUnit), alignListChecker);
    iSdfChecker = nanmean(spike_density_function(alignedRasters, Kernel));
    iNormChecker = max(iSdfChecker);
    
    iNorm = max([iNormTarg, iNormSacc, iNormChecker]);
    
    [alignedRasters, alignmentIndex] = spike_to_raster(trialData.spikeData(iGoTrial, kUnit), alignListGo);
    iSDF = nanmean(spike_density_function(alignedRasters, Kernel)) / iNorm;
    sdf = [sdf; iSDF(alignmentIndex + epochWindow)];
end

%%
for i = 1 : length(Data)
    disp(Data(i).name)
    disp(Data(i).yMax)
    pause
end
%%
td = trialData;
td = td(isnan(td.abortTime),:);

stopTrial = strcmp(td.trialOutcome, 'stopIncorrectTarget') | strcmp(td.trialOutcome, 'stopIncorrectDistractor');
tdStop = td(stopTrial,:);
goTrial = strcmp(td.trialOutcome, 'goCorrectTarget') | strcmp(td.trialOutcome, 'goCorrectDistractor');
tdGo = td(goTrial,:);
%%
fprintf('\n\n%d', iPropIndexL)
fprintf('\nGo: %d',size(Data(kDataIndex, jTarg).(mEpochName).colorCoh(iPropIndexL).goTarg.sdf, 1) + size(Data(kDataIndex, jTarg).(mEpochName).colorCoh(iPropIndexL).goDist.sdf, 1))
fprintf('\nStop: %d', size(rasStopTarg, 1) + size(rasStopDist, 1))
fprintf('\nCancel: %d\n',size(rasStopCorrect, 1))
%%
fprintf('\n\n%d', iPropIndexR)
fprintf('\nGo: %d',size(Data(kDataIndex, jTarg).(mEpochName).colorCoh(iPropIndexR).goTarg.sdf, 1) + size(Data(kDataIndex, jTarg).(mEpochName).colorCoh(iPropIndexR).goDist.sdf, 1))
fprintf('\nStop: %d', size(rasStopTarg, 1) + size(rasStopDist, 1))
fprintf('\nCancel: %d\n',size(rasStopCorrect, 1))

%%
load(fullfile(local_data_path,'broca/broca_behavior2_ssd'))
td = trialData;

opt = ccm_options;
opt.outcome = 'valid';
validTrial = ccm_trial_selection(td, opt);
td = td(validTrial,:);
td(td.responseOnset - td.responseCueOn > 1200, :) = [];
nTrial = size(td, 1);

subj_idx = td.sessionTag;

ss_presented = zeros(nTrial, 1);
ss_presented(~isnan(td.stopSignalOn)) = 1;

inhibited = zeros(nTrial, 1);
inhibited(ss_presented & strcmp(td.trialOutcome, 'stopCorrect')) = 1;
inhibited(~ss_presented) = -999;

ssd = td.ssd;
ssd(~ss_presented) = -999;

rt = td.responseOnset - td.responseCueOn;
rt(ss_presented & strcmp(td.trialOutcome, 'stopCorrect')) = -999;

b = table(subj_idx, ss_presented, inhibited, ssd, rt);

writetable(b, fullfile(local_data_path, 'broca/bayes_ssrt.csv'))

%% Load trialData and SessionData from one of joule's sessions.
% ============================================================
% trialData is a table: trial by trial (each row is a trial) table of various task- related values.
% SessionData is a struct: contains some information about the session that was recorded.
tebaDataPath = '/Volumes/SchallLab/data/';  % You may need to change this line
load(fullfile(tebaDataPath,'joule/jp111n02.mat'));

%% Use the data to create some variables for Jacob to do stuff with (include spike_to_raster function)
% ============================================================

% Gather some variables directly from trialData:
rt              = trialData.rt;  % response times (a vector)
colorCoherence  = trialData.targ1CheckerProp; % the proportion of right target color in the checkerboard stimulus (a vector)
trialOutcome    = trialData.trialOutcome; % trial outcomes (a cell of strings)
trialType       = trialData.trialType; % stop or go trial (a cell of strings)

% Calculate spike rates for a given period of time durin the trial (epoch) aligned on some event (alignList)
nUnit           = length(SessionData.spikeUnitArray);
alignList       = trialData.responseOnset;
epoch           = 101:300;

% Initialize a matrix to be filled in with spike rates for each unit.
spikeRate       = nan(size(trialData.trialOutcome, 1), nUnit);

% Loop through each unit and get the spike rates
for kUnit = 1 : nUnit
    [alignedRasters, alignmentIndex]    = spike_to_raster(trialData.spikeData(:, kUnit), alignList);
    epochRasters                        = alignedRasters(:, alignmentIndex + epoch);
    kSpikeRate                          = sum(epochRasters,2) ./ (length(epoch) / 1000);  % in spikes per second (1000 = ms to sec conversion)
    spikeRate(:, kUnit)                 = kSpikeRate;
end


%%

append = false;

subject = 'joule';
ccm_classify_neuron_pop(subject,projectRoot,projectDate,append)

subject = 'broca';
ccm_classify_neuron_pop(subject,projectRoot,projectDate,append)

%%
opt = ccm_options;
opt.printPlot = false;
subject = 'broca';

sessionID = {'bp095n04'};
unit = {'spikeUnit17b'};
iUnit = [sessionID, unit];
%%
iData               = ccm_session_data(subject, iUnit);

iCat              	= ccm_classify_neuron(iData);
%%
for iVar = 1 : length(variables)
    trialData.(variables{iVar}) = [];
end

%%
tic
[td, S, E] = load_data('joule','jp113n02');
toc

%%
tic
td = load('~/schalllab/scratch/jp113n02');
toc
%%
tic
td = load('~/schalllab/scratch/jp113n02', 'targOn');

toc
%%
tic
variables = {...
    'trialOutcome',...
    'targ1CheckerProp',...
    'responseDirection',...
    'saccToTargIndex',...
    'saccAngle',...
    'rt',...
    'ssd',...
    'targAngle',...
    };
% variables = [variables, {'spikeData'}];
td = load('~/schalllab/scratch/jp113n02', variables{:});
toc


%%
%% Verifying/discovering alignment discrepancy between plexon spike times and reality spike times based on observed shift between electrode channels
sigE = [];
sigL = [];
for i = 1 : 21
    iSigE = Data(i).responseOnset.colorCoh(2).goTarg.sdfMean(Data(i).responseOnset.colorCoh(2).goTarg.alignTime + [-200 : 200]);
    sigE = [sigE; iSigE];
end

for i = 22 : 89
    iSigL = Data(i).responseOnset.colorCoh(2).goTarg.sdfMean(Data(i).responseOnset.colorCoh(2).goTarg.alignTime + [-200 : 200]);
    sigL = [sigL; iSigL];
end

sigEM = nanmean(sigE, 1);
sigLM = nanmean(sigL, 1);

[maxE, iE] = max(sigEM);
[maxL, iL] = max(sigLM);

sigEM = sigEM/maxE;
sigLM = sigLM/maxL;
[acor, lag] = xcorr(sigEM', sigLM')

%%

figure(1), axis equal, axis off

A = [size(classic, 1) size(cancel, 1) size(ddm, 1)];
I = [size(classicCancel, 1) size(classicDdm, 1) size(classicDdmCancel, 1) size(classicDdmCancel, 1)];
venn(A,I,'FaceColor',{'r','y','b'},'FaceAlpha',{.5,0.6,0.2},'EdgeColor','black')

%%
%Compare ErrMinModes
A = [350 300 275]; I = [100 80 60 40];
figure
subplot(1,3,1), h1 = venn(A,I,'ErrMinMode','None');
axis image,  title ('No 3-Circle Error Minimization')
subplot(1,3,2), h2 = venn(A,I,'ErrMinMode','TotalError');
axis image,  title ('Total Error Mode')
subplot(1,3,3), h3 = venn(A,I,'ErrMinMode','ChowRodgers');
axis image, title ('Chow-Rodgers Mode')
set([h1 h2], 'FaceAlpha', 0.6)

%Using the same areas as above, display the error optimization at each
%   iteration. Get the output structure.
F = struct('Display', 'iter');
[H,S] = venn(A,I,F,'ErrMinMode','ChowRodgers','FaceAlpha', 0.6);

%Now label each zone
for i = 1:7
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), ['Zone ' num2str(i)])
end




%%
subject = 'broca';
sessionSet = 'behavior1';
sessionList = task_session_array(subject, 'ccm',sessionSet);

% subject = 'joule';
% sessionList = unique(neuronTypes.sessionID);
% deleteSession = {'jp054n02', 'jp060n02', 'jp061n02'};
% sessionList(ismember(sessionList, deleteSession)) = [];


for i = 1 : length(sessionList)
    fprintf('\n\n\n\n%s\n\n\n\n', sessionList{i});
    %     [td s] = load_data(subject, sessionList{i}, ccm_min_vars);
    data = ccm_session_behavior(subject, sessionList{i});
    data = ccm_inhibition_rt(subject, sessionList{i}, opt);
end

%%

[td s] = load_data('joule', 'jp113n02',[ccm_min_vars, 'trialOnset', 'trialDuration','rewardDuration']);
%%
sessionDuration = floor((td.trialOnset(end) + td.trialDuration(end)) / 1000/60);
rewardRate = nan(sessionDuration, 1);

rewardDuration = cellfun(@sum, td.rewardDuration);

maMinute = 5;
maWindow = maMinute * 60 * 1000;

for i = maMinute : sessionDuration
    iMsBegin = (i - maMinute) * 60000;
    earliestTrial = find(td.trialOnset >= iMsBegin, 1);
    iMsEnd = i * 60000;
    latestTrial = find(td.trialOnset + td.trialDuration <= iMsEnd, 1, 'last');
    
    rewardRate(i) = sum(rewardDuration(earliestTrial : latestTrial)) / maMinute;
    
end


%%
subject = 'xena';
sessionList = task_session_array(subject, 'ccm','behavior1');

tOpt = plexon_translate_datafile_mac;
tOpt.whichData = 'behavior';

for i = 1 : length(sessionList)
    plexon_translate_datafile_mac('xena', sessionList{i},tOpt);
end

%%
target = 7000;
pricePerCustomer = 250;
nCustomerPerOptin = .08;
nOptinPerDay = 5;
nDays = 21;

revenue = pricePerCustomer * nCustomerPerOptin * nOptinPerDay * nDays;
% cost = spendPerDay * nDays;

spendPerDay = (target - revenue) / nDays
%%
spendPerDay = 14;
cost = spendPerDay * nDays;


%%
subject = 'broca';
sessionSet = {'bp121n04', 'bp100n01', 'bp092n02', 'bp228n02'};

opt = ccm_options;
opt.printFlag = 1;

for i = 1 : length(sessionSet)
    fprintf('\n\n\n\n%s\n\n\n\n', sessionSet{i});
    %     [td s] = load_data(subject, sessionList{i}, ccm_min_vars);
    data = ccm_session_behavior(subject, sessionSet{i});
    data = ccm_inhibition_rt(subject, sessionSet{i}, opt);
end

subject = 'joule';
sessionSet = {'jp114n04', 'jp106n02', 'jp104n02', 'jp083n02', 'jp110n02', 'jp125n04'};

opt = ccm_options;
opt.printFlag = 1;

for i = 1 : length(sessionSet)
    fprintf('\n\n\n\n%s\n\n\n\n', sessionSet{i});
    %     [td s] = load_data(subject, sessionList{i}, ccm_min_vars);
    data = ccm_session_behavior(subject, sessionSet{i});
    data = ccm_inhibition_rt(subject, sessionSet{i}, opt);
end


%%
%%matlab

subject = 'joule';

projectDate = '2017-01-11';
projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';

addpath(genpath(fullfile(projectRoot,'src/code',projectDate)));
dataPath = fullfile(projectRoot,'data',projectDate,subject);

categoryList = {'presacc','presaccNoVis','presaccRamp','visPresacc'};
categoryList = {'presacc'};

for i = 1 : length(categoryList)
    
    load(fullfile(dataPath, ['ccm_',categoryList{i},'_neurons']))
    
    % load the population of cancel time anlysis
    load(fullfile(dataPath, ['ccm_canceled_vs_go_neuronTypes']))
    
    % Build a new table of the relevant neurons, and a list of the session/Unit
    cancelData = table();
    for i = 1 : size(neurons, 1)
        % find the indices in cancelTypes that correspond to this unit
        iInd = strcmp(neurons.sessionID(i), cancelTypes.sessionID) & strcmp(neurons.unit(i), cancelTypes.unit);
        cancelData = [cancelData; cancelTypes(iInd,:)];
        
        
    end
end


%%
subject = 'joule';
sessionList = {'jp116n01', 'jp119n01', 'jp121n01', 'jp123n01', 'jp124n01', 'jp125n01'};
sessionList = {'jp106n01'};

tOpt = plexon_translate_datafile_mac;
tOpt.hemisphere = 'left';

for i = 1 : length(sessionList)
    plexon_translate_datafile_mac(subject, sessionList{i},tOpt);
end


%%
subject = 'broca';
sessionSet = {'bp178n02'};
sessionSet = {'bp228n02'};

opt = ccm_options;
opt.printFlag = 1;

for i = 1 : length(sessionSet)
    fprintf('\n\n\n\n%s\n\n\n\n', sessionSet{i});
    %     [td s] = load_data(subject, sessionList{i}, ccm_min_vars);
    data = ccm_session_behavior(subject, sessionSet{i});
    %     data = ccm_inhibition_rt(subject, sessionSet{i}, opt);
end
%%  CREATING PLOTS FOR SSRT ANALYSES

subject = 'broca';
sessionSet = 'behavior2';
% sessionList = task_session_array(subject, 'ccm',sessionSet);
sessionList = 'bp178n02';
sessionList = 'bp183n02';
% sessionList = 'bp228n02';

opt = ccm_options;
opt.plotFlag = 1;
opt.printPlot = 1;

% for i = 1 : length(sessionList)
%     fprintf('\n\n\n\n%s\n\n\n\n', sessionList{i});
% data = ccm_inhibition(subject, sessionList, opt);

% weibullParams = cell2mat(data.weibullParams);

data = ccm_inhibition_rt(subject, sessionList, opt);


%% matlab
subject = 'broca';

projectDate = '2017-01-11';
accreRoot = '/gpfs22';
accreHome = '/home/middlepg';
accreScratch = '/scratch/middlepg';
if isdir(fullfile(accreScratch))
    matRoot = fullfile(accreRoot,accreHome,'m-files'); % Edit this if running on Accre
    projectRoot = fullfile(accreScratch,'perceptualchoice_stop_model');
    environ = 'accre';
else
    matRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/schallab';
    projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';
    environ = 'local';
end

addpath(genpath(fullfile(matRoot,'ccm')));
addpath(genpath(fullfile(projectRoot,'src/code',projectDate)));


projectDate = '2017-01-11';
projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';

addpath(genpath(fullfile(projectRoot,'src/code',projectDate)));
dataPath = fullfile(projectRoot,'data',projectDate,subject);


load(fullfile(dataPath, 'ccm_neuronTypes'))

sessionList = unique(neuronTypes.sessionID);
% deleteSession = {'jp054n02', 'jp060n02', 'jp061n02'};
% deleteSession = {'bp080n01', 'xxxx'};
% sessionList(ismember(sessionList, deleteSession)) = [];


%%

for i = 1 : 1%length(sessionList)
    fprintf('\n\n\n\n%s\n\n\n\n', sessionList{i});
    %         data = ccm_session_behavior(subject, sessionList{i});
    data = ccm_inhibition_ssd_metrics(subject, sessionList{i});
    %     data = ccm_inhibition_rt(subject, sessionList{i});
    
    
    
    
end

%%    SSRT across color coherence within each SSD - Population
opt = ccm_options;
opt.plotFlag = 0;
opt.printPlot = 0;



for i = 1 : 1%length(sessionList)
    fprintf('\n\n\n\n%s\n\n\n\n', sessionList{i});
    data = ccm_session_behavior(subject, sessionList{i});
    data = ccm_inhibition_rt(subject, sessionList{i});
    
    
    
    
    
end





%%
figure(5);
hold on
plot(1:length(data.weibullParams), weibullParams(:,1), '.k')
plot(1:length(data.weibullParams), weibullParams(:,2), '.b')
%             [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);


%%
optInh = ccm_options;
optInh.printPlot = 0;
optInh.plotFlag = 0;
inh = ccm_inhibition('broca', 'bp228n02', optInh);

chron = ccm_chronometric('broca','bp228n02', optInh);
%%
inhTargRT = cellfun(@(x) nanmean(x), inh.goTargRT(:,1), 'uni', false);
cellfun(@(x) nanmean(x), chron.goLeftToTarg, 'uni', false)
cellfun(@(x) nanmean(x), chron.goRightToTarg, 'uni', false)


%%
ssdArray = unique(cell2mat(cancelData.stopStopSsd));


for i = 1 :length(ssdArray)
    hardInd = cellfun(@(x) x == .58, cancelData.stopStopCoh, 'uni', false);
    ssdInd =  cellfun(@(x) x == ssdArray(i), cancelData.stopStopSsd, 'uni', false);
    
    iHardCoh = cellfun(@(x,y,z) x(y & z), cancelData.stopStopCoh, hardInd, ssdInd, 'uni', false);
    iSSD = cellfun(@(x,y,z) x(y & z), cancelData.stopStopSsd, hardInd, ssdInd, 'uni', false);
    iNeuralCancelTime = cellfun(@(x,y,z,k) x(y & z) - k(y & z), cancelData.cancelTime2Std, hardInd, ssdInd, cancelData.stopStopSsd, 'uni', false);
    % iSsrt = cellfun(@(x) x - cancel
end

%% matlab
subject = 'joule';

projectDate = '2017-01-11';
accreRoot = '/gpfs22';
accreHome = '/home/middlepg';
accreScratch = '/scratch/middlepg';
if isdir(fullfile(accreScratch))
    matRoot = fullfile(accreRoot,accreHome,'m-files'); % Edit this if running on Accre
    projectRoot = fullfile(accreScratch,'perceptualchoice_stop_model');
    environ = 'accre';
else
    matRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/schallab';
    projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';
    environ = 'local';
end

addpath(genpath(fullfile(matRoot,'ccm')));
addpath(genpath(fullfile(projectRoot,'src/code',projectDate)));


projectDate = '2017-01-11';
projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';

addpath(genpath(fullfile(projectRoot,'src/code',projectDate)));
dataPath = fullfile(projectRoot,'data',projectDate,subject);


load(fullfile(dataPath, 'ccm_neuronTypes'))

sessionList = unique(neuronTypes.sessionID);
% deleteSession = {'jp054n02', 'jp060n02', 'jp061n02'};
% deleteSession = {'bp080n01', 'xxxx'};
% sessionList(ismember(sessionList, deleteSession)) = [];


%%
subject = 'xena';
sessionSet = 'behavior1';
sessionList = task_session_array(subject, 'ccm',sessionSet)
%%
subject = 'broca';
sessionSet = 'behavior2';
sessionList = task_session_array(subject, 'ccm',sessionSet)
%%
data = ccm_inhibition_population(subject, sessionList)


%%
[td s] = load_data('joule', 'jp125n04', [ccm_min_vars,'spikeData']);

%%              Memory-guided spiking data for Kaleb
% =======================================================================
% Which computer are you on?

if isdir('/Volumes/HD-1/Users/paulmiddlebrooks/')
    projectRoot = '~/memory_guided_saccades';
elseif isdir('/Volumes/Macintosh HD/Users/elseyjg/')
    projectRoot = '/Volumes/Macintosh HD/Users/elseyjg/Memory-Guided-Saccade-Project';
else
    disp('You need to add another condition or the file path is wrong.')
end

dataRoot = fullfile(projectRoot, 'data');

%%
cd('~/schalllab')
subject = 'joule';

load(fullfile(dataRoot, subject, 'mem_units.mat'))
sessionID = units.sessionID;
unit = units.unit;

opt = mem_options;
opt.printPlot = true;


sdfVis = cell(length(sessionID), 1);
sdfMov = cell(length(sessionID), 1);
visEpochWindow = -300:500;
movEpochWindow = -500:300;
poolID = parpool(4);
parfor i = 1 : length(sessionID)
    fprintf('%s\t%s\n',sessionID{i}, unit{i})
    iUnit = [sessionID(i), unit(i)];
    iData = mem_session_data(subject, iUnit, opt);
    if isnan(iData.rf)
        if strcmp(units.hemisphere{i}, 'left')
            iData.rf = 'right';
        else
            iData.rf = 'left';
        end
    end
    sdfVis{i} = iData.([iData.rf,'Targ']).targOn.signalMean(iData.([Data.rf,'Targ']).targOn.alignTime + visEpochWindow);
    sdfMov{i} = iData.([iData.rf,'Targ']).responseOnset.signalMean(iData.([Data.rf,'Targ']).responseOnset.alignTime + movEpochWindow);
end
delete(poolID)
%%

session = 'jp110n01';
session = 'jp121n01';
session = 'jp124n01';

unitInd = strcmp(sessionID, session);
unitList = unit(unitInd);

opt = mem_options;
opt.printPlot = true;
opt.multiUnit = true;

for i = 1 : 32
    
    iUnitName = sprintf('spikeUnit%.2d', i);
    iUnit = {session, iUnitName};
        iData = mem_session_data(subject, iUnit, opt);
end
 
%%
iUnit = {'jp110n01', 'spikeUnit01'};
iData = mem_session_data(subject, iUnit, opt)
%%
sdf = [sdfVis, sdfMov];
time = [repmat({visEpochWindow}, length(sessionID), 1), repmat({movEpochWindow}, length(sessionID), 1)];
save([local_data_path, 'mem_data_joule'], 'sdf', 'time', 'sessionID', 'unit')




%%              Delaye saccade spiking data for Kaleb

load(fullfile(dataRoot, subject, 'mem_units.mat'))
sessionList = unique(units.sessionID);
subject = 'joule';
% sessionList = {'jp116n01', 'jp119n01', 'jp121n01', 'jp123n01', 'jp124n01', 'jp125n01'};
% sessionList = {'jp106n01'};

tOpt = plexon_translate_datafile_mac;
tOpt.hemisphere = 'left';

for i = 1 : length(sessionList)
    plexon_translate_datafile_mac(subject, sessionList{i},tOpt);
end

%%
% =======================================================================
% Which computer are you on?

if isdir('/Volumes/HD-1/Users/paulmiddlebrooks/')
    projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/memory_guided_saccades';
elseif isdir('/Volumes/Macintosh HD/Users/elseyjg/')
    projectRoot = '/Volumes/Macintosh HD/Users/elseyjg/Memory-Guided-Saccade-Project';
else
    disp('You need to add another condition or the file path is wrong.')
end

dataRoot = fullfile(projectRoot, 'data');

%%

subject = 'broca';
%%
load(fullfile(dataRoot, subject, 'del_units.mat'))
sessionID = units.sessionID;
unit = units.unit;

opt = mem_options;
opt.printPlot = true;
opt.task = 'del';


sdfVis = cell(length(sessionID), 1);
sdfMov = cell(length(sessionID), 1);
visEpochWindow = -300:500;
movEpochWindow = -500:300;
poolID = parpool(4);
parfor i = 1 : length(sessionID)
    fprintf('%s\t%s\n',sessionID{i}, unit{i})
    iUnit = [sessionID(i), unit(i)];
    iData = mem_session_data(subject, iUnit, opt);
    if isnan(iData.rf)
        if strcmp(units.hemisphere{i}, 'left')
            iData.rf = 'right';
        else
            iData.rf = 'left';
        end
    end
    sdfVis{i} = iData.([iData.rf,'Targ']).targOn.signalMean(iData.([Data.rf,'Targ']).targOn.alignTime + visEpochWindow);
    sdfMov{i} = iData.([iData.rf,'Targ']).responseOnset.signalMean(iData.([Data.rf,'Targ']).responseOnset.alignTime + movEpochWindow);
end
delete(poolID)

%%
sdf = [sdfVis, sdfMov];
time = [repmat({visEpochWindow}, length(sessionID), 1), repmat({movEpochWindow}, length(sessionID), 1)];
save([local_data_path, 'del_data_broca'], 'sdf', 'time', 'sessionID', 'unit')


%%
lowRate = [1 10 50];
highRate = [1 10 50 100];
nTrial = [5 50 500];

pSignRank = nan(length(lowRate), length(highRate), length(nTrial));
pTTest = nan(length(lowRate), length(highRate), length(nTrial));
for i = 1 : length(lowRate)
    iHighRateInd = highRate >= lowRate(i);
%     iHighRate = highRate(highRate >= lowRate(i));
    for j = find(iHighRateInd) : length(highRate)
        for k = 1 : length(nTrial);
        jLowRand = poissrnd(lowRate(i), nTrial(k), 1);
        jHighRand = poissrnd(highRate(j), nTrial(k), 1);
        
        [p,h,stats]   = signrank(jLowRand, jHighRand);
        pSignRank(i,j,k) = p;

        [h,p,ci,stats] = ttest(jLowRand, jHighRand);
        pTTest(i,j,k) = p; 
        end
    end
end

        
  
%%
[alignedRasters, alignmentIndex] = spike_to_raster(trialData.spikeData(trial, spikeDataInd), trialData.responseOnset(trial));
Kernel.method = 'postsynaptic potential';
Kernel.growth = 1;
Kernel.decay = 20;
sdf = spike_density_function(alignedRasters, Kernel);
sdf = sdf(:,alignmentIndex+(-200:0));

spikeData = cellfun(@(in1,in2,in3) sort(in1(in1 > in2 & in1 < in3)), trialData.spikeData(trial,spikeDataInd), saccadeWinBegin, saccadeWinEnd, 'uni', false);

%%
clf

for i = 1 : size(sdf, 1)
 figure(9)
   plot(sdf(i,:))
    Data.saccadeSpikeInterval{i}
 figure(10)
   
    plot(1:length(Data.saccadeSpikeInterval{i}), Data.saccadeSpikeInterval{i})
    spikeData{i}
    pause
end

%%
opt = ccm_options;
opt.multiUnit = true;

Data = ccm_interspike_interval('joule','jp125n04','spikeUnit17',opt);
%%

[trialData s] = load_data('joule', 'jp060n01', [mem_min_vars,'spikeData'], true);

outcome = {'saccToTarget'};
side = {'right'};
trials = mem_trial_selection(trialData, outcome, side);

%%
unitIndex = 1;
alignEvent = 'responseOnset';
alignList = trialData.(alignEvent)(trials);
 
    
% Get the rasters (and what index they align to)
[alignedRasters, alignmentIndex] = spike_to_raster(trialData.spikeData(trials, unitIndex), alignList);


%%

for i = 1 : size(classicDdmCancel, 1)
options = ccm_options;

options.multiUnit = true;
options.plotFlag = true;
options.printPlot = true;
options.ms2Std = 75;



Data = ccm_neuron_stop_vs_go('joule', classicDdmCancel.sessionID{i},  classicDdmCancel.unit(i), options);

options.unitArray = classicDdmCancel.unit(i);
options.doStops = false;

Data = ccm_session_data('joule', classicDdmCancel.sessionID{i},  options);
    
end
    
%%
opt = ccm_neuron_stop_vs_go;
opt.multiUnit = true;
opt.minTrialPerCond     = 10;
opt.plotFlag     = true;

% data = ccm_neuron_stop_vs_go('joule', 'jp125n04', {'spikeUnit26'}, opt);
% data = ccm_neuron_stop_vs_go('joule', 'jp125n04', {'spikeUnit32'}, opt);
data = ccm_neuron_stop_vs_go('broca', 'bp244n02', {'spikeUnit27'}, opt);
% data = ccm_neuron_stop_vs_go('broca', 'bp247n02', {'spikeUnit12'}, opt);
  



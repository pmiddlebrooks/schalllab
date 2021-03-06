%%
figure(10);
hold all
cellfun(@(x,y) plot(x,y, 'k'), goEyeX{1,1}, goEyeY{1,1})

% plot(goEyeX{1,1}, goEyeY{1,1})

degX = cellfun(@(x) sqrt(x), cellfun(@(x,y) x^2 + y^2, goEyeX{1,1}, goEyeY{1,1}))


velX{d, i} = cellfun(@(x) [0; diff(x(:))], degX, 'uni', false);


scatter(trialData.rt(goTrial{d, i}), cellfun(@max, goVel{d, i}))


%% Add hemisphere to translated data file
subject = 'broca';
session = {'bp095n04'};
% tebaPath = '/Volumes/SchallLab/data/';

hemisphere = 'right';

for i = 1 : length(session)
    [~, SessionData] = load_data(subject, session{i});
    
    SessionData.hemisphere = hemisphere;
    
    save(fullfile(local_data_path, subject, [session{i}, '.mat']), 'SessionData', '-append')
    %     save(fullfile(tebaPath, subject, [session{i}, '.mat']), 'SessionData', '-append')
    
    
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
for k = unitIndex
    kUnit = unitArrayNew{k}
    alignListTarg = trialData.targOn(iGoTrial);
    [alignedRasters, ~] = spike_to_raster(trialData.(kUnit)(iGoTrial), alignListTarg);
    iSdfTarg = nanmean(spike_density_function(alignedRasters, Kernel));
    iNormTarg = max(iSdfTarg);
    
    alignListSacc = trialData.responseOnset(iGoTrial);
    [alignedRasters, ~] = spike_to_raster(trialData.(kUnit)(iGoTrial), alignListSacc);
    iSdfSacc = nanmean(spike_density_function(alignedRasters, Kernel));
    iNormSacc = max(iSdfSacc);
    
    alignListChecker = trialData.targOn(iGoTrial);
    [alignedRasters, ~] = spike_to_raster(trialData.(kUnit)(iGoTrial), alignListChecker);
    iSdfChecker = nanmean(spike_density_function(alignedRasters, Kernel));
    iNormChecker = max(iSdfChecker);
    
    iNorm = max([iNormTarg, iNormSacc, iNormChecker]);
    
    [alignedRasters, alignmentIndex] = spike_to_raster(trialData.(kUnit)(iGoTrial), alignListGo);
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
[alignedRasters, alignmentIndex] = spike_to_raster(trialData.spikeUnit01a(trials), alignList);


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
opt.printPlot     = false;

% data = ccm_neuron_stop_vs_go('joule', 'jp125n04', {'spikeUnit26'}, opt);
% data = ccm_neuron_stop_vs_go('joule', 'jp106n02', {'spikeUnit18'}, opt);
data = ccm_neuron_stop_vs_go('broca', 'bp092n02', {'spikeUnit17'}, opt);
% data = ccm_neuron_stop_vs_go('broca', 'bp247n02', {'spikeUnit12'}, opt);


%%
%%
options = ccm_neuron_choice;
options.unitArray = {'spikeUnit03'};
options.multiUnit = true;

% unitInfo = ccm_neuron_choice('joule', 'jp054n02', options.unitArray, options);
unitInfo = ccm_neuron_choice('broca', 'bp229n02-mm', options.unitArray, options)

%% Establish options to send to ccm_session_data in the for loop below
close all

opt             = ccm_options;
opt.ms2Std = 75;
opt.howProcess  = 'print';
opt.plotFlag    = true;
opt.printPlot    = true;
opt.dataType    = 'neuron';
opt.collapseTarg 	= true;
opt.collapseSignal 	= true;
opt.doStops 	= false;



%     opt.unitArray = 'spikeData';
opt.unitArray = {'spikeUnit32a'};
opt.multiUnit = false;

opt.unitArray = {'spikeUnit32'};
opt.multiUnit = true;


subject = 'joule';
session = 'jp125n03';

%       iData = ccm_session_data(subject, session, opt);
opt.pairTriplet = 'pair';
Data = ccm_rt_history_neural(subject, session, opt)


%%
tic
unitArray = {'spikeUnit01a','spikeUnit01b','spikeUnit02a','spikeUnit02b','spikeUnit02c','spikeUnit03a','spikeUnit03b','spikeUnit04a','spikeUnit04b'};
trialData = load_data('joule','jp125n03',[ccm_min_vars,unitArray]);
fprintf('\nall single units: %.2f\n', toc)

%%
tic
trialData = load_data('joule','jp125n03',[ccm_min_vars]);
unitArray = {'spikeUnit01a','spikeUnit01b','spikeUnit02a','spikeUnit02b','spikeUnit02c','spikeUnit03a','spikeUnit03b','spikeUnit04a','spikeUnit04b'};
for i = 1 : length(unitArray)
    t = load_data('joule','jp125n03',[ccm_min_vars,unitArray(i)]);
    trialData.(unitArray{i}) = t.(unitArray{i});
end
fprintf('\nlooped single units: %.2f\n', toc)

%%
tic
unitArray = {'spikeUnit01','spikeUnit02','spikeUnit03','spikeUnit04'};
trialData = load_data('joule','jp125n03',[ccm_min_vars,unitArray], 1);
fprintf('\nall multi units: %.2f\n', toc)

%%
tic
trialData = load_data('joule','jp125n03',[ccm_min_vars]);
unitArray = {'spikeUnit01','spikeUnit02','spikeUnit03','spikeUnit04'};
for i = 1 : length(unitArray)
    t = load(fullfile(localDataPath, 'jp125n03'),unitArray{i}, 1);
    trialData.(unitArray{i}) = t.(unitArray{i});
end
fprintf('\nlooped multi units: %.2f\n', toc)


%%
close all
tic

options = ccm_options;

options.multiUnit = true;
options.plotFlag = true;
options.printPlot = true;
options.normalize = true;
opt.minTrialPerCond     = 10;
options.doStops = true;

options.unitArray = {'spikeUnit27'};

options.ANALYZE_CANCELED = true;
options.ANALYZE_NONCANCELED = true;
options.ms2Std = 75;
options.plotSigle = false;

Data = ccm_session_data('joule', 'jp125n04',  options);
% Data = ccm_neuron_stop_vs_go('joule', 'jp054n02',  options.unitArray, options);

toc

%%
options = ccm_neuron_choice;
options.unitArray = {'spikeUnit01'};
options.multiUnit = true;

% unitInfo = ccm_neuron_choice('joule', 'jp054n02', options.unitArray, options);
unitInfo = ccm_ddm_like('joule', 'jp054n02', options.unitArray, options);


%%
td = trialData;
var = fieldnames(trialData);
for i = 1 : length(var)
    if ~strncmp(var{i},'lfp', 3) &&...
            ~strncmp(var{i},'eye', 3) &&...
            ~strncmp(var{i},'targ1Check', 10) &&...
            isa(trialData.(var{i}), 'double')
        td.(var{i}) = int32(trialData.(var{i}));
    end
    if strncmp(var{i}, 'spike', 5)
        td.(var{i}) = cellfun(@(x) int32(x), td.(var{i}), 'uni', false);
    end
end



%%
j = 1;
subject = 'broca';
d(j).name = 'bp247n02.mat';
d(j).name = 'bp093n02.mat';

% noChangeList = {'checkerAmp'...
%     'checkerAngle'...
%     'checkerSize'...
%     'checkerWindow'...
%     'distAmp'...
%     'distAngle'...
%     'distSize'...
%     'distWindow'...
%     'eegData'...
%     'eyeX'...
%     'eyeY'...
%     'fixAmp'...
%     'fixAngle'...
%     'fixSize'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%     'xxxxxx'...
%

localDataPath = ['~/Dropbox/local_data/',lower(subject),'/'];

trialData = load(fullfile(localDataPath,d(j).name));
var = fieldnames(trialData);

for i = 1 : length(var)
    %     if ~strncmp(var{i},'lfp', 3) &&...
    %             ~strncmp(var{i},'eye', 3) &&...
    %             ~strncmp(var{i},'targ1Check', 10) &&...
    %             isa(trialData.(var{i}), 'double')
    %         td.(var{i}) = int32(trialData.(var{i}));
    %     end
    if strncmp(var{i}, 'spike', 5)
        var{i}
        trialData.(var{i}) = cellfun(@(x) int32(x), trialData.(var{i}), 'uni', false);
    end
end

% save(fullfile(local_data_path, subject, [d(i).name(1:end-4),'-new'), '-struct', 'trialData','-v7.3')

%%
tic
[trialData, SessionData, ExtraVariable] = load_data('broca', 'bp093n02', [ccm_min_vars,'spikeData']);
toc
%%
opt = ccm_options;
opt.multiUnit = true;
opt.unitArray = {'spikeUnit17'};
opt.printPlot = false;
opt.plotFlag = false;
tic
data = ccm_session_data('broca', 'bp093n02', opt);
clear data
toc
tic
data = ccm_session_data('broca', 'bp093n02', opt);
clear data
toc
tic
data = ccm_session_data('broca', 'bp093n02', opt);
clear data
toc
tic
data = ccm_session_data('broca', 'bp093n02', opt);
clear data
toc
tic
data = ccm_session_data('broca', 'bp093n02', opt);
clear data
toc

fprintf('\n\n\n\n Now New File:\n\n\n')

tic
data = ccm_session_data('broca', 'bp093n02-new', opt);
clear data
toc
tic
data = ccm_session_data('broca', 'bp093n02-new', opt);
clear data
toc
tic
data = ccm_session_data('broca', 'bp093n02-new', opt);
clear data
toc
tic
data = ccm_session_data('broca', 'bp093n02-new', opt);
clear data
toc
tic
data = ccm_session_data('broca', 'bp093n02-new', opt);
clear data
toc

%%
subject = 'joule';
tebaDataPath = '/Volumes/SchallLab/data/';

trialData = load(fullfile(tebaDataPath,subject,'jp083n02'));

[trialData, SessionData, ExtraVariable] = load_data('broca', 'bp093n02', [ccm_min_vars,'spikeData']);

%%
subject = 'xena'
localDataPath = ['~/Dropbox/local_data/',lower(subject),'/'];
tebaDataPath = '/Volumes/SchallLab/data/';

switch lower(subject)
    case 'joule'
        fileName = [subject, '.mat'];
        tebaDataPath = [tebaDataPath, 'Joule/'];
    case 'broca'
        fileName = [subject, '.mat'];
        tebaDataPath = [tebaDataPath, 'Broca/'];
    case 'xena'
        fileName = [subject, '.mat'];
        tebaDataPath = [tebaDataPath, 'Xena/Plexon/'];
    otherwise
        fprintf('%s is not a valid subject ID, try again?/n', subject)
        return
end

% d = dir(localDataPath);
d = dir(tebaDataPath);


for i = 1 : size(d, 1)
    i
    
    if ~isempty(regexp(d(i).name, '.*n0.*.mat')) && ~strncmp(d(i).name, '._', 2) && isempty(regexp(d(i).name, '.*legacy.*.mat'))
        tic
        disp(d(i).name)
        trialData = load(fullfile(tebaDataPath ,d(i).name));
        
        save(fullfile(local_data_path, subject, d(i).name(1:end-4)), '-struct', 'trialData')
        % Save to teba also
        save(fullfile(tebaDataPath, d(i).name(1:end-4)), '-struct', 'trialData')
        clear trialData
        
        toc
        
    end
    
end



%%    SET VARIABLES FOR RUNNING ANALYSES ON A SET OF UNITS UNDER SOME CATEGORY (SUCH AS PRESACC DDM-MODULATED CANCELING UNITS)

subject = 'broca';

multiUnit = true;
ssrtUse = 'intWeightPerSession';

% Which ddm criteria do we want to use?
ddmType = 'ddmRankMean';  % ddm coherence determined by comparing means of easy vs hard spike rates into RF (and out)
% ddmType = 'ddm';  % ddm coherence determined by ranksum test and slopes of means easy vs hard rates into RF (and out)

% Which cancel time criterium do we want to use?
cancelType = 'meanDifference';
cancelType = 'trialByTrial';
cancelType = 'meanSdf';

saccadeBaseRatio = 2;
saccadeBaseRatio = [];

deleteUnmodulated = true;
deleteSessions = true;

if multiUnit
    addMulti = '_multiUnit';
else
    addMulti = [];
end


%% SET PATHS

projectDate = '2017-01-11';
accreRoot = '/gpfs22';
accreHome = '/home/middlepg';
accreScratch = '/scratch/middlepg';
if isdir(fullfile(accreScratch))
    matRoot = fullfile(accreRoot,accreHome,'m-files'); % Edit this if running on Accre
    projectRoot = fullfile(accreScratch,'perceptualchoice_stop_model');
    environ = 'accre';
else
    matRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/schalllab';
    projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';
    matRoot = '~/schalllab';
    projectRoot = '~/perceptualchoice_stop_spikes_population';
    environ = 'local';
end

addpath(genpath(fullfile(matRoot,'ccm')));
addpath(genpath(fullfile(projectRoot,'src/code',projectDate)));

cd(matRoot)


projectRoot = '~/perceptualchoice_stop_spikes_population';

addpath(genpath(fullfile(projectRoot,'src/code',projectDate)));
dataPath = fullfile(projectRoot,'data',projectDate,subject);


%%     OPEN THE RELEVANT CATEGORY

dataPath = fullfile(projectRoot,'data',projectDate,subject);


% Possible categories:
% ddm cancel
category = 'ccm_presacc_ddmRankMeanStim_cancel_meanSdf_neurons_multiUnit';

% ddm no-cancel
% category = 'ccm_presacc_ddmRankMeanStim_noCancel_meanSdf_neurons_multiUnit';

% no-ddm cancel
% category = 'ccm_presacc_cancel_meanSdf_noddmRankMeanStim_neurons_multiUnit';

% no-ddm no-cancel
% category = 'ccm_presacc_noddmRankMeanStim_noCancel_meanSdf_neurons_multiUnit';

% visual-only neurons
% category = 'ccm_visNoPresacc_neurons_multiUnit';

load(fullfile(dataPath, category))
if strcmp(category, 'ccm_visNoPresacc_neurons_multiUnit')
    neurons = neurons(strcmp(neurons.rf, 'none'),:);
end
neurons
size(neurons)
%% Loop through categories and create ccm_neuron_stop_vs_go figures
subject = 'broca';

dataPath = fullfile(projectRoot,'data',projectDate,subject);
categoryList = {'ccm_presacc_ddmRankMeanStim_cancel_meanSdf_neurons_multiUnit', ...
    'ccm_presacc_ddmRankMeanStim_noCancel_meanSdf_neurons_multiUnit', ...
    'ccm_presacc_cancel_meanSdf_noddmRankMeanStim_neurons_multiUnit', ...
    'ccm_presacc_noddmRankMeanStim_noCancel_meanSdf_neurons_multiUnit'};
options = ccm_neuron_stop_vs_go;
% options = ccm_options;

options.multiUnit = true;
options.plotFlag = true;
options.printPlot = true;
options.plotSingle = true;
options.ANALYZE_CANCELED = true;
options.ANALYZE_NONCANCELED = true;
options.normalize = true;

for c = 1 : length(categoryList)
    options.category = categoryList{c};
    load(fullfile(dataPath, categoryList{c}))
    
    
    for i = 1 : size(neurons, 1)
        %     for i = size(neurons, 1) : size(neurons, 1)
        fprintf('\n%d of %d category: %s\n', i, size(neurons, 1), categoryList{c})
        Data = ccm_neuron_stop_vs_go(subject, neurons.sessionID{i},  neurons.unit(i), options);
    end
    
end
%% Add hemisphere to translated data file
session = {'bp234n02'};
% tebaPath = '/Volumes/SchallLab/data/';

hemisphere = 'left';

for i = 1 : length(session)
    [~, SessionData] = load_data(subject, session{i});
    
    SessionData.hemisphere = hemisphere;
    
    save(fullfile(local_data_path, subject, [sessionID, '.mat']), 'SessionData', '-append')
    %     save(fullfile(tebaPath, subject, [session{i}, '.mat']), 'SessionData', '-append')
    
    
end


%%
plexon_translate_datafile_mac('broca','bp247n02');

%%
subject = 'broca';
session = 'bp093n02';
unitArray = {'spikeUnit17'};
session = 'bp246n02';
unitArray = {'spikeUnit30'};


%%
options = ccm_neuron_stop_vs_go;
% options = ccm_options;

options.normalize = true;
options.multiUnit = true;
options.plotFlag = true;
options.printPlot = false;
options.ANALYZE_NONCANCELED = true;
options.ANALYZE_CANCELED = true;
options.plotSingle          = true;

Data = ccm_neuron_stop_vs_go(subject, session, unitArray, options);


%%
subject = 'joule';
dataPath = fullfile(projectRoot,'data',projectDate,subject);


% Possible categories:
% ddm cancel
category = 'ccm_presacc_ddmRankMeanStim_cancel_meanSdf_neurons_multiUnit';

% cancel no-ddm
% category = 'ccm_presacc_cancel_meanSdf_noddmRankMeanStim_neurons_multiUnit';


load(fullfile(dataPath, category))

optInh              = ccm_options;
optInh.plotFlag     = false;
oldSSRT = nan(size(neurons, 1), 1);
newSSRT = nan(size(neurons, 1), 1);

optInh.INCLUDE_GO_OMISSION = false;
for i = 1 : size(neurons, 1)
    dataInh             = ccm_inhibition(subject, neurons.sessionID{i}, optInh);
    oldSSRT(i) = round(nanmean(dataInh.ssrtIntegrationWeighted));
end

optInh.INCLUDE_GO_OMISSION = true;
for i = 1 : size(neurons, 1)
    dataInh             = ccm_inhibition(subject, neurons.sessionID{i}, optInh);
    newSSRT(i) = round(nanmean(dataInh.ssrtIntegrationWeighted));
end

neurons.oldSSRT = oldSSRT;
neurons.newSSRT = newSSRT;
neurons.ssrtDiff = newSSRT - oldSSRT;

neurons
%%
%% Inhibition and chronometric population plots
subject = 'broca';
dataPath = fullfile(projectRoot,'data',projectDate,subject);
category = 'presacc';
load(fullfile(dataPath, ['ccm_',category,'_neurons']))
sessionSet = unique(neurons.sessionID);

dataInh = ccm_inhibition_population(subject, sessionSet(1:2));

%% Why are there so many that now don't count as cancel?
options = ccm_neuron_stop_vs_go;
% options = ccm_options;

options.multiUnit = true;
options.plotFlag = true;
options.printPlot = true;
options.plotSingle = true;
options.ANALYZE_CANCELED = true;
options.ANALYZE_NONCANCELED = false;


for i = 1 : size(neuronsNotDC, 1)
    fprintf('\n%d of %d \n', i, size(neuronsNotDC, 1))
    Data = ccm_neuron_stop_vs_go(subject, neuronsNotDC.sessionID{i},  neuronsNotDC.unit(i), options);
end
disp('done')


%% Test whether ccm_find_saccade_rf really needs collapss color coherence...

category = 'ccm_presacc_neurons_multiUnit';
subject = 'joule';
dataPath = fullfile(projectRoot,'data',projectDate,subject);
load(fullfile(dataPath, category))

opt             = ccm_options;
opt.multiUnit    = options.multiUnit;
opt.plotFlag    = false;
opt.printFlag    = false;
opt.collapseTarg = true;
opt.doStops   = false;

rf = cell(size(neurons, 1), 1);
rfCollapse = cell(size(neurons, 1), 1);

for i = 1 : size(neurons, 1)
    fprintf('\n%d of %d \n', i, size(neurons, 1))
    
    opt.unitArray   = neurons.unit(i);
    
    opt.collapseSignal = true;
    UnitCollapse                = ccm_session_data(subject, neurons.sessionID{i}, opt);
    rfCollapse{i} = ccm_find_saccade_rf(UnitCollapse);
    
    opt.collapseSignal = false;
    UnitCollapse                = ccm_session_data(subject, neurons.sessionID{i}, opt);
    rf{i} = ccm_find_saccade_rf(UnitCollapse);
    
end
%%
[rf, rfCollapse]
sum(~strcmp(rf, rfCollapse))

%%
subject = 'broca';
session = 'bp247n02';
    [~, SessionData] = load_data(subject, session);



function batch

%%
%%
sessionList = ...
    {'bp229n02-mm', ...
    'bp203n02',...
    'bp204n02',...
    'bp205n02',...
    'bp206n02'};
sessionList = ...
    {'bp198n01', ...
    'bp199n01',...
    'bp204n02',...
    'bp205n02',...
    'bp206n02'};
for i = 1 : length(sessionList)
    iSessionID = sessionList{i}
    plexon_translate_datafile_mac('broca',iSessionID);
end
%%      GET FIGURES FOR CCM: CCM_SESSION_DATA, CCM_DDM_LIKE, AND CCM_NEURON_STOP_VS_GO

% dataPath = 'Broca_sessions.mat';
% load(dataPath)
% nSession = length(sessions.ccm.stop);

subjectID = 'broca';
taskID = 'ccm';
dataIncludeArray = {'spike'};

sessionList = find_sessions_data(subjectID, taskID, dataIncludeArray);


opt = ccm_options;
opt.collapseSignal = true;
% sessionList = {'bp093n02'};
for i = 1 : length(sessionList)
    iSessionID = sessionList{i};
    data = ccm_session_data('broca', iSessionID, opt);
end

%%
dataDir = '/Volumes/SchallLab/data/Broca';
d = dir(dataDir);
for i = 1240 : size(d, 1)
    if regexp(d(i).name, '.*n01.plx')
        disp(d(i).name(1:end-4))
        plexon_translate_datafile_mac('broca',d(i).name(1:end-4));
    end
end

%%
dataDir = '/Volumes/SchallLab/data/Joule';
d = dir(dataDir);
for i = 100 : size(d, 1)
    if regexp(d(i).name, 'jp.*.plx')
        disp(d(i).name(1:end-4))
        plexon_translate_datafile_mac('joule',d(i).name(1:end-4));
    end
end


%%      GET FIGURES FOR MEM/DEL: MEM_SESSION_DATA

% dataPath = 'Broca_sessions.mat';
% load(dataPath)
% nSession = length(sessions.ccm.stop);

subjectID = 'broca';
taskID = 'mem';
taskID = 'del';
dataIncludeArray = {'spike'};

sessionList = find_sessions_data(subjectID, taskID, dataIncludeArray);


% opt = ccm_options;
% opt.collapseSignal = true;
% sessionList = {'bp093n02'};
for i = 1 : length(sessionList)
    iSessionID = sessionList{i};
    data = mem_session_data('broca', iSessionID);
end

%%
for i = 1 : length(sessionList)
    iSessionID = sessionList{i}
    [data, options] = ccm_session_data(subjectID, iSessionID,'neuron', 'normalize', false, 'filterData', true,'collapseSignal', false, 'printPlot', true);
    options.collapseSignal = true;
    ccm_session_data_plot(data, 'neuron', options);
    ccm_neuron_stop_vs_go(subjectID, iSessionID, 'collapseSignal', 1, 'printPlot', 1);
    ddmLike = ccm_ddm_like(subjectID, iSessionID, 'printPlot', 1);
end


%%
sessionArray = ...
    {'bp076n01', ...
    'bp077n01', ...
    'bp080n01', ...
    'bp081n01', ...
    'bp082n02', ...
    'bp084n02', ...
    'bp086n01', ...
    'bp087n02', ...
    'bp088n02', ...
    'bp089n02', ...
    'bp090n02', ...
    'bp091n02-pm', ...
    'bp092n02', ...
    'bp093n02', ...
    'bp095n04', ...
    'bp096n02', ...
    'bp097n03', ...
    'bp099n02-pm', ...
    'bp100n01', ...
    'bp101n01', ...
    'bp103n01', ...
    'bp104n01', ...
    'bp106n01', ...
    'bp107n01', ...
    'bp109n01', ...
    'bp110n01-pm', ...
    'bp111n02-pm', ...
    'bp112n02-pm', ...
    'bp113n01', ...
    'bp114n01', ...
    'bp115n02', ...
    'bp116n02-pm', ...
    'bp117n01', ...
    'bp118n01', ...
    'bp119n02-pm', ...
    'bp120n02-pm', ...
    'bp121n02-pm', ...
    'bp121n04', ...
    'bp122n02', ...
    'bp123n01', ...
    'bp124n02', ...
    'bp124n04', ...
    'bp126n02-pm', ...
    'bp127n02', ...
    'bp128n02', ...
    'bp129n02-pm', ...
    'bp130n04', ...
    'bp131n02', ...
    'bp132n02-pm'};


for i = length(sessionArray)
    iSessionID = sessionArray{i}
    
    [td, SessionData, ExtraVar] = load_data('broca', iSessionID);
    disp([ExtraVar.pSignalArray, unique(sum(td.checkerArray,2))])
    pause
end


%%
sessionArray = ...
    {'andy_fef_086', ...
    'andy_fef_096', ...
    'andy_fef_196', ...
    'andy_fef_215', ...
    'andy_fef_218', ...
    'andy_fef_219', ...
    'andy_fef_220', ...
    'andy_fef_221', ...
    'andy_fef_226', ...
    'andy_fef_228', ...
    'andy_fef_234', ...
    'andy_fef_240'};

subjectID = 'andy';
ssrt = nan(length(sessionArray), 1);
ssdArray  = [];
nTrial = nan(length(sessionArray), 1);
for i = 3 : length(sessionArray)
    sessionArray{i}
    D(i).data = cmd_inhibition(subjectID, sessionArray{i});
    
    ssrt(i) = D(i).data.ssrt.integrationWeighted;
    ssdArray = [ssdArray; D(i).data.ssdArray];
    D(i).data.nTrial
    nTrial(i) = D(i).data.nTrial;
    
    pause
    data = cmd_session_data(subjectID, sessionArray{i});
    pause
    cmd_neuron_stop_vs_go(subjectID, sessionArray{i});
    pause
end
%%
sessionArray = ...
    {'andy_fef_084', ...
    'andy_fef_096', ...
    'andy_fef_154', ...
    'andy_fef_196', ...
    'andy_fef_215', ...
    'andy_fef_218', ...
    'andy_fef_219', ...
    'andy_fef_220', ...
    'andy_fef_226', ...
    'andy_fef_228', ...
    'andy_fef_234', ...
    'andy_fef_240'};

subjectID = 'andy';
ssrt = nan(length(sessionArray), 1);
ssdArray  = [];
nTrial = nan(length(sessionArray), 1);
for i = 1 : length(sessionArray)
    D(i).data = cmd_inhibition(subjectID, sessionArray{i}, 'plotFlag', false)
    
    ssrt(i) = D(i).data.ssrt.integrationWeighted
    ssdArray = [ssdArray; D(i).data.ssdArray]
    nTrial(i) = D(i).data.nTrial
end
%%
sessionArray = ...
    {'chase_fef_375', ...
    'chase_fef_387', ...
    'chase_fef_389', ...
    'chase_fef_395', ...
    'chase_fef_413'};

subjectID = 'chase';
ssrt = nan(length(sessionArray), 1);
ssdArray  = [];
nTrial = nan(length(sessionArray), 1);
for i = 1 : length(sessionArray)
    D(i).data = cmd_inhibition(subjectID, sessionArray{i})
    %    D(i).data = cmd_inhibition(subjectID, sessionArray{i}, 'plotFlag', false)
    pause
    ssrt(i) = D(i).data.ssrt.integrationWeighted
    ssdArray = [ssdArray; D(i).data.ssdArray]
    nTrial(i) = D(i).data.nTrial
end

%% MEMORY GUIDED SACCADE ANALYSES


sessionArray = ...
    {'bp082n01', ...
    'bp084n01', ...
    'bp085n01', ...
    'bp088n01', ...
    'bp089n01', ...
    'bp090n01', ...
    'bp091n01', ...
    'bp092n01', ...
    'bp093n01', ...
    'bp096n01', ...
    'bp097n01', ...
    'bp099n01', ...
    'bp111n01-pm', ...
    'bp112n01', ...
    'bp115n01', ...
    'bp116n01', ...
    'bp119n01', ...
    'bp120n01', ...
    'bp121n01', ...
    'bp121n03', ...
    'bp124n01', ...
    'bp124n03', ...
    'bp126n01', ...
    'bp128n01', ...
    'bp129n01', ...
    'bp130n03-pm', ...
    'bp131n01', ...
    'bp132n01', ...
    'bp156n01', ...
    'bp158n01', ...
    'bp159n01', ...
    'bp159n02', ...
    'bp162n01', ...
    'bp163n01', ...
    'bp164n02', ...
    'bp166n01', ...
    'bp167n01', ...
    'bp168n01', ...
    'bp192n01', ...
    'bp193n01', ...
    'bp195n01', ...
    'bp195n03', ...
    'bp196n01', ...
    'bp197n01'};

sessionArray = ...
    {'bp195n01', ...
    'bp195n03', ...
    'bp196n01', ...
    'bp197n01'};



opt = mem_session_data;
opt.printPlot = true;
opt.plotFlag = true;

om = mem_plot_epoch;
om.printPlot = true;
for i = 1 : length(sessionArray)
    iSessionID = sessionArray{i}
    
    U = mem_session_data('broca',iSessionID,opt);
    
    
    %    D = mem_plot_epoch(U,om);
    
end

%%
% 85, 87, 89, 90, 91, 92, 93, 95, 96, 97, 99, 111, 112, 115, 116, 119, 120, 121, 124, 126, 128, 129, 131, 132,
% 151, 156, 158, 162, 163, 166, 167, 168, 192, 194, 195, 196, 197
%
% why not 88, 110, (121 is twice), 130, 160, 164
%



sessionArray = ...
    {'bp085n02', ...
    'bp087n02', ...
    'bp088n02', ...
    'bp089n02', ...
    'bp090n02', ...
    'bp091n02', ...
    'bp092n02-pm', ...
    'bp093n02', ...
    'bp096n02', ...
    'bp097n02', ...
    'bp097n03', ...
    'bp099n02-pm', ...
    'bp110n01', ...
    'bp111n02-pm', ...
    'bp112n02-pm', ...
    'bp115n02', ...
    'bp116n02-pm', ...
    'bp119n02-pm', ...
    'bp120n02-pm', ...
    'bp121n02-pm', ...
    'bp121n04', ...
    'bp124n02', ...
    'bp126n02-pm', ...
    'bp128n02', ...
    'bp129n02-pm', ...
    'bp130n04-pm', ...
    'bp131n02', ...
    'bp132n02-pm', ...
    'bp151n02', ...
    'bp156n02', ...
    'bp158n02-pm', ...
    'bp160n01', ...
    'bp162n03-pm', ...
    'bp163n02-01', ...
    'bp164n01-pm', ...
    'bp166n02', ...
    'bp167n04', ...
    'bp168n04', ...
    'bp192n02-tp', ...
    'bp193n02-tp', ...
    'bp194n02', ...
    'bp195n04-tp', ...
    'bp196n02-tp', ...
    'bp197n02-tp'};

validSession = {};
nGo = {};
psych = {};
for i = 1 : length(sessionArray)
    i
    iSessionID = sessionArray{i}
    
    [dataFile, localDataPath, localDataFile] = data_file_path('broca', iSessionID, 'monkey');
    
    % If the file hasn't already been copied to a local directory, do it now
    if exist(dataFile, 'file') ~= 2
        plexon_translate_datafile_mac('broca', iSessionID)
    end
    
    U = ccm_session_behavior('broca',iSessionID,'plotFlag',false);
    
    U.nGo
    
    
    criteria = all(U.nGo > 60);
    if criteria
        validSession = [validSession; iSessionID];
        nGo = [nGo; {U.nGo}];
        psych = [psych; {U.nGoRight ./ U.nGo}];
    end
end

%%
figure()
for i = 1 : length(validSession)
    clf
    plot(1:length(psych{i}), psych{i})
    set(gca, 'ylim', [0 1])
    title(sprintf('%s', validSession{i}))
    pause
end

%%
sessionArray = ...
    {'bp216n02', ...
    'bp218n02', ...
    'bp220n02', ...
    'bp221n02', ...
    'bp222n02', ...
    'bp224n02', ...
    'bp226n02', ...
    'bp227n02', ...
    'bp228n02', ...
    'bp229n02', ...
    'bp230n02'};

opt = ccm_session_data;
opt.dataType = 'lfp';
for i = 1 : length(sessionArray)
    i
    iSessionID = sessionArray{i}
    
    [dataFile, localDataPath, localDataFile] = data_file_path('broca', iSessionID, 'monkey');
    
    % If the file hasn't already been translated to teba, do it now
    if exist(dataFile, 'file') ~= 2
        plexon_translate_datafile_mac('broca', iSessionID)
    end
    
    data = ccm_session_data('broca',iSessionID);
    dataOpt = ccm_session_data('broca',iSessionID,opt);
    clear('data','dataOpt')
end


%%
sessionArray = ...
    {'bp216n02', ...
    'bp218n02', ...
    'bp220n02', ...
    'bp221n02', ...
    'bp222n02', ...
    'bp224n02', ...
    'bp226n02', ...
    'bp227n02', ...
    'bp228n02', ...
    'bp229n02', ...
    'bp230n02'};
validSession = {};
nGo = {};
psych = {};
for i = 9 : length(sessionArray)
    i
    iSessionID = sessionArray{i}
    
    [dataFile, localDataPath, localDataFile] = data_file_path('broca', iSessionID, 'monkey');
    
    % If the file hasn't already been translate, do it now
    if exist(dataFile, 'file') ~= 2
        plexon_translate_datafile_mac('broca', iSessionID)
    end
    
    U = ccm_session_behavior('broca',iSessionID,'plotFlag',true);
    
    U.nGo
    
    
    criteria = all(U.nGo > 60);
    if criteria
        validSession = [validSession; iSessionID];
        nGo = [nGo; {U.nGo}];
        psych = [psych; {U.nGoRight ./ U.nGo}];
    end
    
    pause
    clear U
end





%% print joule's session behavior for each ccm session
subject = 'joule';

projectDate = '2016-08-12';
projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';

dataPath = fullfile(projectRoot,'data',projectDate,subject);


% Open the sessions file and makes lists of the entries
fid=  fopen(fullfile(dataPath,['ccm_sessions_',subject,'.csv']));


nCol = 5;
formatSpec = '%s';
mHeader = textscan(fid, formatSpec, nCol, 'Delimiter', ',');

mData = textscan(fid, '%s %s %d %d %d', 'Delimiter', ',','TreatAsEmpty',{'NA','na'});

sessionList     = mData{1};
hemisphereList  = mData{2};

for i = 70 : length(sessionList)
    data = ccm_session_behavior(subject, sessionList{i});
    clear data
end

%% print broca's session behavior for each ccm session
subject = 'broca';

projectDate = '2016-08-12';
projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';

dataPath = fullfile(projectRoot,'data',projectDate,subject);


% Open the sessions file and makes lists of the entries
fid=  fopen(fullfile(dataPath,['ccm_sessions_',subject,'.csv']));


nCol = 5;
formatSpec = '%s';
mHeader = textscan(fid, formatSpec, nCol, 'Delimiter', ',');

mData = textscan(fid, '%s %s %d %d %d', 'Delimiter', ',','TreatAsEmpty',{'NA','na'});

sessionList     = mData{1};
hemisphereList  = mData{2};

for i = 1 : length(sessionList)
    data = ccm_session_behavior(subject, sessionList{i});
    clear data
end

%% List number of neurons in each category for broca and joule
subject = 'joule';
% subject = 'joule';

projectDate = '2016-08-12';
projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';
dataPath = fullfile(projectRoot,'data',projectDate,subject);

neuronNumber = table();

% Total channels recorded, total modulation neurons
load(fullfile(dataPath, 'ccm_neuronTypes'))
neuronNumber.sessions = length(unique(neuronTypes.sessionID));
neuronNumber.channels = size(neuronTypes, 1);
neuronNumber.modulated = sum(neuronTypes.fix | ...
    neuronTypes.vis | ...
    neuronTypes.checker | ...
    neuronTypes.presacc | ...
    neuronTypes.postsacc | ...
    neuronTypes.reward | ...
    neuronTypes.intertrial);

categories = {'fix', 'visNoPresacc', 'visPresacc', 'presaccNoVis', 'postSaccNoPresacc'};
for i = 1 : length(categories)
    load(fullfile(dataPath, ['ccm_',categories{i},'_neurons']))
    neuronNumber.(categories{i}) = length(neurons.sessionID);
end


%% Inhibition and chronometric population plots
subject = 'broca';
dataPath = fullfile(projectRoot,'data',projectDate,subject);
category = 'presacc';
load(fullfile(dataPath, ['ccm_',category,'_neurons']))
sessionSet = unique(neurons.sessionID);

dataInh = ccm_inhibition_population(subject, sessionSet);
% dataChron =  ccm_chronometric_population(subject, sessionSet);
% ccm_psychometric_population(subject, sessionSet);

%% Population behavioral measures
% subject = 'joule';
subject = 'broca';

% category = 'vis';
% load(fullfile(dataPath, ['ccm_',category,'_neurons']))
% sessionSet = unique(neurons.sessionID);


sessionSet = 'behavior1';
% ccm_chronometric_population(subject, sessionSet);
% ccm_psychometric_population(subject, sessionSet);
% ccm_rt_distribution_population(subject, sessionSet);
ccm_inhibition_population(subject, sessionSet);

%%
subject = 'broca';
dataPath = fullfile(projectRoot,'data',projectDate,subject);
category = 'ddm';
load(fullfile(dataPath, ['ccm_',category,'_neurons']))

% interest = table();

opt = ccm_options;

for i = 1 : length(neurons.sessionID)
    fprintf('%d\tSession: %s\t Unit: %s\n', i, neurons.sessionID{i}, neurons.unit{i})
    iTable = table();
    iTable.session  = neurons.sessionID(i);
    iTable.unit  = neurons.unit(i);
    % opt.unitArray = neurons.unit(i);
    data = ccm_neuron_stop_vs_go(subject, neurons.sessionID{i}, neurons.unit(i));
    prompt = 'add to list?';
    addToList = input(prompt);
    if addToList
        
        interest = [interest; iTable];
        
    end
    clear data
end

%% Total number of sessions, trial numbers for all recorded sessions
subject = 'broca';
% subject = 'joule';

projectDate = '2016-08-12';
projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';
dataPath = fullfile(projectRoot,'data',projectDate,subject);

number = table();

% Total channels recorded, total modulation neurons
load(fullfile(dataPath, 'ccm_neuronTypes'))
sessions = unique(neuronTypes.sessionID);
for i = 1 : length(sessions)
    iNumber = table();
    iNumber.session = sessions(i);
    [td, S] = load_data(subject, sessions{i});
    opt = ccm_options;
    
    % Go trials
    opt.outcome = {...
        'goCorrectTarget', 'goCorrectDistractor', ...
        };
    
    iNumber.goTrials = length(ccm_trial_selection(td, opt));
    % Stop trials
    opt.outcome = {...
        'stopCorrect', ...
        'stopIncorrectTarget', 'stopIncorrectDistractor'};
    
    iNumber.stopTrials = length(ccm_trial_selection(td, opt));
    number = [number; iNumber];
    
end
disp('done')
%% Total number of sessions, trial numbers for all recorded sessions
subject = 'broca';
% subject = 'joule';

projectDate = '2016-08-12';
projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';
dataPath = fullfile(projectRoot,'data',projectDate,subject);

number = table();

% Total channels recorded, total modulation neurons
load(fullfile(dataPath, 'ccm_neuronTypes'))
sessions = unique(neuronTypes.sessionID);
for i = 1 : length(sessions)
    iNumber = table();
    iNumber.session = sessions(i);
    [td, S] = load_data(subject, sessions{i});
    opt = ccm_options;
    
    % Go trials
    opt.outcome = {...
        'goCorrectTarget', 'goCorrectDistractor', ...
        };
    
    iNumber.goTrials = length(ccm_trial_selection(td, opt));
    % Stop trials
    opt.outcome = {...
        'stopCorrect', ...
        'stopIncorrectTarget', 'stopIncorrectDistractor'};
    
    iNumber.stopTrials = length(ccm_trial_selection(td, opt));
    number = [number; iNumber];
    
end
disp('done')
%% Total number of sessions, trial numbers for all recorded sessions. 32 channel neuronexus, u-probe, etc
subject = 'broca';
% subject = 'joule';

projectDate = '2016-08-12';
projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';
dataPath = fullfile(projectRoot,'data',projectDate,subject);

number = table();

% Total channels recorded, total modulation neurons
load(fullfile(dataPath, 'ccm_neuronTypes'))
sessions = unique(neuronTypes.sessionID);
for i = 1 : length(sessions)
    iNumber = table();
    
    iNumber.session = sessions(i);
    % find the index of the first spike unit from this session
    firstInd = find(strcmp(sessions{i}, neuronTypes.sessionID), 1, 'first');
    lastInd = find(strcmp(sessions{i}, neuronTypes.sessionID), 1, 'last');
    
    % How many spike units recorded on each channel?
    unitInd = cellfun(@(x) str2double(regexp(x,'\d*','Match')), neuronTypes.unit(firstInd:lastInd));
    
    % How many channels/electrodes were recorded/used during this session?
    iChannels = unique(unitInd);
    iNumber.channels = length(iChannels);
    
    iNumber.unitPerChannel = cell(1,1);
    iNumber.unitPerChannel{1} = nan(1, iNumber.channels);
    
    for j = 1 : iNumber.channels
        iNumber.unitPerChannel{1}(j) = sum(iChannels(j) == unitInd);
    end
    
    iNumber.unitPerSession = sum(iNumber.unitPerChannel{1});
    number = [number; iNumber];
    
end
disp(number)

channelType = unique(number.channels);
electrode = nan(length(channelType),1);
meanSignal = nan(length(channelType),1);
for i = 1 : length(channelType)
    
    % How many sessions with this electrode configuration?
    iChannelInd = number.channels == channelType(i);
    electrode(i) = sum(iChannelInd);
    meanSignal(i) = mean(number.unitPerSession(iChannelInd));
    
    fprintf('%d channel:\t %d sessions\t Mean single/multi unit per channel: %.1f\n', channelType(i), electrode(i), meanSignal(i))
end
%%
cEpoch = 'presacc';
presaccSessions = c(c.(cEpoch) & ~dg.ddm, :);
opt = ccm_options;
presaccNotDDM = table();

for i = 1 : size(presaccSessions, 1)
    iInterest = table;
    iInterest.sessionID = presaccSessions.sessionID(i);
    iInterest.unit = presaccSessions.unit(i);
    
    fprintf('%s\tRF:%s\n', presaccSessions.sessionID{i}, presaccSessions.rf{i})
    pdfName1 = [presaccSessions.sessionID{i},'_ccm_',presaccSessions.unit{i},'_neuron.pdf'];
    %     pdfName2 = [presaccSessions.sessionID{i},'_',presaccSessions.unit{i},'_ccm_neuron.pdf'];
    if exist(fullfile(local_figure_path,subject,pdfName1))
        open(fullfile(local_figure_path,subject,pdfName1))
    elseif exist(fullfile(local_figure_path,subject,pdfName2))
        open(fullfile(local_figure_path,subject,pdfName2))
    else
        opt.unitArray = presaccSessions.unit(i);
        ccm_session_data(subject, presaccSessions.sessionID{i}, opt);
    end
    
    
    prompt = 'add to list?';
    addToList = input(prompt);
    if addToList
        
        presaccNotDDM = [presaccNotDDM; iInterest];
        
    end
end

%%
opt = ccm_options;
presaccNotDDM = table();
for i = 1 : size(presaccList, 1)
    iInterest = table;
    iInterest.sessionID = presaccList.sessionID(i);
    iInterest.unit = presaccList.unit(i);
    
    fprintf('%s\t%s\n', presaccList.sessionID{i}, presaccList.unit{i})
    opt.unitArray = presaccList.unit(i);
    pdfName1 = [presaccList.sessionID{i},'_ccm_',presaccList.unit{i},'_neuron.pdf'];
    %     pdfName2 = [presaccList.sessionID{i},'_',presaccList.unit{i},'_ccm_neuron.pdf'];
    if exist(fullfile(local_figure_path,subject,pdfName1))
        open(fullfile(local_figure_path,subject,pdfName1))
    elseif exist(fullfile(local_figure_path,subject,pdfName2))
        open(fullfile(local_figure_path,subject,pdfName2))
    else
        opt.unitArray = presaccList.unit(i);
        ccm_session_data(subject, presaccList.sessionID{i}, opt);
    end
    
    prompt = 'add to list?';
    addToList = input(prompt);
    if addToList
        
        presaccNotDDM = [presaccNotDDM; iInterest];
        
    end
    
end

%%

opt = ccm_options;
ddmNoPresacc = table();
for i = 1 : size(dg, 1)
    iInterest = table;
    iInterest.sessionID = dg.sessionID(i);
    iInterest.unit = dg.unit(i);
    
    fprintf('%s\t%s\n', dg.sessionID{i}, dg.unit{i})
    opt.unitArray = dg.unit(i);
    pdfName1 = [dg.sessionID{i},'_ccm_',dg.unit{i},'_neuron.pdf'];
    %     pdfName2 = [dg.sessionID{i},'_',dg.unit{i},'_ccm_neuron.pdf'];
    if exist(fullfile(local_figure_path,subject,pdfName1))
        open(fullfile(local_figure_path,subject,pdfName1))
    elseif exist(fullfile(local_figure_path,subject,pdfName2))
        open(fullfile(local_figure_path,subject,pdfName2))
    else
        opt.unitArray = dg.unit(i);
        ccm_session_data(subject, dg.sessionID{i}, opt);
    end
    
    prompt = 'add to list?';
    addToList = input(prompt);
    if addToList
        
        ddmNoPresacc = [ddmNoPresacc; iInterest];
        
    end
    
end

%%
opt = ccm_options;
for i = 1 : size(ddmNoPresacc, 1)
    pdfName1 = [ddmNoPresacc.sessionID{i},'_ccm_',ddmNoPresacc.unit{i},'_neuron.pdf'];
    %     pdfName2 = [presaccList.sessionID{i},'_',presaccList.unit{i},'_ccm_neuron.pdf'];
    if exist(fullfile(local_figure_path,subject,pdfName1))
        open(fullfile(local_figure_path,subject,pdfName1))
    else
        opt.unitArray = presaccList.unit(i);
        ccm_session_data(subject, presaccList.sessionID{i}, opt);
    end
    pause
end









%% Total number of sessions, trial numbers for all recorded sessions
subject = 'broca';
% subject = 'joule';

projectDate = '2016-08-12';
projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';
dataPath = fullfile(projectRoot,'data',projectDate,subject);


% Total channels recorded, total modulation neurons
load(fullfile(dataPath, 'ccm_neuronTypes'))
sessions = unique(neuronTypes.sessionID);
for i = 1 : length(sessions)
    
    [td, S] = load_data(subject,sessions{i});
    nUnit = length(S.spikeUnitArray);
    
    for j = 1 : nUnit
        jUnitName = S.spikeUnitArray{j};
        saveFileName = [sessions{i}, '_', jUnitName];
        
        spikeData = td.spikeData(:, j);
        save(fullfile(local_data_path, subject, saveFileName), 'spikeData')
    end
    
end

%% Correlate SSRT with mean RT per session
subject = 'broca';
dataPath = fullfile(projectRoot,'data',projectDate,subject);

% Open the sessions file and makes lists of the entries
fid=  fopen(fullfile(dataPath,['ccm_sessions_',subject,'.csv']));


nCol = 5;
formatSpec = '%s';
mHeader = textscan(fid, formatSpec, nCol, 'Delimiter', ',');

mData = textscan(fid, '%s %s %d %d %d', 'Delimiter', ',','TreatAsEmpty',{'NA','na'});

sessionList     = mData{1};

opt = ccm_options;
opt.plotFlag = false;
opt.printFlag = false;

ssrt = nan(length(sessionList), 1);
rt = nan(length(sessionList), 1);
poolID = parpool(3);
parfor i = 1 : length(sessionList)
    
    dataInh = ccm_inhibition(subject, sessionList{i}, opt);
    
    if ~isempty(dataInh)
        rt(i) = nanmean(dataInh.allRT);
        ssrt(i) = dataInh.ssrtCollapseIntegrationWeighted;
    end
end
delete(poolID)
ridNan = isnan(rt);
rt(ridNan) = [];
ssrt(ridNan) = [];

[p, s] = polyfit(rt(:), ssrt(:), 1);
% [y, delta] = polyval(p, rt(:), s);
stats = regstats(rt(:), ssrt(:))
fprintf('p-value for regression: %.4f\n', stats.tstat.pval(2))

% R = corrcoef(rt(:), ssrt(:));
Rsqrd = stats.rsquare;

cov(rt(:), ssrt(:));
xVal = min(rt(:)) : max(rt(:));
yVal = p(1) * xVal + p(2);

figure()
hold all;
plot(rt, ssrt, '.',   'color', 'k', 'markersize', 30)
plot(xVal, yVal, 'color', 'k', 'linewidth', 2)


%% Population sdfs for various categories of neurons
subject = 'broca';
sessionList = unique(presacc.sessionID);
rtCorrect = cell(size(sessionList, 1), 4);
rtError = cell(size(sessionList, 1), 4);
chronOpt = ccm_chronometric;
chronOpt.collapseTarg = true;
chronOpt.plotFlag = false;
chronOpt.USE_TWO_COLORS = true;
chronOpt.doStops = false;

for i =  size(sessionList, 1) : size(sessionList, 1)
    i
    iData = ccm_chronometric(subject, sessionList{i}, chronOpt);
    
    rtCorrect(i, 1) = iData.goLeftToTarg(1);
    rtCorrect(i, 2) = iData.goLeftToTarg(2);
    rtCorrect(i, 3) = iData.goRightToTarg(1);
    rtCorrect(i, 4) = iData.goRightToTarg(2);
    
    rtError(i, 1) = iData.goRightToDist(1);
    rtError(i, 2) = iData.goRightToDist(2);
    rtError(i, 3) = iData.goLeftToDist(1);
    rtError(i, 4) = iData.goLeftToDist(2);
    
end

rtCorrect1 = cell2mat(rtCorrect(:, 1));
rtCorrect2 = cell2mat(rtCorrect(:, 2));
rtCorrect3 = cell2mat(rtCorrect(:, 3));
rtCorrect4 = cell2mat(rtCorrect(:, 4));
rtError1 = cell2mat(rtError(:, 1));
rtError2 = cell2mat(rtError(:, 2));
rtError3 = cell2mat(rtError(:, 3));
rtError4 = cell2mat(rtError(:, 4));

accuracy1 = length(rtCorrect1) / (length(rtCorrect1) + length(rtError1));
accuracy2 = length(rtCorrect2) / (length(rtCorrect2) + length(rtError2));
accuracy3 = length(rtCorrect3) / (length(rtCorrect3) + length(rtError3));
accuracy4 = length(rtCorrect4) / (length(rtCorrect4) + length(rtError4));

meanCorrect1 = nanmean(rtCorrect1);
meanCorrect2 = nanmean(rtCorrect2);
meanCorrect3 = nanmean(rtCorrect3);
meanCorrect4 = nanmean(rtCorrect4);
meanError1 = nanmean(rtError1);
meanError2 = nanmean(rtError2);
meanError3 = nanmean(rtError3);
meanError4 = nanmean(rtError4);

varCorrect1 = nanvar(rtCorrect1);
varCorrect2 = nanvar(rtCorrect2);
varCorrect3 = nanvar(rtCorrect3);
varCorrect4 = nanvar(rtCorrect4);
varError1 = nanvar(rtError1);
varError2 = nanvar(rtError2);
varError3 = nanvar(rtError3);
varError4 = nanvar(rtError4);

%%
% Open original popoulation with all the RFs

dataPath = fullfile(projectRoot,'data',projectDate,subject);
o = load(fullfile(dataPath, 'ccm_neuronTypes'));

%Open full list of dingGold neuronTypes
dg = load(fullfile(dataPath, ['ccm_ding_gold_neuronTypes']));

rfList = cell(size(dg.neuronTypes, 1), 1);
% Walk through dingGood neuronTypes, add RF data
for i = 1 : size(o.neuronTypes, 1)
    dgInd = strcmp(dg.neuronTypes.sessionID ,o.neuronTypes.sessionID(i)) & strcmp(dg.neuronTypes.unit ,o.neuronTypes.unit(i));
    rfList(dgInd) = o.neuronTypes.rf(i);
end
neuronTypes = dg.neuronTypes;
neuronTypes.rf = rfList;
save(fullfile(dataPath, 'ccm_ding_gold_neuronTypes'), 'neuronTypes')
%%
% Open original popoulation with all the RFs

dataPath = fullfile(projectRoot,'data',projectDate,subject);
o = load(fullfile(dataPath, 'ccm_neuronTypes'));


%Open full list of dingGold neuronTypes
dg = load(fullfile(dataPath, ['ccm_ding_gold_ddm_Stim_neurons']));

rfList = cell(size(dg.neurons, 1), 1);
% Walk through dingGood neuronTypes, add RF data
for i = 1 : size(o.neuronTypes, 1)
    dgInd = strcmp(dg.neurons.sessionID ,o.neuronTypes.sessionID(i)) & strcmp(dg.neurons.unit ,o.neuronTypes.unit(i));
    rfList(dgInd) = o.neuronTypes.rf(i);
end
neurons = dg.neurons;
neurons.rf = rfList;
save(fullfile(dataPath, 'ccm_ding_gold_ddm_Stim_neurons'), 'neurons')
%%
% Open original popoulation with all the RFs

dataPath = fullfile(projectRoot,'data',projectDate,subject);
o = load(fullfile(dataPath, 'ccm_neuronTypes'));

% load the population of cancel time anlysis
load(fullfile(dataPath, ['ccm_canceled_vs_go_neuronTypes']))
neuronTypes.('hemishphere') = [];

% rfList = cell(size(neuronTypes, 1), 1);
hemisphereList = cell(size(neuronTypes, 1), 1);
% Walk through dingGood neuronTypes, add RF data
for i = 1 : size(o.neuronTypes, 1)
    dgInd = strcmp(neuronTypes.sessionID ,o.neuronTypes.sessionID(i)) & strcmp(neuronTypes.unit ,o.neuronTypes.unit(i));
    % rfList(dgInd) = o.neuronTypes.rf(i);
    hemisphereList(dgInd) = o.neuronTypes.hemisphere(i);
end

% neuronTypes.rf = rfList;
neuronTypes.hemisphere = hemisphereList;
save(fullfile(dataPath, 'ccm_canceled_vs_go_neuronTypes'), 'neuronTypes')



%%
categoryName = 'ding_gold_ddm_';
epoch = 'Stim';

load(fullfile(dataPath, ['ccm_',categoryName, epoch,'_neurons']))



for i = 1 : size(neurons, 1)
    
    % Find corresponding receptive field
    unitInfo = table();
    unitInfo.sessionID  = neuronTypes.sessionID(sessionInd(i));
    unitInfo.unit       = neuronTypes.unit(sessionInd(i));
    unitInfo.hemisphere  = neuronTypes.hemisphere(sessionInd(i));
    
    fprintf('%02d of %d\t%s\t%s\n',i,length(sessionInd), neuronTypes.sessionID{sessionInd(i)},neuronTypes.unit{sessionInd(i)})
    fprintf('Hem: %s\n',neuronTypes.hemisphere{sessionInd(i)})
    
    opt.unitArray = unitInfo.unit;
    opt.hemisphere = neuronTypes.hemisphere{sessionInd(i)};
    opt.printPlot = true;
    
    pdfName = [neuronTypes.sessionID{sessionInd(i)},'_ccm_',neuronTypes.unit{sessionInd(i)},'_neuron.pdf'];
    if exist(fullfile(local_figure_path,subject,pdfName))
        open(fullfile(local_figure_path,subject,pdfName))
    else
        iData = ccm_session_data(subject, neuronTypes.sessionID{sessionInd(i)}, opt);
    end
    
    
    prompt = 'add to list?';
    addToList = input(prompt);
    if addToList
        
        neurons = [neurons; unitInfo];
        save(fullfile(dataPath, ['ccm_',categoryName, epoch, '_neurons']), 'neurons')
        
    end
    clear iData
end
%%
subject = 'broca';
categoryList = {'fix', 'visNoPresacc', 'presaccNoVis', 'visPresacc', 'postsaccNoPresacc'};
for i = 1 : length(categoryList)
    clear Data
    load(fullfile(dataPath, ['ccm_',categoryList{i},'_neurons']))
    
    fprintf('%s\t%d\n',categoryList{i}, size(neurons, 1))
end

%%
categoryList = {'checker'};%, 'fix', 'checker', 'visNoPresacc', 'visPresacc', 'presacc', 'presaccNoVis', 'postsaccNoPresacc'}
opt = ccm_population_neuron_plot;
opt.easyOnly = true;
opt.doStops = false;

for i = 1 : length(categoryList)
    opt.categoryName = categoryList{i};
    
    ccm_population_neuron_plot(subject,projectRoot,projectDate,opt)
end



%%  Compile a list of known duplicate units to
duplicates = table();

duplicates.sessionID = {...
    'bp160n01',...
    'bp229n02-mm',...
    'bp229n02-mm',...
    };

duplicates.Unit = {...
    'spikeUnit17a',...
    'spikeUnit03a',...
    'spikeUnit03b',...
    };

save(fullfile(dataPath, 'ccm_duplicate_neurons'), 'duplicates')

%% comparing cmd and ccm SSRTs
subject = 'joule';
fid=  fopen(fullfile(local_data_path,subject,['cmd_sessions_',subject,'.csv']));
nCol = 5;
formatSpec = '%s';
mHeader = textscan(fid, formatSpec, nCol, 'Delimiter', ',');

mData = textscan(fid, '%s %s %d %d %d', 'Delimiter', ',','TreatAsEmpty',{'NA','na'});
cmdSessionList     = mData{1};


optT = plexon_translate_datafile_mac;
optT.whichData = 'behavior';
for i = 1 : length(cmdSessionList)
    i
    cmdSessionList{i}
    plexon_translate_datafile_mac(subject, cmdSessionList{i}, optT);
end
%%  JOULE CMD
subject = 'joule';
cmdSessionList = {...
    'jp093n01',...
    'jp102n01',...
    'jp103n01',...
    'jp105n01'...
    'jp108n01'...
    'jp109n01'...
    };
opt = cmd_options;
opt.plotFlag = false;
cmdSsrt = nan(length(cmdSessionList), 1);
cmdRt = nan(length(cmdSessionList), 1);
cmdInh = cell(length(cmdSessionList), 1);
cmdInhX = cell(length(cmdSessionList), 1);
for i = 1 : length(cmdSessionList)
    i
    cmdSessionList{i}
    cmd = cmd_inhibition(subject, cmdSessionList{i}, opt);
    cmdSsrt(i) = cmd.ssrtIntegrationWeighted;
    cmdRt(i) = nanmean(cmd.goRT);
    cmdInh{i} = cmd.inhibitionFn;
    cmdInhX{i} = min(cmd.ssdArray) : max(cmd.ssdArray);
end

%%  JOULE CCM
subject = 'joule';
ccmSessionList = {...
    'jp093n02',...
    'jp102n02',...
    'jp103n02',...
    'jp105n02'...
    'jp108n02'...
    'jp109n02'...
    };

opt = ccm_options;
opt.plotFlag = false;
ccmSsrt = nan(length(ccmSessionList), 1);
ccmRt = nan(length(ccmSessionList), 1);
ccmInh = cell(length(ccmSessionList), 1);
ccmInhX = cell(length(ccmSessionList), 1);
for i = 1 : length(ccmSessionList)
    i
    ccmSessionList{i}
    ccm = ccm_inhibition(subject, ccmSessionList{i}, opt);
    ccmSsrt(i) = ccm.ssrtCollapseIntegrationWeighted;
    ccmRt(i) = nanmean(cell2mat(ccm.goTotalRT(:,1)));
    ccmInh{i} = ccm.inhibitionFnGrand;
    ccmInhX{i} = min(cell2mat(ccm.ssd)) : max(cell2mat(ccm.ssd));
end

%% BROCA CMD
subject = 'broca';
fid=  fopen(fullfile(local_data_path,subject,['cmd_sessions_',subject,'.csv']));
nCol = 5;
formatSpec = '%s';
mHeader = textscan(fid, formatSpec, nCol, 'Delimiter', ',');

mData = textscan(fid, '%s %s %d %d %d', 'Delimiter', ',','TreatAsEmpty',{'NA','na'});
cmdSessionList     = mData{1};


optT = plexon_translate_datafile_mac;
optT.whichData = 'behavior';

opt = cmd_options;
opt.plotFlag = false;
cmdSsrt = nan(length(cmdSessionList), 1);
cmdRt = nan(length(cmdSessionList), 1);
cmdInh = cell(length(cmdSessionList), 1);
cmdInhX = cell(length(cmdSessionList), 1);
% for i = 1 : length(cmdSessionList)
%     i
%     cmdSessionList{i}
%     plexon_translate_datafile_mac(subject, cmdSessionList{i}, optT);
% end
for i = 1 : length(cmdSessionList)
    i
    cmdSessionList{i}
    
    cmd = cmd_inhibition(subject, cmdSessionList{i}, opt);
    cmdSsrt(i) = cmd.ssrtIntegrationWeighted;
    cmdRt(i) = nanmean(cmd.goRT);
    cmdInh{i} = cmd.inhibitionFn;
    cmdInhX{i} = min(cmd.ssdArray) : max(cmd.ssdArray);
end
%% BROCA CCM

subject = 'broca';
fid=  fopen(fullfile(local_data_path,subject,['ccm_sessions_',subject,'.csv']));
nCol = 5;
formatSpec = '%s';
mHeader = textscan(fid, formatSpec, nCol, 'Delimiter', ',');

mData = textscan(fid, '%s %s %d %d %d', 'Delimiter', ',','TreatAsEmpty',{'NA','na'});
ccmSessionList     = mData{1};


optT = plexon_translate_datafile_mac;
optT.whichData = 'behavior';

opt = ccm_options;
opt.plotFlag = false;
ccmSsrt = nan(length(ccmSessionList), 1);
ccmRt = nan(length(ccmSessionList), 1);
ccmInh = cell(length(ccmSessionList), 1);
ccmInhX = cell(length(ccmSessionList), 1);
for i = 1 : length(ccmSessionList)
    i
    ccmSessionList{i}
    ccm = ccm_inhibition(subject, ccmSessionList{i}, opt);
    if ~isempty(ccm)
        ccmSsrt(i) = ccm.ssrtCollapseIntegrationWeighted;
        ccmRt(i) = nanmean(cell2mat(ccm.goTotalRT(:,1)));
        ccmInh{i} = ccm.inhibitionFnGrand;
        ccmInhX{i} = min(cell2mat(ccm.ssd)) : max(cell2mat(ccm.ssd));
    end
end

%%
subject = 'JOULE';
figureHandle = 56;
cmdColor = [0 0 0];
ccmColor = [.5 .8 .5];
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(1, 3, figureHandle);
clf
ax(1, 1) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth/1.5 axisHeight]);
hold(ax(1, 1), 'on')
set(ax(1, 1), 'xlim', [.5 2.5], 'ylim', [50 380], 'xtick', [1 2], 'xticklabel', {'CMD', 'CCM'})

plot([2], [nanmean(ccmSsrt)], 'd', 'color', ccmColor, 'markersize', 15)
plot([1], [nanmean(cmdSsrt)], 'd', 'color', cmdColor, 'markersize', 15)
errorbar([2], [nanmean(ccmSsrt)],  [nanstd(ccmSsrt)/sqrt(length(ccmSsrt))], 'linestyle' , 'none', 'color', ccmColor, 'linewidth' , 2)
errorbar([1], [nanmean(cmdSsrt)],  [nanstd(cmdSsrt)/sqrt(length(cmdSsrt))], 'linestyle' , 'none', 'color', cmdColor, 'linewidth' , 2)

plot([2], [nanmean(ccmRt)], 'o', 'color', ccmColor, 'markersize', 15)
plot([1], [nanmean(cmdRt)], 'o', 'color', cmdColor, 'markersize', 15)
errorbar([2], [nanmean(ccmRt)],  [nanstd(ccmRt)/sqrt(length(ccmRt))], 'linestyle' , 'none', 'color', ccmColor, 'linewidth' , 2)
errorbar([1], [nanmean(cmdRt)],  [nanstd(cmdRt)/sqrt(length(cmdRt))], 'linestyle' , 'none', 'color', cmdColor, 'linewidth' , 2)

title(subject)

ax(1, 2) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth*1.7 axisHeight]);
hold(ax(1, 2), 'on')
set(ax(1, 2), 'xlim', [50 650], 'ylim', [0 1])
cellfun(@(x,y) plot(x,y, '-', 'color', ccmColor, 'linewidth', 2), ccmInhX, ccmInh)
cellfun(@(x,y) plot(x,y, '-', 'color', cmdColor, 'linewidth', 2), cmdInhX, cmdInh)
print(figureHandle,fullfile(local_figure_path, [subject, '_ssrt']),'-depsc', '-r300')


%%  Transform tables to variables in save files
subject = 'joule';
% dataDir = '/Volumes/SchallLab/data/Joule';
dataDir = ['/Volumes/SchallLab/data/', subject];
d = dir(dataDir);
for i = 1 : size(d, 1)
    if regexp(d(i).name, 'jp.*mat')
        tic
        
        disp(i)
        disp(d(i).name(1:end-4))
        
        load(fullfile(dataDir, d(i).name))
        trialData = table2struct(trialData, 'ToScalar',true);
        trialData.SessionData = SessionData;
        
        save(fullfile(dataDir, d(i).name), '-struct', 'trialData','-v7.3')
        save(fullfile(local_data_path, subject, d(i).name), '-struct', 'trialData','-v7.3')
        clear trialData SessionData
        
        disp(toc)
    end
end
%%  Transform tables to variables in save files
subject = 'broca';
% dataDir = '/Volumes/SchallLab/data/Joule';
dataDir = ['/Volumes/SchallLab/data/', subject];
d = dir(dataDir);
for i = 1405 : size(d, 1)
    if ~isempty(regexp(d(i).name, 'bp.*mat', 'once')) && isempty(regexp(d(i).name, 'bp.*_legacy.mat', 'once'))
        tic
        disp(i)
        disp(d(i).name(1:end-4))
        
        load(fullfile(dataDir, d(i).name))
        if exist('trialData', 'var')
            switch class(trialData)
                case 'dataset'
                    trialData = dataset2struct(trialData, 'AsScalar',true);
                case 'table'
                    trialData = table2struct(trialData, 'ToScalar',true);
            end
            trialData.SessionData = SessionData;
            
            save(fullfile(dataDir, d(i).name), '-struct', 'trialData','-v7.3')
            save(fullfile(local_data_path, subject, d(i).name), '-struct', 'trialData','-v7.3')
            clear trialData SessionData
        end
        disp(toc)
    end
end

%%  Transform tables to variables in save files
subject = 'xena';
% dataDir = '/Volumes/SchallLab/data/Joule';
dataDir = ['/Volumes/SchallLab/data/', subject, '/Plexon/'];
d = dir(dataDir);
for i = 1 : size(d, 1)
    if ~isempty(regexp(d(i).name, 'xp.*mat', 'once')) && isempty(regexp(d(i).name, 'xp.*_legacy.mat', 'once'))
        tic
        disp(i)
        disp(d(i).name(1:end-4))
        
        load(fullfile(dataDir, d(i).name))
        if exist('trialData', 'var')
            switch class(trialData)
                case 'dataset'
                    trialData = dataset2struct(trialData, 'AsScalar',true);
                case 'table'
                    trialData = table2struct(trialData, 'ToScalar',true);
            end
            trialData.SessionData = SessionData;
            
            save(fullfile(dataDir, d(i).name), '-struct', 'trialData','-v7.3')
            save(fullfile(local_data_path, subject, d(i).name), '-struct', 'trialData','-v7.3')
            clear trialData SessionData
        end
        disp(toc)
    end
end



%%
if isdir('/Volumes/HD-1/Users/paulmiddlebrooks/')
    projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/memory_guided_saccades';
elseif isdir('/Volumes/Macintosh HD/Users/elseyjg/')
    projectRoot = '/Volumes/Macintosh HD/Users/elseyjg/Memory-Guided-Saccade-Project';
else
    disp('You need to add another condition or the file path is wrong.')
end

dataRoot = fullfile(projectRoot, 'data');

%% BROCA MEM and DEL
subject = 'broca';
fid=  fopen(fullfile(dataRoot,subject,['del_sessions_',subject,'.csv']));
nCol = 5;
formatSpec = '%s';
mHeader = textscan(fid, formatSpec, nCol, 'Delimiter', ',');

mData = textscan(fid, '%s %s %d %d %d', 'Delimiter', ',','TreatAsEmpty',{'NA','na'});
cmdSessionList     = mData{1};


optT = plexon_translate_datafile_mac;
optT.hemisphere = 'left';

delSessionList = cmdSessionList(find(strcmp(cmdSessionList, 'bp228n01')) : end);
for i = 1 : length(delSessionList)
    delSessionList{i}
    plexon_translate_datafile_mac(subject, delSessionList{i}, optT);
end



%%
opt = ccm_options;
opt.collapseSignal = true;

opt.unitArray = {'spikeUnit14a', 'spikeUnit14b'}; 
session =  'bp243n02';
for i = 1 : length(sessionList)
    data = ccm_session_data('broca', session, opt);
end


%%
%%matlab

categoryName = 'presacc';
    load(fullfile(dataPath, ['ccm_',categoryName,'_neurons']))


% Establish options to send to ccm_session_data in the for loop below
opt             = ccm_options;
opt.howProcess  = 'print';
opt.plotFlag    = true;
opt.printPlot    = true;
opt.dataType    = 'neuron';
opt.collapseTarg 	= true;
opt.collapseSignal 	= true;
opt.doStops 	= false;
    opt.multiUnit = multiUnit;



for i = 1 : size(neurons, 1)
    opt.unitArray = neurons.unit(i);
    opt.hemisphere = neurons.hemisphere{i};
    
    pdfName = [neurons.sessionID{i},'_ccm_',neurons.unit{i},'_neuron_collapse.pdf'];
    if exist(fullfile(local_figure_path,subject,'sessionCollapseChoice',pdfName))
    else
      iData = ccm_session_data(subject, neuronTypes.sessionID{sessionInd(i)}, opt);
    end
    
    clear iData
end


fileName = fullfile(dataPath, ['ccm_',categoryName,'_neurons', addMulti]);

save(fileName, 'neurons')

%%
opt = ccm_options;
opt.collapseSignal = true;
opt.doStops = false;
sessionList = ...
    {'bp088n02',...
    'bp234n02',...
    'bp247n02'};
sessionList = ...
    {'bp196n02-tp'};
for i = 1 : length(sessionList)
    iSessionID = sessionList{i};
    data = ccm_session_data('broca', iSessionID, opt);
end



%% Population behavioral measures Xena behavior modeled sessions
opt = ccm_options;
opt.plotFlag = false;
opt.printPlot = false;
opt.saveName = 'behavior1';
subject = 'xena';

sessionSet = 'behavior1';
% ccm_chronometric_population(subject, sessionSet, opt);
ccm_psychometric_population(subject, sessionSet, opt);
% ccm_rt_distribution_population(subject, sessionSet, opt);
% data = ccm_inhibition_population(subject, sessionSet, opt);

%% Population behavioral measures: Broca behavior modeled sessions
opt = ccm_options;
opt.plotFlag = false;
opt.printPlot = false;
subject = 'broca';

sessionSet = 'behavior2';
% ccm_chronometric_population(subject, sessionSet, opt);
ccm_psychometric_population(subject, sessionSet, opt);
% ccm_rt_distribution_population(subject, sessionSet, opt);
% data = ccm_inhibition_population(subject, sessionSet, opt);

%% Population behavioral measures:  Broca neural modeled sessions
opt = ccm_options;
opt.plotFlag = false;
opt.printPlot = false;
opt.saveName = 'neural_model';
subject = 'broca';
sessionSet = {...
    'bp242n02';...
    'bp244n02';...
    'bp245n02';...
    'bp246n02';...
    'bp247n02'};
% ccm_chronometric_population(subject, sessionSet, opt);
ccm_psychometric_population(subject, sessionSet, opt);
% ccm_rt_distribution_population(subject, sessionSet, opt);
% data = ccm_inhibition_population(subject, sessionSet, opt);


%% Population behavioral measures: Joule neural modeled sessions
opt = ccm_options;
opt.plotFlag = false;
opt.printPlot = false;
opt.saveName = 'neural_model';
subject = 'joule';

sessionSet = {...
    'jp110n02';...
    'jp114n04';...
    'jp121n02';...
    'jp124n04';...
    'jp125n04'};
% ccm_chronometric_population(subject, sessionSet, opt);
ccm_psychometric_population(subject, sessionSet, opt);
% ccm_rt_distribution_population(subject, sessionSet, opt);
% data = ccm_inhibition_population(subject, sessionSet, opt);


%% Population behavioral measures: Broca all neural sessions
subject = 'broca';
opt = ccm_options;
opt.plotFlag = false;
opt.printPlot = false;
opt.saveName = 'neural_all';

% Figure out which sessions to use
projectDate = '2017-01-11';
projectRoot = '~/perceptualchoice_stop_spikes_population';

addpath(genpath(fullfile(projectRoot,'src/code',projectDate)));
dataPath = fullfile(projectRoot,'data',projectDate,subject);

multiUnit = true;
deleteSessions = true;

if multiUnit
    addMulti = '_multiUnit';
else
    addMulti = [];
end

sessionRemove = ccm_exclude_sessions(subject);

% LOAD ALL RELEVANT DATA AND LISTS

% Load full popuulation of neurons
load(fullfile(dataPath, ['ccm_neuronTypes', addMulti]));
if deleteSessions
    neuronTypes = neuronTypes(~ismember(neuronTypes.sessionID, sessionRemove),:);
end



sessionSet = unique(neuronTypes.sessionID);
% ccm_chronometric_population(subject, sessionSet, opt);
% ccm_psychometric_population(subject, sessionSet, opt);
% ccm_rt_distribution_population(subject, sessionSet, opt);
% data = ccm_inhibition_population(subject, sessionSet, opt);

%% Population behavioral measures: Joule all neural sessions
subject = 'joule';
opt = ccm_options;
opt.plotFlag = false;
opt.printPlot = false;
opt.saveName = 'neural_all';

% Figure out which sessions to use
projectDate = '2017-01-11';
projectRoot = '~/perceptualchoice_stop_spikes_population';

addpath(genpath(fullfile(projectRoot,'src/code',projectDate)));
dataPath = fullfile(projectRoot,'data',projectDate,subject);

multiUnit = true;
deleteSessions = true;

if multiUnit
    addMulti = '_multiUnit';
else
    addMulti = [];
end

sessionRemove = ccm_exclude_sessions(subject);

% LOAD ALL RELEVANT DATA AND LISTS

% Load full popuulation of neurons
load(fullfile(dataPath, ['ccm_neuronTypes', addMulti]));
if deleteSessions
    neuronTypes = neuronTypes(~ismember(neuronTypes.sessionID, sessionRemove),:);
end



sessionSet = unique(neuronTypes.sessionID);
% ccm_chronometric_population(subject, sessionSet, opt);
% ccm_psychometric_population(subject, sessionSet, opt);
% ccm_rt_distribution_population(subject, sessionSet, opt);
% data = ccm_inhibition_population(subject, sessionSet, opt);


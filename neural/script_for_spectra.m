%%
% Load the data
        addVar = 'lfpData';
variables = [ccm_min_vars, addVar, 'trialOnset', 'trialDuration'];

[td, S, E] = load_data('broca','bp131n02', variables);

pSignalArray = E.pSignalArray;
nTrial = size(td, 1);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 KDATA SCRIPTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleHz = 1000;
samplePerMS = 1;

lfpData = cellfun(@(x) x(1:end-1), td.lfpData, 'uni', false);
lfpData = cell2mat(lfpData);

alignTime = td.trialOnset + td.targOn;


% Select trials
opt                 = ccm_trial_selection;
opt.outcome         = {'goCorrect', 'targetHoldAbort'};
opt.rightCheckerPct = 100 * pSignalArray(pSignalArray < .5);
goCorrTrialL        = ccm_trial_selection(td, opt);


opt.rightCheckerPct = 100 * pSignalArray(pSignalArray > .5);
goCorrTrialR        = ccm_trial_selection(td, opt);



Fs      = 1000; % sampling frequency
T       = 1/Fs; % sample time
L       = 250;
sigT    = 500;
% L_s   = Step_Offset_time_of_interest_s-Step_Onset_time_of_interest+1; % length of signal
% L_l   = Step_Offset_time_of_interest_l-Step_Onset_time_of_interest+1; % length of signal

t_sig   = (0 : sigT - 1) * T; % Time Vector
% t_s   = (0:L_s-1)*T; % Time Vector
% t_l   = (0:L_l-1)*T; % Time Vector



NFFT    = 2^nextpow2(L); %Next power of 2 from length of
f       = Fs/2*linspace(0,1,NFFT/2+1);

hzMax           = 120; %
epochDuration   = 2000; % how many ms to analyze after event of inerest?
window          = 100; % ms
windowStep      = 10; % step in windowStep ms increments
nStep           = 1 + (epochDuration - window) / windowStep; % How many steps are possible given the epochDuration, window width, and windowStep?
nOverlap        = window - windowStep; %Number of overlapping window frames

% Initialize spec matrix
spec = nan(nTrial, hzMax, nStep);

for i = 1 : length(goCorrTrialL)
    iTrial = goCorrTrialL(i);
    
    clear S F T p;
    
    [S,F,T,P]= spectrogram(td.lfpData{iTrial}(td.targOn(iTrial) : td.targOn(iTrial) + epochDuration-1), window, nOverlap, NFFT, Fs);
    spec(i,:,:) = P(1:hzMax, :);
    
    plotTrials = 1;
    if plotTrials
    % Pllot
     figure(4); hold on;
     cla
   surf(T,F(1:40),10*log10(P(1:40, :)),'edgecolor','none')
    x = [td.checkerOn(iTrial) - td.targOn(iTrial) td.checkerOn(iTrial) - td.targOn(iTrial)] / Fs;
    rt = [td.responseOnset(iTrial) - td.targOn(iTrial) td.responseOnset(iTrial) - td.targOn(iTrial)] / Fs;
    y = [F(1) F(40)];
   	z = [1 1];
    plot3(x,y,z, 'k', 'linewidth', 3);
    plot3(rt,y,z, '--k', 'linewidth', 3);
    axis tight, view(0,90)
    xlabel 'Time (s)', ylabel 'Frequency (Hz)', title 'TargOn'
    pause
    end
    
end

return
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 FIEDLTRIP SCRIPTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The trial definition "trl" is an Nx3 matrix, N is the number of trials.
% The first column contains the sample-indices of the begin of each trial
% relative to the begin of the raw data, the second column contains the
% sample-indices of the end of each trial, and the third column contains
% the offset of the trigger with respect to the trial. An offset of 0
% means that the first sample of the trial corresponds to the trigger. A
% positive offset indicates that the first sample is later than the trigger,
% a negative offset indicates that the trial begins before the trigger.

nTrial = size(td, 1);
trl = nan(nTrial, 3);

trigger = -td.targOn;

trl(:,1) = td.trialOnset;
trl(:,2) = td.trialDuration + td.trialOnset;
trl(:,3) = trigger;


lfpData = cellfun(@(x) x(1:end-1), td.lfpData, 'uni', false);
lfpData = cell2mat(lfpData);




%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 CHRONUS SCRIPTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


duration = 2000;
datalfp = td.lfpData(:,1);
% datalfp = align_signals(datalfp, td.fixWindowEntered, 101:600);
datalfp1 = align_signals(datalfp, td.targOn, -duration:-1);
datalfp2 = align_signals(datalfp, td.targOn, 1:duration);

time = 1:duration;
test = sin(10*2*pi/1000*time);

datalfp3 = test;




datalfp1 = datalfp1';
datalfp2 = datalfp2';


params=struct('tapers',[5 9],'pad',2,'Fs',1000,'fpass',[0 100],'err',0,'trialave',0);
params.err=[1 0.05];


%% Compute spectrum for first trial
% [S,f,Serr]=mtspectrumc(detrend(datalfp(:,1)),params);
[S1,f1,Serr1]=mtspectrumc(datalfp1(:,1),params);
[S2,f2,Serr2]=mtspectrumc(datalfp2(:,1),params);
[S3,f3,Serr3]=mtspectrumc(datalfp3,params);

% plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:)));
% plot(f,S.^2, f,Serr(1,:).^2, f,Serr(2,:).^2);
figure(1); hold on;
% plot(f1,S1, f1,Serr1(1,:), f1,Serr1(2,:));
% plot(f2,S2, f2,Serr2(1,:), f2,Serr2(2,:));
plot(f1,S1 ,'r');
plot(f2,S2, 'g');
% plot(f3,S3, 'k');

figure(3); hold on;
plot(datalfp1(:,1), 'r');
plot(datalfp2(:,1), 'g');
% plot(datalfp3, 'k');




%% Get the left and right trials
opt = ccm_trial_selection;
opt.outcome = {'goCorrect', 'targetHoldAbort'};
opt.rightCheckerPct = 100 * pSignalArray(pSignalArray < .5);
goCorrTrialL = ccm_trial_selection(td, opt);


opt.rightCheckerPct = 100 * pSignalArray(pSignalArray > .5);
goCorrTrialR = ccm_trial_selection(td, opt);


%% Set up plot
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(2, 2, 50);
clf
ax(1,1) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
cla
hold(ax(1,1), 'on')
ax(1,2) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
cla
hold(ax(1,2), 'on')

ax(2,1) = axes('units', 'centimeters', 'position', [xAxesPosition(2,1) yAxesPosition(2,1) axisWidth axisHeight]);
cla
hold(ax(2,1), 'on')

ax(2,2) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2) yAxesPosition(2, 2) axisWidth axisHeight]);
cla
hold(ax(2,2), 'on')


%% Plot the raw LFP
plot(ax(1,1), mean(datalfp1(:,goCorrTrialL), 2), 'b', 'lineWidth', 2)
plot(ax(1,1), mean(datalfp1(:,goCorrTrialR), 2), 'r', 'lineWidth', 2)

plot(ax(1,2), mean(datalfp2(:,goCorrTrialL), 2), 'b', 'lineWidth', 2)
plot(ax(1,2), mean(datalfp2(:,goCorrTrialR), 2), 'r', 'lineWidth', 2)


%% Compute average spectrum for correct trials to left:
params.trialave = 1;
[SL1,fL1,SerrL1]=mtspectrumc(datalfp1(:, goCorrTrialL), params);
[SR1,fR1,SerrR1]=mtspectrumc(datalfp1(:, goCorrTrialR), params);

[SL2,fL2,SerrL2]=mtspectrumc(datalfp2(:, goCorrTrialL), params);
[SR2,fR2,SerrR2]=mtspectrumc(datalfp2(:, goCorrTrialR), params);

plot(ax(2,1), fL1,SL1, 'b', 'linewidth', 2)%, fL,SerrL(1,:).^2, fL,SerrL(2,:).^2);
plot(ax(2,1), fR1,SR1, 'r', 'linewidth', 2)%, fL,SerrL(1,:).^2, fL,SerrL(2,:).^2);

plot(ax(2,2), fL2,SL2, 'b', 'linewidth', 2)%, fL,SerrL(1,:).^2, fL,SerrL(2,:).^2);
plot(ax(2,2), fR2,SR2, 'r', 'linewidth', 2)%, fL,SerrL(1,:).^2, fL,SerrL(2,:).^2);


disp('done')

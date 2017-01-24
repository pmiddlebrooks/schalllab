
%% current translation code

sessionID = 'bp067n01';
fileDir = 'local_data/';
load([fileDir, sessionID])

nTrial = size(trialData.trialOutcome, 1);

figure(912)
for iTrial = 4 : nTrial
    cla
    hold on
    fprintf('%s \nSignal: %d \tTarg: %d \tChecker: %d \tResponse: %d\n', trialData.trialOutcome{iTrial}, round(trialData.targ1CheckerProp(iTrial)*100), trialData.targOn(iTrial), trialData.checkerOn(iTrial), trialData.responseOnset(iTrial))
    trialData.saccBegin{iTrial}
    
    plot(trialData.eyeX{iTrial}, trialData.eyeY{iTrial}, 'k')
    if ~isnan(trialData.saccToTargIndex(iTrial));
        fprintf('Response: %d\n', trialData.saccBegin{iTrial}(trialData.saccToTargIndex(iTrial)))
        saccSpan = trialData.saccBegin{iTrial}(trialData.saccToTargIndex(iTrial)) : trialData.saccBegin{iTrial}(trialData.saccToTargIndex(iTrial)) + trialData.saccDuration{iTrial}(trialData.saccToTargIndex(iTrial));
        plot(trialData.eyeX{iTrial}(saccSpan), trialData.eyeY{iTrial}(saccSpan), 'b', 'lineWidth', 3)
    end
    nSacc = length(trialData.saccBegin{iTrial});
    for jSacc = 1 : nSacc
        fprintf('Saccade: %d\n', trialData.saccBegin{iTrial}(jSacc))
        rectangle('Position', [-trialData.fixWindow(iTrial)/2, -trialData.fixWindow(iTrial)/2, trialData.fixWindow(iTrial), trialData.fixWindow(iTrial)])
        targX = trialData.targAmp(iTrial) * cosd(trialData.targAngle(iTrial));
        targY = trialData.targAmp(iTrial) * sind(trialData.targAngle(iTrial));
        rectangle('Position', [targX - trialData.targWindow(iTrial)/2, targY - trialData.targWindow(iTrial)/2, trialData.targWindow(iTrial), trialData.targWindow(iTrial)], 'edgecolor', 'b')
        rectangle('Position', [-targX - trialData.targWindow(iTrial)/2, -targY - trialData.targWindow(iTrial)/2, trialData.targWindow(iTrial), trialData.targWindow(iTrial)])
        saccSpan = trialData.saccBegin{iTrial}(jSacc) : trialData.saccBegin{iTrial}(jSacc) + trialData.saccDuration{iTrial}(jSacc);
        plot(trialData.eyeX{iTrial}(saccSpan), trialData.eyeY{iTrial}(saccSpan), 'k', 'lineWidth', 2)
        plot(trialData.eyeX{iTrial}(trialData.saccBegin{iTrial}(jSacc)), trialData.eyeY{iTrial}(trialData.saccBegin{iTrial}(jSacc)), '.g', 'markersize', 25)
        plot(trialData.eyeX{iTrial}(trialData.saccBegin{iTrial}(jSacc) + trialData.saccDuration{iTrial}(jSacc)), trialData.eyeY{iTrial}(trialData.saccBegin{iTrial}(jSacc) + trialData.saccDuration{iTrial}(jSacc)), '.r', 'markersize', 25)
        pause
    end
end
%% Old translation cod
sessionID = 'bp067n01';
fileDir = 'local_data/';
load([fileDir, sessionID])

nTrial = size(trialData.trialOutcome, 1);

figure(912)
for iTrial = 1 : nTrial
    cla
    hold on
    fprintf('%s \nSignal: %d \tTarg: %d \tChecker: %d \tResponse: %d\n', trialData.trialOutcome{iTrial}, round(trialData.target1CheckerProportion{iTrial}*100), trialData.targetOnset{iTrial}, trialData.choiceStimulusOnset{iTrial}, trialData.responseOnset{iTrial})
    trialData.saccadeOnset{iTrial}
    
    plot(trialData.eyePositionX{iTrial}, trialData.eyePositionY{iTrial}, 'k')
    if ~isnan(trialData.saccadeToTargetIndex{iTrial});
        fprintf('Response: %d\n', trialData.saccadeOnset{iTrial}(trialData.saccadeToTargetIndex{iTrial}))
        saccSpan = trialData.saccadeOnset{iTrial}(trialData.saccadeToTargetIndex{iTrial}) : trialData.saccadeOnset{iTrial}(trialData.saccadeToTargetIndex{iTrial}) + trialData.saccadeDuration{iTrial}(trialData.saccadeToTargetIndex{iTrial});
        plot(trialData.eyePositionX{iTrial}(saccSpan), trialData.eyePositionY{iTrial}(saccSpan), 'b', 'lineWidth', 3)
    end
    nSacc = length(trialData.saccadeOnset{iTrial});
    for jSacc = 1 : nSacc
        fprintf('Saccade: %d\n', trialData.saccadeOnset{iTrial}(jSacc))
        rectangle('Position', [-trialData.fixationWindow{iTrial}/2, -trialData.fixationWindow{iTrial}/2, trialData.fixationWindow{iTrial}, trialData.fixationWindow{iTrial}])
        saccSpan = trialData.saccadeOnset{iTrial}(jSacc) : trialData.saccadeOnset{iTrial}(jSacc) + trialData.saccadeDuration{iTrial}(jSacc);
        plot(trialData.eyePositionX{iTrial}(saccSpan), trialData.eyePositionY{iTrial}(saccSpan), 'k', 'lineWidth', 2)
        plot(trialData.eyePositionX{iTrial}(trialData.saccadeOnset{iTrial}(jSacc)), trialData.eyePositionY{iTrial}(trialData.saccadeOnset{iTrial}(jSacc)), '.g', 'markersize', 25)
        plot(trialData.eyePositionX{iTrial}(trialData.saccadeOnset{iTrial}(jSacc) + trialData.saccadeDuration{iTrial}(jSacc)), trialData.eyePositionY{iTrial}(trialData.saccadeOnset{iTrial}(jSacc) + trialData.saccadeDuration{iTrial}(jSacc)), '.r', 'markersize', 25)
        pause
    end
end

%%
%%
sessionID = 'xp054n02';
fileDir = 'local_data/';
load([fileDir, sessionID])

trialData = cell_to_mat(trialData);

figure(238)
clf
hold all
plot(trialData.fixOn, '.-k', 'markersize', 10)
plot(trialData.targOn, '.-g', 'markersize', 10)
plot(trialData.stopSignalOn, '.-r', 'markersize', 10)
plot(trialData.checkerOn, '.-b', 'markersize', 10)




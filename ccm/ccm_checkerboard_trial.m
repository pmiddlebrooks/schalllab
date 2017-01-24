%%
subjectID = 'Broca';
% sessionID = 'bp087n02';
% sessionID = 'bp088n02';
% sessionID = 'bp089n02';
sessionID = 'bp090n02';
sessionID = 'bp091n02';
sessionID = 'bp092n02';
sessionID = 'bp093n02';
sessionID = 'bp107n01';
% sessionID = 'bp095n02';
% sessionID = '0724saccade';
% subjectID = 'xb';



% Load the data
[trialData, SessionData, targ1PropArray, ssdArray] = load_data(subjectID, sessionID);

%%
nCheckerRow = 10;
nCheckerCol = 10;
        if iscell(trialData.checkerArray)
            % reshape cell contents if necessary
            if sum(cellfun(@(x) size(x, 1)>1, trialData.checkerArray))
                trialData.checkerArray = cellfun(@(x) x', trialData.checkerArray, 'uniformoutput', false);
            end
            trialData.checkerArray = cell2mat(trialData.checkerArray);
        end

figure(77)
for i = 1 : size(trialData.trialOutcome, 1)
        colormap(flipud(cool(64)));
    if isnan(trialData.checkerOn(i))
        continue
    end
        c = reshape(trialData.checkerArray(i,:), nCheckerRow, nCheckerCol)';
        
        horiz = sum(c, 1)
        vert = sum(c, 2)
        c(c == 0) = size(colormap, 1);
image(c)
pause
end

%%
nCheckerRow = 10;
nCheckerCol = 10;
        if iscell(trialData.checkerArray)
            % reshape cell contents if necessary
            if sum(cellfun(@(x) size(x, 1)>1, trialData.checkerArray))
                trialData.checkerArray = cellfun(@(x) x', trialData.checkerArray, 'uniformoutput', false);
            end
            trialData.checkerArray = cell2mat(trialData.checkerArray);
        end

horiz = zeros(1, nCheckerCol);
vert = zeros(nCheckerRow, 1);
nTrial = size(trialData.trialOutcome(~isnan(trialData.checkerOn),:), 1);
for i = 1 : size(trialData.trialOutcome, 1)
        colormap(flipud(cool(64)));
    if isnan(trialData.checkerOn(i))
        continue
    end
        c = reshape(trialData.checkerArray(i,:), nCheckerRow, nCheckerCol)';
        
        horiz = horiz + sum(c, 1);
        vert = vert + sum(c, 2);
end
horiz ./ nTrial
vert ./ nTrial

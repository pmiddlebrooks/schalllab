function unitInfo = ccm_classify_neuron_ding_gold(Data)
%
%
% function neuronTable = ccm_classify_neuron_ding_gold(Data);
%
% Classify neurons according to the criteria of Ding & Gold 2012 Cereb Ctx.
% They query "choice" and "coherence" dependence within 4 defined epochs:
% Stim, Sacc, Post, Reward.


%%

% Define the four epochs
epoch = Data.epoch; % {'Stim', 'Sacc', 'Post', 'Reward'};
switch epoch
    case 'Stim'
        epochWindow = 81:150;
        epochAlign = 'checkerOn';
    case 'Sacc'
        epochWindow = -49 : 0;
        epochAlign = 'responseOnset';
    case 'Post'
        epochWindow = 41:140;
        epochAlign = 'responseOnset';
    case 'Reward'
        epochWindow = 51:250;
        epochAlign = 'rewardOn';
end
rasterWindow = -499: 900;



% These will determine the trial-to-trial epochs used for analyses:
preSaccadeBuffer = 10;
minEpochDuration = 10; % Only include trials for which the determined epoch is this long

% Initialize the table and add session information
% unitInfo = table();
unitInfo = cell(1, 12);







% Initialize table row
%     iUnitInfo = table();
%     iUnitInfo.sessionID  = {Data.sessionID};
%     iUnitInfo.unit       = {Data.name};
%     iUnitInfo.hemisphere       = {Data.hemisphere};
%     iUnitInfo.rf       = {Data.rf};
unitInfo(1)  = {Data.sessionID};
unitInfo(2)      = {Data.name};
unitInfo(3)       = {Data.hemisphere};
unitInfo(4)       = {Data.rf};

% Loop through color coherence (singal strength) values and build
% vectors of relevant things to send to ding_gold_ddm_like.m
iSpikeRate       = [];
iColorCoherence  = [];
iAlignedRasters  = {};
epochAlignInd = -rasterWindow(1);

for j = 1 : length(Data.pSignalArray)
    
    nTrial = length(Data.checkerOn.colorCoh(j).goTarg.rt);
    jAlign = Data.(epochAlign).colorCoh(j).goTarg.alignTime;
    jRaster = num2cell(Data.(epochAlign).colorCoh(j).goTarg.raster(:, jAlign + rasterWindow), 2);
    
    switch epoch
        % If we're dealing with the Stim epoch, want to define the end of
        % the epoch with respect to median RT for a given side and
        % trial-by-trial RTs for RTs less than the median RT
        case 'Stim'
            if Data.pSignalArray(j) < .5
                medianRT = round(nanmean(Data.checkerOn.colorCoh(1).goTarg.rt));
            else
                medianRT = round(nanmean(Data.checkerOn.colorCoh(end).goTarg.rt));
            end
            % Initialize epochEnd as xx ms (preSaccadeBuffer) before
            % median RT
            epochEnd = medianRT * ones(nTrial, 1);
            epochBegin = ceil(.5 * epochEnd);
            
            
            % Replace epoch-cutoffs for trials with rts shorter than the median RT
            earlySaccTrial = Data.checkerOn.colorCoh(j).goTarg.rt < epochEnd + preSaccadeBuffer;
            epochEnd(earlySaccTrial) = Data.checkerOn.colorCoh(j).goTarg.rt(earlySaccTrial) - preSaccadeBuffer;
            
            
            % Adjust for the alignment index
            epochEnd = epochEnd + epochAlignInd;
            epochBegin = epochBegin + epochAlignInd;
            
            % If there are trials with negative epochs because of the
            negativeEpochTrial = epochEnd < epochBegin + minEpochDuration;
            epochBegin(negativeEpochTrial) = epochEnd(negativeEpochTrial) - minEpochDuration;
            
            
            % Make sure epoch doesn't extend beyond length of rasters
            epochEnd(epochEnd > length(jRaster{1})) = length(jRaster{1});
            
            nSpike = cellfun(@(x,y,z) sum(x(y:z)), jRaster, num2cell(epochBegin), num2cell(epochEnd));
            jSpikeRate = nSpike .* 1000 ./ (epochEnd - epochBegin);
            
        otherwise
            epochBegin = epochAlignInd + epochWindow(1) * ones(nTrial, 1);
            epochEnd = epochAlignInd + epochWindow(end) * ones(nTrial, 1);
            nSpike = cellfun(@(x,y,z) sum(x(y:z)), jRaster, num2cell(epochBegin), num2cell(epochEnd));
            jSpikeRate = nSpike .* 1000 ./ (epochWindow(end) - epochWindow(1));
    end
    
    
    iColorCoherence = [iColorCoherence; Data.pSignalArray(j) * ones(nTrial, 1)];
    iAlignedRasters = [iAlignedRasters; jRaster];
    iSpikeRate = [iSpikeRate; jSpikeRate];
end




% Call ding_gold_ddm_like for this epoch for this set of color
% coherence values
iLeftTrial = iColorCoherence < 0.5;
iRightTrial = iColorCoherence > 0.5;


Unit.spikeRate      = iSpikeRate;
Unit.leftTrial      = iLeftTrial;
Unit.rightTrial     = iRightTrial;
Unit.signalP        = iColorCoherence;
Unit.alignedRasters = iAlignedRasters;
Unit.alignIndex     = epochAlignInd;
Unit.epochOffset    = epochWindow(1);


[ddmData] = ding_gold_ddm_like(Unit);


unitInfo{5}         = epoch;
unitInfo{6}        = ddmData.choiceDependent;
unitInfo{7}     = ddmData.coherenceDependent;
unitInfo{8}           = ddmData.ddmLike; % true if this activity counts as ddm-like, false if not
unitInfo{9}       = ddmData.tChoice; % if ddm, this is time of coherence difference onset a la ding_gold_2012
unitInfo{10}      = ddmData.leftIsIn; % true of left direction counts as into RF, false if not
unitInfo{11}       = ddmData.coeffIn;  % regression coefficients into RF
unitInfo{12}      = ddmData.coeffOut; % regression coefficients out of RF

%     unitInfo = [unitInfo; iUnitInfo];



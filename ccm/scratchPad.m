%% how many trials comprise each condition?
jTarg = 1;
kDataIndex = 1;
mEpochName = 'checkerOn';
iSSD = 3;

sumGoTarg = 0;
sumGoDist = 0;
for iPropIndex = 1 : 4
    
    nGoTarg = size(data(kDataIndex, jTarg).signalStrength(iPropIndex).goTarg.(mEpochName).raster, 1);
    nGoDist = size(data(kDataIndex, jTarg).signalStrength(iPropIndex).goDist.(mEpochName).raster, 1);

    sprintf('*************** Unit %i, Condition %i *************',kDataIndex,iPropIndex)
    sprintf('goTarg condiiton %i, %i trials', iPropIndex, size(data(kDataIndex, jTarg).signalStrength(iPropIndex).goTarg.(mEpochName).raster, 1))
    sprintf('goDist condiiton %i, %i trials', iPropIndex, size(data(kDataIndex, jTarg).signalStrength(iPropIndex).goDist.(mEpochName).raster, 1))
    
    
    sprintf('stopTarg condiiton %i, %i trials', iPropIndex, size(data(kDataIndex, jTarg).signalStrength(iPropIndex).stopTarg.ssd(iSSD).(mEpochName).raster, 1))
    sprintf('stopDist condiiton %i, %i trials', iPropIndex, size(data(kDataIndex, jTarg).signalStrength(iPropIndex).stopDist.ssd(iSSD).(mEpochName).raster, 1))
    sprintf('stopStop condiiton %i, %i trials', iPropIndex, size(data(kDataIndex, jTarg).signalStrength(iPropIndex).stopStop.ssd(iSSD).(mEpochName).raster, 1))

    sumGoTarg = sumGoTarg + nGoTarg;
    sumGoDist = sumGoDist + nGoDist;
    
end
sumGoTarg
sumGoDist
%%


sefContra = load('rexroutines/data/population/epochAnalysis/sef_saccade_contra_dataset_Decision Visual Extended.mat')

shuffContra = sefContra.data(strncmp(sefContra.data.sessionID, 'sp', 1),:)


chContra = cellfun(@(x) nanmean(x, 1), shuffContra.sdfCH, 'uni', false);
clContra = cellfun(@(x) nanmean(x, 1), shuffContra.sdfCL, 'uni', false);
ihContra = cellfun(@(x) nanmean(x, 1), shuffContra.sdfIH, 'uni', false);
ilContra = cellfun(@(x) nanmean(x, 1), shuffContra.sdfIL, 'uni', false);


Y = nan(size(shuffContra, 1),size(epochWindow),8);

Y(:,:,1) = cell2mat(chContra);
Y(:,:,2) = cell2mat(clContra);
Y(:,:,3) = cell2mat(ihContra);
Y(:,:,4) = cell2mat(ilContra);

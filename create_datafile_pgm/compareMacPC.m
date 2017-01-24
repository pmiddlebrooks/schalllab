%%
[td, S] = load_data('broca','bp093n01');
[nd, S] = load_data('broca','bp093n01_test');
figure(4)
clf
hold all

i = 1;
plot(td.eyeX{i}, td.eyeY{i},'b')
plot(nd.eyeX{i}, nd.eyeY{i},'r')

xGain = td.eyeX{i} ./ nd.eyeX{i};
yGain = td.eyeY{i} ./ nd.eyeY{i};

%%
[td, S] = load_data('broca','bp093n01');
[nd, S] = load_data('broca','bp093n01_test');
figure(5)
clf
hold all

i = 1;

ras1 = spike_to_raster(td.spikeData{i,1});
ras2 = spike_to_raster(nd.spikeData{i,1});
sdf1 = spike_density_function(ras1);
sdf2 = spike_density_function(ras2);
plot(sdf1,'b')
plot(sdf2,'r')

%%
eyeChannelX = 61;
eyeChannelY = 64;
plx = readPLXFileC('local_data/CalibrationEyeGain.plx','fullread',     'all');

eyeX = plx.ContinuousChannels(eyeChannelX).Values;
eyeX    = eyeX ./ plx.ContinuousChannels(eyeChannelX).ADGain;
eyeX    = eyeX ./ plx.ContinuousChannels(eyeChannelX).PreAmpGain;
eyeY = plx.ContinuousChannels(eyeChannelY).Values;
eyeY    = eyeY ./ plx.ContinuousChannels(eyeChannelY).ADGain;
eyeY    = eyeY ./ plx.ContinuousChannels(eyeChannelY).PreAmpGain;


figure(5)
clf
plot(eyeX, eyeY)
hold all;
for i = 1 : 300 : length(eyeX)
    
    plot(eyeX(i), eyeY(i), 'og', 'markersize', 30, 'markerfacecolor', 'g')
    sprintf('eyeX: %.3f\teyeY: %.3f\n',eyeX(i),eyeY(i))
    pause
end

%%
% [td, S] = load_data('broca','bp093n01');
% [nd, S] = load_data('broca','bp093n01_test');
[td, S] = load_data('broca','bp093n02');
[nd, S] = load_data('broca','bp093n02_test');

vScale = 2.4414;

figure(5)
clf
plot(eyeX, eyeY)
for i = 1 : size(td, 1)
    clf
    plot(td.eyeX{i}, td.eyeY{i}, 'k')
    hold on;
    plot(nd.eyeX{i}, nd.eyeY{i}, 'g')
    plot(nd.eyeX{i} * vScale, nd.eyeY{i} * vScale, 'r')
    sprintf('eyeX: %.3f\teyeY: %.3f\n',eyeX(i),eyeY(i))
    pause
end
%%
% [td, S] = load_data('broca','bp093n01');
% [nd, S] = load_data('broca','bp093n01_test');
[td, S] = load_data('broca','bp093n02');
[nd, S] = load_data('broca','bp093n02_test');
%%

figure(5)
clf
for i = 1 : size(td, 1)
    [tdAR, alignmentIndex] = spike_to_raster(td.spikeData(i, 1));
    [ndAR, alignmentIndex] = spike_to_raster(nd.spikeData(i, 1));
    tdSDF = spike_density_function(tdAR);
    ndSDF = spike_density_function(ndAR);
    
    
    clf
    plot(tdSDF, 'k')
    hold on;
    plot(ndSDF, 'b')
    sprintf('Old spikes: %d\tNew spikes: %d\n',length(td.spikeData{i, 1}),length(nd.spikeData{i, 1}))
    pause
end
%%

figure(5)
clf
nOldSpike = nan(size(td, 1),1);
nNewSpike = nan(size(td, 1),1);
for i = 1 : size(td, 1)
    nOldSpike(i) = length(td.spikeData{i, 1});
    nNewSpike(i) = length(nd.spikeData{i, 1});
    
%     clf
%     plot(tdSDF, 'k')
%     hold on;
%     plot(ndSDF, 'b')
%     sprintf('Old spikes: %d\tNew spikes: %d\n',length(td.spikeData{i, 1}),length(nd.spikeData{i, 1}))
%     pause
end
hist(nOldSpike, 'k')
hold on
hist(nNewSpike, 'b')
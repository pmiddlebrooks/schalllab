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

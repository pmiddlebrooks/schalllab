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


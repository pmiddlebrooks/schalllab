function [adValues ADname] = plexon_ADchannels_translation(path,...
    file,...
    channel,...
    eyeXGain,...
    eyeYGain,...
    eyeXOffset,...
    eyeYOffset,...
    dateNumber)

% see: plx2mat_ADchannels.m
if nargin < 4
    eyeXGain      = 1;
    eyeYGain      = 1;
    eyeXOffset    = 0;
    eyeYOffset    = 0;
    dateNumber = 0;
end

% this is rig and date specific
% conditional statements may be added
if dateNumber < 734597
    xChannel = 49;
    yChannel = 50;
elseif dateNumber >= 734597 % 04-Apr-2011 (moved eye channels to top)
    xChannel = 63;
    yChannel = 64;
end


% name the channel 
if channel < 10
    ADname = ['AD0',num2str(channel)];
else
    ADname = ['AD',num2str(channel)];
end

if channel == xChannel
    ADname = 'eyeX';
elseif channel == yChannel
    ADname = 'eyeY';
end


% get the data for the channel
[adfreq, n, startTimes, frag_ns, adValues] = plx_ad_v([path file],channel-1);

% get starting time for each fragment of analog data. (if data collection
% was interrupted there will be multiple fragments.)

% actual time stamps
startTimes = round(startTimes * 1000);
endTimes   = frag_ns + startTimes - 1;

startFrags = ones(1, length(startTimes));  % Initialize to one: first value is always one
for iFrag = 2 : length(startTimes)
    startFrags(iFrag) = startFrags(iFrag - 1) + frag_ns(iFrag - 1);
end

if ~isempty(startFrags)
    endFrags(1 : length(startFrags ) - 1) = startFrags(2 : end) - 1;
    endFrags(end + 1) = length(adValues);
else
    endFrags = length(adValues);
end

% now fill in zeros for all the time when we weren't actually collecting
% data.
timeCorrectedData = zeros(endTimes(end), 1);
for iStartTime = 1 : length(startTimes)
    timeCorrectedData(startTimes(iStartTime)+1 : endTimes(iStartTime)+1) = adValues(startFrags(iStartTime) : endFrags(iStartTime));
end
adValues = timeCorrectedData;
clear timeCorrectedData


% if it is an eye channel, do some post processing
if strcmp(ADname,'eyeX') ||...
        strcmp(ADname,'eyeY')
    % convert the channel from voltage to degrees
    if channel == xChannel
        adValues = (adValues * eyeXGain) - eyeXOffset; %this may need to be done to take vector into account
    elseif channel == yChannel
        adValues = (adValues * eyeYGain) - eyeYOffset; %this may need to be done to take vector into account
        adValues = -adValues;
    end    
    
    % gaussian polynomial for convolving eye traces
    polyWidth = 12;
    polyn = gaussmf(polyWidth/-2 : polyWidth/2, [polyWidth/4,0]);
    polyn = polyn/sum(polyn); % now it sums to 1
    
    %smooth out analog line noise
    adValues = conv_2009a(adValues, polyn, 'same'); 
end   





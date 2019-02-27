function cMap = ccm_colormap(signalArray, grayscale)

% If percentages were input instead of proportions, change it to
% proportions
if sum(signalArray < 1) == 0
    signalArray = signalArray .* 100;
end
if nargin == 1
    grayscale = false;
end
nSignalLevel = length(signalArray);



pureCyan = [0 1 1];
pureMagenta = [1 0 1];

cyanMax = 0.85;
magMax = 1;
minAll = .1;
lightestGray = 1;
switch nSignalLevel
    case 2
        if grayscale
            cMap = [0 0 0; .5 .5 .5];
        else
            cMap = [0 .67 .67; .67 0 .67];
        end
    otherwise
        fiftyP = ismember(signalArray, .5);
        nFifty = sum(fiftyP);
        
        leftArray = signalArray(signalArray < .5);
        rightArray = signalArray(signalArray > .5);
        
        if grayscale
            leftGrad = linspace(0, lightestGray, length(leftArray)+1)';
            rightGrad = linspace(lightestGray, 0, length(rightArray)+1)';
            
            % Createa color map that assumes two .5 signal levels:
            cMap = [leftGrad leftGrad leftGrad; rightGrad rightGrad rightGrad];
        else
            leftGrad = linspace(cyanMax, minAll, length(leftArray)+1)';
            leftNoGun = zeros(length(leftArray)+1, 1);
            rightGrad = linspace(minAll, magMax, length(rightArray)+1)';
            rightNoGun = zeros(length(rightArray)+1, 1);
            
            % Createa color map that assumes two .5 signal levels:
            cMap = [leftNoGun leftGrad leftGrad; rightGrad rightNoGun rightGrad];
        end
        
        
        % take out members of the color map w.r.t. how many .5 levels the
        % signal array actually has
        switch nFifty
            case 0
                cMap = [cMap(1 : length(leftArray), :); cMap(length(leftArray)+3 : end, :)];
            case 1
                cMap = [cMap(1 : length(leftArray), :); cMap(length(leftArray)+2 : end, :)];
                cMap(ceil(size(cMap, 1)/2),:) = [.15 .15 .3];
            case 2
                % do nothing
        end
end
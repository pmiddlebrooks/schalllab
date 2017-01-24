function cMap = maskbet_colormap(soaArray)

% If percentages were input instead of proportions, change it to
% proportions


switch soaArray
    case 'collapse'
        nSOA = 1;
    otherwise
        nSOA = length(soaArray);
end

targ = [0 .3 1];
dist = [1 0 .3];
high = [0 .3 1];
dist = [1 0 .3];


switch nSOA
    case 2
        cMap = [0 .33 .67; .67 0 .33];
    otherwise

        
        leftArray = soaArray(soaArray < .5);
        rightArray = soaArray(soaArray > .5);
        
        leftGrad = linspace(1, 0, length(leftArray)+1)';
        leftNoGun = zeros(length(leftArray)+1, 1);
        rightGrad = linspace(0, 1, length(rightArray)+1)';
        rightNoGun = zeros(length(rightArray)+1, 1);
        
        % Createa color map that assumes two .5 signal levels:
        cMap = [leftNoGun leftGrad leftGrad; rightGrad rightNoGun rightGrad];
        

end
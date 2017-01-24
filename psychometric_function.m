function [signalStrength, probabilityChoice1] = psychometric_function(Parameter, modelFit, signalStrengthArray, probabilityChoice1Data)

% Parameter:          stucture containing relevant parameters for the particular modelFit
% 
% modelFit:           string that determines the equation used to transform the parameters to a continuos function
% 
% signalStrengthArray: 
% probabilityChoice1Data: 
% 
% 
% 
% 
% signalStrength:     x-values of output continuous function
% probabilityChoice1: y-values of output continuous function
% 

%%
thresholdA = 250;
thresholdB = 250;
k = .1;
% signalStrength = [-.15 -.8 -.4 0 .4 .8 .15];
signalStrength = [-.15 : .005 : .15];
driftRate = k * signalStrength;


probabilityT1 = (exp(2 * driftRate * thresholdB) - 1) ./ (exp(2 * driftRate * thresholdB) - exp(-2 * driftRate * thresholdA))
figure(45445)
plot(signalStrength, probabilityT1, '-k')
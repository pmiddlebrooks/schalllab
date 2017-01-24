function Inh_Funct = weibull_curve(params, timepoints)

% Returns the Probabiliy values of the Weibull-fitted inhibition function
% Inh_Funct = gamma - ((exp(-((xdata./alpha).^beta))).*(gamma-delta))

alpha = params(1);
beta  = params(2);
gamma = params(3);
delta = params(4);

% 2. Compute Values
Inh_Funct = gamma - ((exp( -((timepoints ./ alpha) .^ beta))) .* (gamma - delta));

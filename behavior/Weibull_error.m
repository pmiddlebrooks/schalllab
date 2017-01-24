% THIS HAS TO BE IN ANOTHER FILE (STUPID)
function [f] = Weibull_error(params,xdata)
%This subfuction looks at the current data and parameters and figures out
%the sum of squares error.  

%get params
alpha = params(1); %alpha: time at which inhition function reaches 67% probability
beta  = params(2); %beta : slope
gamma = params(3); %maximum probability value
delta = params(4); %minimum probability value

%generate predictions
[f] = gamma - ((exp(-((xdata./alpha).^beta))).*(gamma-delta));

end
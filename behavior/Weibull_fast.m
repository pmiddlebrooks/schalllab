function [best_params best_fit] = Weibull_fast(xdata, ydata, weights, nFit)
%   [best_params best_fit] = Weibull(xdata,ydata,weights,nFit)
%   
%   Given data which describe points on the x and y axes, Weibull_fast uses
%   fast implimentation of simplex (nlinfit.m) to find parameters which 
%   minimize sum of squares error based on the Weibull function:
% 
%         ydata = gamma - ((exp(-((xdata./alpha).^beta))).*(gamma-delta))
% 
%   The starting parameters are set to provide good fits to inhibiton 
%   function data in a data driven way.
%   see Hanes, Patterson, and Schall. JNeurophysiol. 1998.
% 
%   Written by david.c.godlove@vanderbilt.edu 10-23-10
%   
%   INPUT:            NOTE: all input is in vector format
%       xdata       = points on the x axis. (SSDs in the case of an
%                     inhibion function)
%       ydata       = points on the y axis. (p(noncanceled|SSD) for
%                     inhibition functions)
% 
%   OPTIONAL INPUT:
%       weights     = the number of observations at each point
%       nFit       = number of times to complete fitting routine to
%                     ensure robustness of fitting
%       
%   OUTPUT:
%       best_params = four value vector containing optimum
%                     coeffecients such that:
%                     alpha = best_fit_params(1);
%                     beta  = best_fit_params(2);
%                     gamma = best_fit_params(3);
%                     delta = best_fit_params(4);
%       best_fit    = sum of squares error of xdata and ydata at the
%                     best_fit_params value.
% 
%   see also get_SSRT, ga, and gaoptimset

if nargin < 4, nFit = 100;  end
if nargin < 3, weights = []; end


seed = 1;  rand('seed',seed); % ensure reproducability

% _________________________________________________________________________
%1) specify initial params by guesstimation.

% 1a) guess on beta by fitting a straight line from min to max and
% calculating slope
last_min  = xdata(find(ydata == min(ydata),1,'last'));
first_max = xdata(find(ydata == max(ydata),1,'first'));
rise      = max(ydata) - min(ydata);
run       = (first_max - last_min)/1000; % in seconds
beta      = rise/run;

% 1b) guess on alpha by fitting a straight line from min to max and
% finding 67% probablility
if min(ydata) >= 0.67
    alpha = min(xdata);
elseif max(ydata) <= 0.67
    alpha = max(xdata);
else
    x_line = min(xdata):(max(xdata)-min(xdata))/100:max(xdata);
    y_line = min(ydata):(max(ydata)-min(ydata))/100:max(ydata);
    alpha  = min(x_line(y_line >= 0.67));
end

% 1c) guess on gamma and delta
gamma = max(ydata);
delta = min(ydata);


% _________________________________________________________________________
%2) weight Data Points if called for
if ~isempty(weights)
    x_weighted = [];
    y_weighted = [];
    for iSSD=1:length(xdata)
        CurrWeighted_x = repmat(xdata(iSSD),weights(iSSD),1);
        CurrWeighted_y = repmat(ydata(iSSD),weights(iSSD),1);
            x_weighted = [x_weighted; CurrWeighted_x];
            y_weighted = [y_weighted; CurrWeighted_y];
    end
    xdata = x_weighted;
    ydata = y_weighted;
end


% _________________________________________________________________________
%3) fit a bunch of times and take the best fitting answer (often this is
%not needed, but sometimes a fit fails)
best_fit = inf;

for ii = 1:nFit
    
    params=[alpha beta gamma delta];
    [params,r] = nlinfit(xdata,ydata,'Weibull_error',params);
    
    if sum(r) < best_fit
        best_params = params;
        best_fit = sum(r);
        
    end
    
end

end

% % THIS HAS TO BE IN ANOTHER FILE (STUPID)
% function [f] = Weibull_error(params,xdata)
% %This subfuction looks at the current data and parameters and figures out
% %the sum of squares error.  
% 
% %get params
% alpha = params(1); %alpha: time at which inhition function reaches 67% probability
% beta  = params(2); %beta : slope
% gamma = params(3); %maximum probability value
% delta = params(4); %minimum probability value
% 
% %generate predictions
% [f] = gamma - ((exp(-((xdata./alpha).^beta))).*(gamma-delta));
% 
% end
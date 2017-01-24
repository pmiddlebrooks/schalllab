function [bestFitParams,minDiscrepancyFn] = Weibull(xData, yData, weights)
%   [bestFitParams,minDiscrepancyFn] = Weibull(xData,yData,weights,pops)
%
%   Given data which describe points on the x and y axes, Weibull uses a
%   genetic algorithm approach to find parameters which minimize sum of
%   squares error based on the Weibull function:
%
%         yData = gamma - ((exp(-((xData./alpha).^beta))).*(gamma-delta))
%
%   The starting parameters and the upper and lower bounds are set to
%   provide good fits to inhibiton function data.
%   see Hanes, Patterson, and Schall. JNeurophysiol. 1998.
%
%   Written by david.c.godlove@vanderbilt.edu 10-23-10
%
%   INPUT:
%       xData              = points on the x axis. (SSDs in the case of an
%                            inhibion function)
%       yData              = points on the y axis. (p(noncanceled|SSD) for
%                            inhibition functions)
%       weights            = the number of observations at each point
%
%   OPTIONAL INPUT:
%       pops               = two value vector describing the starting
%                            number of individuals in each population and
%                            the starting number of populations. (see ga.m)
%                            default = [60 3];
%
%   OUTPUT:
%       bestFitParams    = four value vector containing optimum
%                            coeffecients such that:
%                            alpha = bestFitParams(1);
%                            beta  = bestFitParams(2);
%                            gamma = bestFitParams(3);
%                            delta = bestFitParams(4);
%       minDiscrepancyFn         = sum of squares error of xData and yData at the
%                            bestFitParams value.
%
%   see also get_SSRT, ga, and gaoptimset



if nargin < 4, pops = [60 3]; end
if nargin < 3, weights = [];  end

if length(xData) ~= length(yData)
  fprintf('Weibull: xData has length %d) and yData has lenght %d: they need to be the same', ...
    length(xData), length(yData));
  return
end

xData = reshape(xData, length(xData), 1);
yData = reshape(yData, length(yData), 1);
% Get rid of NaNs in the data
% ----------------------------
nanData         = isnan(xData) | isnan(yData);
xData(nanData)  = [];
yData(nanData)  = [];

% Sort the data
[xData, iX]     = sort(xData);
yData           = yData(iX);


% Might want to force the maximum to 1 and/or minimum to 0
MAX_TO_1_FLAG = 0;
MIN_TO_0_FLAG = 0;
if yData(end) > .9
    MAX_TO_1_FLAG = 1;
    minGamma = .9;
end
if yData(end) == 1
    minGamma = 1;
end
if yData(1) == 0
    MIN_TO_0_FLAG = 1;
    maxDelta = 0;
end


%1) specify initial param.
alpha = 200; %alpha: time at which inhition function reaches 67% probability
beta  = 1;   %beta : slope
gamma = 1;   %maximum probability value
delta = 0.5;   %minimum probability value

param=[alpha beta gamma delta]; %must be in this format for ge.m

lower_bounds = [1       1       0.5      0.0];  %bounds for parameters
upper_bounds = [1000     25      1.0      0.5];
if MAX_TO_1_FLAG
    lower_bounds = [1       1     minGamma    0.0];  %bounds for parameters
end
if MIN_TO_0_FLAG
    upper_bounds = [1000     25      1.0      maxDelta];
end

%2) weight Data Points if called for
if ~isempty(weights)
    x_weighted = [];
    y_weighted = [];
    for iSSD=1:length(xData)
        CurrWeighted_x = repmat(xData(iSSD),weights(iSSD),1);
        CurrWeighted_y = repmat(yData(iSSD),weights(iSSD),1);
        x_weighted = [x_weighted; CurrWeighted_x];
        y_weighted = [y_weighted; CurrWeighted_y];
    end
    xData = x_weighted;
    yData = y_weighted;
end

%3) set ga options
% pop_number = pops(1);%length(pop_options)=number of populations, values = size of populations
% pop_size = pops(2);  %more/larger populations means more thorough search of param space, but
% %also longer run time.  [30 30 30] is probably bare minimum.
% pop_options(1:pop_number) = pop_size;
% 
% hybrid_options=@fmincon;%run simplex after ga to refine parameters
% % ga_options=gaoptimset('PopulationSize',pop_options,'HybridFcn',hybrid_options,'display','off','UseParallel','always');
% ga_options=gaoptimset('PopulationSize',pop_options,'HybridFcn',hybrid_options,'display','off');

%4) run GA
%fit model
% [bestFitParams,minDiscrepancyFn]=ga(...
%     @(param) Weibull_error(xData,yData,param),...
%     length(param),...
%     [],[],[],[],...
%     lower_bounds,...
%     upper_bounds,...
%     [],...
%     ga_options);

    options = optimset('MaxIter', 100000,'MaxFunEvals', 100000,'useparallel','always');
    %     options = optimset('MaxIter', 100,'MaxFunEvals', 100,'useparallel','always');
        [bestFitParams minDiscrepancyFn exitflag output] = fminsearchbnd(@(param) Weibull_error(xData, yData, param),param,lower_bounds,upper_bounds,options);

        
%5)Debugging: test-plot
% alpha=bestFitParams(1);
% beta=bestFitParams(2);
% gamma=bestFitParams(3);
% delta=bestFitParams(4);
% ypred = gamma - ((exp(-((xData./alpha).^beta))).*(gamma-delta));
%
% hold on
% plot(xData,yData,'marker','o','linestyle','none')
% plot(xData,ypred,'k')



end



function discrepancyFn = Weibull_error(xData,yData,param)
%This subfuction looks at the current data and parameters and figures out
%the sum of squares error.  The genetic fitting algorithm above tries to
%find param values to minimize SSE.

%get param
alpha = param(1);
beta  = param(2);
gamma = param(3);
delta = param(4);


% Sum of squared errors method (SSE):
%generate predictions
ypred = gamma - ((exp(-((xData./alpha).^beta))).*(gamma-delta));
% % If we need a decreasing Weibull, do that here
% if mean(diff(yData)) < 0
%     ypred = 1-ypred;
% end

%compute SSE
SSE=sum((ypred-yData).^2);
discrepancyFn = SSE;


end


function [bestFitParams, maxLikelihood] = psychometric_weibull_fit(xData, yData)
%   [bestFitParams, maxLikelihood] = psychometric_weibull_fit(xData, yData, nData, pops)
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
%       nData            = the number of observations at each point
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
%       lowest_SSE         = sum of squares error of xData and yData at the
%                            bestFitParams value.
%
%   see also get_SSRT, ga, and gaoptimset

xData(isnan(xData)) = [];
yData(isnan(yData)) = [];


MAX_TO_1_FLAG = 0;
MIN_TO_0_FLAG = 0;
% if yData(end) == 1
%     MAX_TO_1_FLAG = 1;
% end
% if yData(1) == 0
%     MIN_TO_0_FLAG = 1;   
% end

% if nargin < 4, pops = [60 3]; end
% if nargin < 3, nData = [];  end

%1) specify initial params.
% alpha = 200; %alpha: time at which inhition function reaches 67% probability
% beta  = 1;   %beta : slope
alpha = .5; %alpha: time at which inhition function reaches 67% probability
beta  = 6;   %beta : slope
gamma = .9;   %maximum probability value
delta = 0.2;   %minimum probability value

param=[alpha beta gamma delta]; %must be in this format for ge.m

lower_bounds = [0       1       0.5      0.0];  %bounds for parameters
upper_bounds = [1     25      1.0      0.5];
if MAX_TO_1_FLAG
    lower_bounds = [0       1       1      0.0];  %bounds for parameters
end
if MIN_TO_0_FLAG
    upper_bounds = [1     25      1.0      0.0];
end


% %2) weight Data Points if called for
% if ~isempty(nData)
%     x_weighted = [];
%     y_weighted = [];
%     for iSSD=1:length(xData)
%         CurrWeighted_x = repmat(xData(iSSD),nData(iSSD),1);
%         CurrWeighted_y = repmat(yData(iSSD),nData(iSSD),1);
%         x_weighted = [x_weighted; CurrWeighted_x];
%         y_weighted = [y_weighted; CurrWeighted_y];
%     end
%     xData = x_weighted;
%     yData = y_weighted;
% end

% %3) set ga options
% pop_number = pops(1);%length(pop_options)=number of populations, values = size of populations
% pop_size = pops(2);  %more/larger populations means more thorough search of param space, but
% %also longer run time.  [30 30 30] is probably bare minimum.
% pop_options(1:pop_number) = pop_size;
% 
% hybrid_options=@fmincon;%run simplex after ga to refine parameters
% % ga_options=gaoptimset('PopulationSize',pop_options,'HybridFcn',hybrid_options,'display','off','UseParallel','always');
% ga_options=gaoptimset('PopulationSize',pop_options,'HybridFcn',hybrid_options,'display','off');
% 
% %4) run GA
% %fit model
% [bestFitParams, minDiscrepancyFn]=ga(...
%     @(params) Weibull_error(xData, yData, params),...
%     length(params),...
%     [],[],[],[],...
%     lower_bounds,...
%     upper_bounds,...
%     [],...
%     ga_options);


    %options = optimset('MaxIter', 1000000,'MaxFunEvals', 1000000);
    options = optimset('MaxIter', 100000,'MaxFunEvals', 100000,'useparallel','always');
    %     options = optimset('MaxIter', 100,'MaxFunEvals', 100,'useparallel','always');
        [bestFitParams minDiscrepancyFn exitflag output] = fminsearchbnd(@(param) Weibull_error(xData, yData, param),param,lower_bounds,upper_bounds,options);

        
% Re-negate the minimum to obtain maximum likelihood
maxLikelihood = -minDiscrepancyFn;
end



function discrepancyFn = Weibull_error(xData, yData, param)
%This subfuction looks at the current data and parameters and figures out
%the sum of squares error.  The genetic fitting algorithm above tries to
%find param values to minimize SSE.

%get param
alpha = param(1);
beta  = param(2);
gamma = param(3);
delta = param(4);

%generate predictions
yPredict = gamma - ((exp(-((xData./alpha).^beta))).*(gamma-delta));
% pull the values of w away from 0 and 1
yPredict = yPredict * .99 + .005;
% size(yData)
% size(yPredict)
% size(xData)
% Maximum Likelihood method (logMaxLhood):
logLhood = sum(yData.*log(yPredict) + (1-yData).*log(1-yPredict));
discrepancyFn = -logLhood;


% % % Sum of squared errors method (SSE):
% % %compute SSE
% % SSE=sum((yPredict-yData).^2);
% % discrepancyFn = SSE;
% 
% 
% % Maximum Likelihood method (logMaxLhood):
% logLhood = zeros(1, length(xData));
% for i = 1 : length(xData)
%     logLhood(i) = ...
%         log10(nchoosek(nData(i), nData(i)*yData(i))) + ...
%         yData(i) * nData(i) * log10(yPredict(i)) + ...
%         (1 - yData(i)) * nData(i) * log10(1 - yPredict(i));
% end
% % We're trying to find maximum likelihood, but our fitting algorithm is
% % trying to minimize the function- so we need to negate the likelihood to
% % find best fitting parameters.
% discrepancyFn = -sum(logLhood);

end


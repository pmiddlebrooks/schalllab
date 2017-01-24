function [best_fit_params,lowest_SSE] = ccm_chronometric_fit(t1SignalStrength, t2SignalStrength, t1RT, t2RT, weights, pops)
%   [best_fit_params,lowest_SSE] = Weibull(xdata,ydata,weights,pops)
%
%   Given data which describe points on the x and y axes, Weibull uses a
%   genetic algorithm approach to find parameters which minimize sum of
%   squares error based on the Weibull function:
%
%         ydata = gamma - ((exp(-((xdata./alpha).^beta))).*(gamma-delta))
%
%   The starting parameters and the upper and lower bounds are set to
%   provide good fits to inhibiton function data.
%   see Hanes, Patterson, and Schall. JNeurophysiol. 1998.
%
%   Written by david.c.godlove@vanderbilt.edu 10-23-10
%
%   INPUT:
%       xdata              = points on the x axis. (SSDs in the case of an
%                            inhibion function)
%       ydata              = points on the y axis. (p(noncanceled|SSD) for
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
%       best_fit_params    = four value vector containing optimum
%                            coeffecients such that:
%                            alpha = best_fit_params(1);
%                            beta  = best_fit_params(2);
%                            gamma = best_fit_params(3);
%                            delta = best_fit_params(4);
%       lowest_SSE         = sum of squares error of xdata and ydata at the
%                            best_fit_params value.
%
%   see also get_SSRT, ga, and gaoptimset


MAX_TO_1_FLAG = 0;

if nargin < 6, pops = [60 3]; end
if nargin < 5, weights = [];  end

%1) specify initial params.
A = 30;  % A: target 1 bound
B = 30;  % B: target 2 bound
k = 10;  % k: drift rate constant (drift = k * signalStrength)
t1ND = 80; % t1ND: target 1 nondecision time
t2ND = 80; % t2ND: target 2 nondecision time

params = [A, B, k, t1ND, t2ND]; %must be in this format for ge.m

t1SignalStrength = t1SignalStrength(2:end) - .5;
t2SignalStrength = t2SignalStrength(1:end-1) - .5;
t1RT = t1RT(2:end);
t2RT = t2RT(1:end-1);

if MAX_TO_1_FLAG
%     lower_bounds = [1       1       1      0.0];  %bounds for parameters
%     upper_bounds = [500     25      1.0      0.5];
else
    lower_bounds = [1       1       .001     30.0   30.0];  %bounds for parameters
    upper_bounds = [500     500     1000     150.0  150.0];
end

%2) weight Data Points if called for
if ~isempty(weights)
%     t1SSweighted = [];
%     t2SSweighted = [];
%     t1RTweighted = [];
%     t2RTweighted = [];
%     y_weighted = [];
%     for iSS=1 : length(t1SignalStrength)
%         iWeightedSS1 = repmat(t1SignalStrength(iSS),weights(iSS),1);
%         iWeightedSS2 = repmat(t2SignalStrength(iSS),weights(iSS),1);
%         CurrWeighted_y = repmat(ydata(iSS),weights(iSS),1);
%         x_weighted = [x_weighted; CurrWeighted_x];
%         y_weighted = [y_weighted; CurrWeighted_y];
%     end
%     xdata = x_weighted;
%     ydata = y_weighted;
end

%3) set ga options
pop_number = pops(1);%length(pop_options)=number of populations, values = size of populations
pop_size = pops(2);  %more/larger populations means more thorough search of param space, but
%also longer run time.  [30 30 30] is probably bare minimum.
pop_options(1:pop_number) = pop_size;

hybrid_options=@fmincon;%run simplex after ga to refine parameters
% ga_options=gaoptimset('PopulationSize',pop_options,'HybridFcn',hybrid_options,'display','off','UseParallel','always');
ga_options=gaoptimset('PopulationSize',pop_options,'HybridFcn',hybrid_options,'display','off');

%4) run GA
%fit model

[best_fit_params,lowest_SSE]=ga(...
    @(params) Weibull_error(t1SignalStrength, t2SignalStrength, t1RT, t2RT, params),...
    length(params),...
    [],[],[],[],...
    lower_bounds,...
    upper_bounds,...
    [],...
    ga_options);

%5)Debugging: test-plot
% alpha=best_fit_params(1);
% beta=best_fit_params(2);
% gamma=best_fit_params(3);
% delta=best_fit_params(4);
% ypred = gamma - ((exp(-((xdata./alpha).^beta))).*(gamma-delta));
%
% hold on
% plot(xdata,ydata,'marker','o','linestyle','none')
% plot(xdata,ypred,'k')

%get params
A = best_fit_params(1);
B = best_fit_params(2);
k = best_fit_params(3);
t1ND = best_fit_params(4);
t2ND = best_fit_params(5);

t1TimePoints = t1SignalStrength(1) : .001 : t1SignalStrength(end);
t2TimePoints = t2SignalStrength(1) : .001 : t2SignalStrength(end);
t1Drift = k * t1TimePoints;
t2Drift = k * t2TimePoints; 
t1 = t1ND + ((A + B) .* coth(t1Drift .* (A+B)) ./ t1Drift) - (B .* coth(t1Drift .* B) / t1Drift);
t2 = t2ND + ((A + B) .* coth(t2Drift .* (A+B)) ./ t2Drift) - (A .* coth(t2Drift .* A) / t2Drift);

plotFlag = 1;
if plotFlag
figure(99)
hold on;
plot(t1TimePoints, t1, 'r', t2TimePoints, t2, 'b')
plot(t1SignalStrength, t1RT, 'or')
plot(t2SignalStrength, t2RT, 'ob')
end
end



function [SSE] = Weibull_error(t1SignalStrength, t2SignalStrength, t1RT, t2RT, params)
%This subfuction looks at the current data and parameters and figures out
%the sum of squares error.  The genetic fitting algorithm above tries to
%find params values to minimize SSE.

%get params
A = params(1);
B = params(2);
k = params(3);
t1ND = params(4);
t2ND = params(5);

%generate predictions
t1Drift = k * t1SignalStrength;
t2Drift = k * t2SignalStrength; 
t1Predict = t1ND + ((A + B) .* coth(t1Drift .* (A+B)) ./ t1Drift) - (B .* coth(t1Drift .* B) / t1Drift);
t2Predict = t2ND + ((A + B) .* coth(t2Drift .* (A+B)) ./ t2Drift) - (A .* coth(t2Drift .* A) / t2Drift);

%compute SSE
SSE=sum(([t1Predict t2Predict] - [t1RT t2RT]).^2);
end


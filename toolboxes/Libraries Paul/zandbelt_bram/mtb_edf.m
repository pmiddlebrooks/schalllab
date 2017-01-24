function F = mtb_edf(t,y)
%% Empirical distribution function (edf)
%
%% Syntax
% F = mtb_edf(t,y)
%
%% Description
% |F = edf(t,y)| takes an ordered column vector |y| and a column vector of
% points |t| for which the edf is desired. It returns a column vector |F| equal
% in length to the vector |t| containing the points of the edf.
%
%% References
% Van Zandt, T. (2002). Analysis of response time distributions. In: Pashler, H.
% Stevens' handbook of experimental psychology. New York: John Wiley & Sons.
%
% This function is a modified version of Trisha van Zandt's edf.m, accessed at
% http://maigret.psy.ohio-state.edu/~trish/Downloads/matlab/EDF.m on January 6,
% 2012

% Bram Zandbelt, January 2012

% Check if t and y are column vectors
if ~iscolumn(t)
   t = t';
   if ~iscolumn(t)
      error('t should be a column vector');
   end
end

if ~iscolumn(y)
   y = y';
   if ~iscolumn(y)
      error('y should be a column vector');
   end
end

% Make sure y is ordered
if ~issorted(y)
   y = sort(y);
end

% Compute EDF
F=ones(length(t),1);
for i=1:length(t)
    F(i) = sum(y<=t(i))/length(y);
end

% Replace any nans in data with nan in F
F(isnan(y)) = nan;
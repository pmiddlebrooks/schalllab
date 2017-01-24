function F = cmtb_edf(t,y)
% CMTB_EDF Empirical distribution function
%  
% DESCRIPTION 
% Computes the cumulative distribution function associated with observed data.
%
% This function is a modified version of Trisha van Zandt's edf.m, accessed at
% http://maigret.psy.ohio-state.edu/~trish/Downloads/matlab/EDF.m on January 6,
% 2012
%  
% SYNTAX 
% F = CMTB_EDF(t,y); 
%
% t             - Mx1 ordered vector of time points for which to compute cumulative probabilities
% y             - Nx1 ordered vector of observed data
% F             - Mx1 ordered vector of cumulative probabilities
%  
% EXAMPLES
% t = 0:2000;                       % Time points
% y = sort(randi(2000,[1000,1]));   % Sample and sort integers from a uniform distribution
% F = cmtb_edf(t,y);                % Compute cumulative probabilities
% plot(t,F);
%
% REFERENCES 
% Van Zandt, T. (2002). Analysis of response time distributions. In: Pashler, H.
% Stevens' handbook of experimental psychology. New York: John Wiley & Sons.
%
% .........................................................................................................................
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Tue 11 Mar 2014 15:53:55 CDT by bram 
% $Modified: Tue 11 Mar 2014 16:32:59 CDT by bram

% Make sure that inputs are vectors
if ~all([isvector(t), isvector(y)])
  error('All inputs must be (column) vectors');
end

% Make sure that vectors are ordered columns vectors
t = sort(t(:));
y = sort(y(:));

% Compute empirical distribution function
F = arrayfun(@(inp1) sum(y<=t(inp1))/numel(y),1:numel(t));

% Replace any NaNs in data with NaN in F
F(isnan(y)) = NaN;
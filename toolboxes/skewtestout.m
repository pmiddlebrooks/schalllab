function [output] = skewtestout(x)
% PURPOSE:
%     Performs a Skewness Test on normality with the transformation of
%     the null distribution of skewness to normality after D´Agostino (1970).
% 
% USAGE:
%     results = skewtest(data)
% 
% INPUTS:
%     x - A set of data from a presumed normal distribution
% 
% OUTPUTS:
% 
% COMMENTS:
% 
%    D´Agostino (1970), Transformation to normality of the null distribution
%    of g1, Biometrika 57, p. 679-681
%
% Author: Oliver Schwindler
% oliver.schwindler@sowi.uni-bamberg.de
% Version 1.0   Date: 02/26/2005

%
% Last line in output.
%

output = 'Copyright Oliver Schwindler';

%
% Ensure the sample data is a VECTOR.
%

[rows , columns]  =  size(x);

if (rows ~= 1) & (columns ~= 1) 
    error('stats:jbtest:VectorRequired','Input sample X must be a vector.');
end

%
% Remove missing observations indicated by NaN's.
%

x  =  x(~isnan(x));

if length(x) == 0
   error('stats:jbtest:NotEnoughData',...
         'Input sample X has no valid data (all NaN''s).');
end

x  =  x(:);               % Ensure a column vector.

% Set default ALPHA levels

alpha1  =  0.05;
alpha2  =  0.025;
alpha3  =  0.005;

%
% Transformation of Skewness to normality
%

n  =  length(x);                                                 % Sample size.
x  =  (x - mean(x)); 
M2 =  sum(x.^2) / n;                                             % Sample Moment 2
M3 =  sum(x.^3) / n;                                             % Sample Moment 3
M4 =  sum(x.^4) / n;                                             % Sample Moment 4
S  =  M3 / M2.^(3/2);                                            % Sample Skewness
K  =  M4 / M2.^2;                                                % Sample Kurtosis

Y  = S*(((n +1)*(n+3))/(6*(n-2))).^(1/2);
B2 = (3*(n.^2+27*n-70)*(n+1)*(n+3))/((n-2)*(n+5)*(n+7)*(n+9));
W2 = sqrt(2*(B2-1))-1;
W  = sqrt(W2);
d  = 1/sqrt(log(W));
a  = sqrt(2/(W2-1));

b1  = d*log(Y/a+sqrt((Y/a)^2+1));

%
% Compute the P-values and critical values with the
% NORMINV function. Under the null hypothesis of composite 
% normality, the test statistic of the standardized data is 
% normal distributed.
%

cValpha1two  =  norminv([0.95]);
cValpha2two  =  norminv([0.975]);
cValpha3two  =  norminv([0.995]);

if b1 <=0
    cValpha1one  =  norminv([0.10]);
    cValpha2one  =  norminv([0.05]);
    cValpha3one  =  norminv([0.01]);
else
    cValpha1one  =  norminv([0.90]);
    cValpha2one  =  norminv([0.95]);
    cValpha3one  =  norminv([0.99]);    
end;

if b1 <= 0
    pVb1  =  normcdf(b1);
else
    pVb1  = 1 - normcdf(b1);
end;

% Test output.

fprintf('Skewness Test with D´Agostino Transformation\n');
fprintf('--------------------------------------------\n');
fprintf('          Skewness Test D´Agostino          \n');
fprintf('--------------------------------------------\n');
fprintf('Sample size:        %5.0f\n', n);
fprintf('Skewness:        %1.5f\n', S);
fprintf('Kurtosis:         %1.5f\n', K);
fprintf('X(b1):            %3.4f\n', b1);
fprintf('pV:                %3.4f\n', pVb1);
fprintf('--------------------------------------------\n');
fprintf('                 one-sided                  \n');
fprintf('--------------------------------------------\n');
fprintf('cV one-sided 1%%:      %3.4f\n', cValpha3one);
fprintf('cV one-sided 5%%:      %3.4f\n', cValpha2one);
fprintf('cV one-sided 10%%:     %3.4f\n', cValpha1one);
if b1 >=0
    if b1 >= cValpha3one;
        fprintf('The sample population is not normal\ndistributed at 1%% significance with a one-sided test.\n');
    elseif b1 >= cValpha2one;
        fprintf('The sample population is not normal\ndistributed at 5%% significance with a one-sided test.\n');
    elseif b1 >= cValpha1one;
        fprintf('The sample population is not normal\ndistributed at 10%% significance with a one-sided test.\n');
    else
        fprintf('The sample population is normal\ndistributed with a one-sided test.\n');
    end;
else
    if b1 <= cValpha3one;
        fprintf('The sample population is not normal\ndistributed at 1%% significance with a one-sided test.\n');
    elseif b1 <= cValpha2one;
        fprintf('The sample population is not normal\ndistributed at 5%% significance with a one-sided test.\n');
    elseif b1 <= cValpha1one;
        fprintf('The sample population is not normal\ndistributed at 10%% significance with a one-sided test.\n');
    else
        fprintf('The sample population is normal\ndistributed with a one-sided test.\n');
    end;
end;
fprintf('--------------------------------------------\n');
fprintf('                 two-sided                  \n');
fprintf('--------------------------------------------\n');
fprintf('cV two-sided 1%%:      %3.4f\n', cValpha3two);
fprintf('cV two-sided 1%%:     -%3.4f\n', cValpha3two);
fprintf('cV two-sided 5%%:      %3.4f\n', cValpha2two);
fprintf('cV two-sided 5%%:     -%3.4f\n', cValpha2two);
fprintf('cV two-sided 10%%:     %3.4f\n', cValpha1two);
fprintf('cV two-sided 10%%:    -%3.4f\n', cValpha1two);
if b1 >=0
    if b1 >= cValpha3two;
        fprintf('The sample population is not normal\ndistributed at 1%% significance with a two-sided test.\n');
    elseif b1 >= cValpha2two;
        fprintf('The sample population is not normal\ndistributed at 5%% significance with a two-sided test.\n');
    elseif b1 >= cValpha1two;
        fprintf('The sample population is not normal\ndistributed at 10%% significance with a two-sided test.\n');
    else
        fprintf('The sample population is normal\ndistributed with a two-sided test.\n');
    end;
else
    if abs(b1) >= cValpha3two;
        fprintf('The sample population is not normal\ndistributed at 1%% significance with a two-sided test.\n');
    elseif abs(b1) >= cValpha2two;
        fprintf('The sample population is not normal\ndistributed at 5%% significance with a two-sided test.\n');
    elseif abs(b1) >= cValpha1two;
        fprintf('The sample population is not normal\ndistributed at 10%% significance with a two-sided test.\n');
    else
        fprintf('The sample population is normal\ndistributed with a two-sided test.\n');
    end;
end; 
fprintf('--------------------------------------------\n');
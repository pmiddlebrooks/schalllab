function C = assignat(C,i,val)
% ASSIGNAT Explicit assignments of values to individual cells in an array
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% C     - cell array
% i     - indices within each cell that need to be assigned value val
% val   - value to be assigned
%
% C = ASSIGNAT(C,i,val); 
%  
% EXAMPLES 
% C = num2cell(logical(binornd(1,0.5,2,10)));
% i = 1;
% val = true;
% C = ASSIGNAT(C,i,val); 
%
% REFERENCES 
% This function was suggested by Walter Robinson on Matlab Central
% http://www.mathworks.com/matlabcentral/answers/27981-redefining-values-within-cell-array-using-cellfun
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Fri 11 Jul 2014 13:52:29 CDT by bram 
% $Modified: Fri 11 Jul 2014 13:52:29 CDT by bram 

C(i) = val;
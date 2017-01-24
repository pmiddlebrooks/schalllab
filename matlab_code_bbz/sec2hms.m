function hms = sec2hms(t)
% SEC2HMS Converts a double in seconds to a string in hrs:min:sec
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% hms = SEC2HMS(t); 
% t     - time in seconds (1x1 double)
%
% EXAMPLES 
% t   = 88000;
% hms = SEC2HMS(t)
%
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Thu 15 May 2014 09:24:09 CDT by bram 
% $Modified: Thu 15 May 2014 09:24:09 CDT by bram 

hours   = floor(t / 3600);
t       = t - hours * 3600;
mins    = floor(t / 60);
secs    = t - mins * 60;
hms     = sprintf('%02d:%02d:%05.2f', hours, mins, secs);
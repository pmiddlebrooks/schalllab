function chiSquare = sam_chi_square(pMObs,pMPrd,nObs)
% SAM_CHI_SQUARE <Synopsis of what this function does> 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% chiSquare = SAM_CHI_SQUARE(pMObs,pMPrd,fObs); 
%  
% pMObs       observed probability mass
% pMPrd       predicted probability mass
% nObs        observed frequency
%
% EXAMPLES 
%    
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Wed 18 Sep 2013 17:03:57 CDT by bram 
% $Modified: Wed 18 Sep 2013 17:06:56 CDT by bram

% chiSquare = sum((((pMObs.*fObs) - (pMPrd.*fObs)).^2) ./ (pMPrd.*fObs));
chiSquare = nObs.*sum(((pMObs - pMPrd).^2)./pMPrd);
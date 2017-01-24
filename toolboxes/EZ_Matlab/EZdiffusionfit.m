function [v,a,Ter] = EZdiffusionfit(Pc,VRT,MRT,s)

%EZdiffusionfit  --  Fit the EZ-diffusion model to P(correct) and RTs
%
% [v,a,Ter] = EZdiffusionfit(Pc,VRT,MRT)
% [v,a,Ter] = EZdiffusionfit(Pc,VRT,MRT,s)
%
% Pc is probability correct. Can be a matrix of conditions.
% VRT is the variance of response time for correct decisions (only!).
% MRT is the mean response time for correct decisions (only!).
% The error RT distribution is assumed identical to the correct RT distrib.
% size(Pc)=size(VRT)=size(MRT)=size(v)=size(a)=size(Ter). Scalars work too.
% s is a scaling parameter. Default s=0.1
% Edge corrections are required for cases with Pc=0 or Pc=1. (Pc=.5 is OK)
%
% EZdiffusionfit transforms the triad [Pc,VRT,MRT] into the triad [v,a,Ter]:
% V is a (matrix of) mean drift rate (per condition)
% A is a (matrix of) boundary separation (per condition)
% Ter is a (matrix of) non-decision time
%
% Assumptions of the EZ-diffusion model:
%  * The error RT distribution is identical to the correct RT distrib.
%  * z=a/2 -- starting point is equidistant from the response boundaries
%  * eta=0 -- across-trial variability in drift rate is negligible
%  * sz=0  -- across-trial variability in starting point is negligible
%  * st=0  -- across-trial range in nondecision time is negligible
% See Wagenmakers et al. (2007) for three EZ checks for misspecification.
%
% Example:  At end of Appendix in Wagenmakers et al. (2007)
%  [v,a,Ter] = EZdiffusionfit(.802,.112,.723)
%  -->  v = 0.0999    a = 0.1400    Ter = 0.3000
%
% Example:  Map out the correspondence [Pc,VRT] --> [v,a]
%  Pc1=[.50:.01:.99]; VRT1=[0:.010:.500]'; [Pc,VRT]=meshgrid(Pc1,VRT1);
%  [v,a,Ter]=EZdiffusionfit(Pc,VRT,1);
%  subplot(1,2,1);contourf(Pc1,VRT1,v,15);colorbar;grid on;
%  xlabel('Pc');ylabel('VRT');title('Mean drift rate v');
%  subplot(1,2,2);contourf(Pc1,VRT1,a,15);colorbar;grid on;
%  xlabel('Pc');ylabel('VRT');title('Boundary separation a');
%
% Reference:
% * Wagenmakers, E.-J., van der Maas, H. Li. J., & Grasman, R. (2007).
%   An EZ-diffusion model for response time and accuracy.
%   Psychonomic Bulletin & Review, 14 (1), 3-22.
%  
% See also EZdiffusion, dprime.

% Original coding by Alexander Petrov, Ohio State University.
% $Revision: 1.0 $  $Date: 2008-07-19 $
%
% Part of the utils toolbox version 1.2 for MATLAB version 5 and up.
% http://alexpetrov.com/softw/utils/
% Copyright (c) Alexander Petrov 1999-2008, http://alexpetrov.com
% Please read the LICENSE and NO WARRANTY statement in ../utils_license.m

% 1.0 2008-07-19 AP: Wrote it, using Appendix in Wagenmakers et al.

%-- The scaling parameter s has a ubiquitous default of 0.1
if (nargin<4) ; s = 0.1 ; end
s2 = s.^2 ;

%-- Bypass the problematic cases
idx = find(Pc>0 & Pc~=.5 & Pc<1) ;
Pci = Pc(idx) ;

%-- Allow for scalar VRT and MRT
if (numel(VRT)==1)
    VRTi = VRT ;   % and thus can divide anything
else
    VRTi = VRT(idx) ;
end

if (numel(MRT)==1)
    MRTi = MRT ;   % and thus can be subtracted from anything
else
    MRTi = MRT(idx) ;
end

%-- Log-odds
Li = log(Pci ./ (1-Pci)) ;

%-- Calculate mean drift rate v from Equation 7 in Wagenmakers et al.
xi = Li .* (Li.*Pci.^2 - Li.*Pci + Pci - .5) ./ VRTi ;
vi = sign(Pci-.5) .* s .* xi.^(.25) ;
v = NaN(size(Pc)) ;
v(idx) = vi ;

%-- Calculate boundary separation from Eq. 5
ai = s2 .* Li ./ vi ;
a = NaN(size(Pc)) ;
a(idx) = ai ;

%-- Mean decision time from Eq. 9
yi = -vi .* ai ./ s2 ;
MDTi = (ai./(2.*vi)) .* (1-exp(yi)) ./ (1+exp(yi)) ;

%-- Non-decision time is just the residual (Eq. 8)
Ter = NaN(size(Pc)) ;
Ter(idx) = MRTi - MDTi ;

%-- Handle special case Pc=.5
idx = find(Pc==.5) ;
if (~isempty(idx))
    %- Chance-level Pc implies zero drift
    v(idx) = 0 ;
    %- Boundary separation
    if (numel(VRT)==1) ; VRTi = VRT ; else VRTi = VRT(idx) ; end
    ai = s .* (24.*VRTi).^(.25) ;   % Eq. on top of p. 9, right
    a(idx) = ai ;
    %- Nondecision time
    if (numel(MRT)==1) ; MRTi = MRT ; else MRTi = MRT(idx) ; end
    Ter(idx) = MRTi - ai.^2./(4*s2) ;
end

%---   Return v, a, Ter
%%%%%% End of file EZdiffusionfit.m

function [Pc,VRT,MRT] = EZdiffusion(v,a,Ter,s)

%EZdiffusion  --  EZ-diffusion model for response time and accuracy
%
% [Pc,VRT,MRT] = EZdiffusion(v,a,Ter)
% [Pc,VRT,MRT] = EZdiffusion(v,a,Ter,s)
%
% V is a (matrix of) mean drift rate (per condition)
% A is a (matrix of) boundary separation (per condition)
% Ter is a (matrix of) non-decision time
% size(v)=size(a)=size(Ter)=size(Pc)=size(VRT)=size(MRT). Scalars work too.
% s is a scaling parameter. Default s=0.1
%
% EZdiffusion transforms the triad [Pc,VRT,MRT] into the triad [v,a,Ter]:
% Pc is probability correct.
% VRT is the variance of response times, both for correct and errors.
% MRT is the mean response time, both for correct and errors.
% The error RT distribution is identical to the correct RT distrib.
%
% Assumptions of the EZ-diffusion model:
%  * z=a/2 -- starting point is equidistant from the response boundaries
%  * eta=0 -- across-trial variability in drift rate is negligible
%  * sz=0  -- across-trial variability in starting point is negligible
%  * st=0  -- across-trial range in nondecision time is negligible
%
% Example:  At end of Appendix in Wagenmakers et al. (2007)
%  [Pc,VRT,MRT] = EZdiffusion(.1,.14,.300)
%  -->  Pc = 0.802    VRT = 0.112    MRT = 0.723
%
% Example:  Map out the correspondence [v,a] --> [Pc,VRT]
%  v1=[-.4:.025:.4]; a1=[0:.025:.4]'; [v,a]=meshgrid(v1,a1);
%  [Pc,VRT,MRT] = EZdiffusion(v,a,.300) ;
%  subplot(1,2,1);contourf(v1,a1,Pc,15);colorbar;grid on;
%  xlabel('Drift v');ylabel('Boundary a');title('Probability correct Pc');
%  subplot(1,2,2);contourf(v1,a1,sqrt(VRT),15);colorbar;grid on;
%  xlabel('Drift v');ylabel('Boundary separation a');title('Std. dev. RT');
%
% Reference:
% * Wagenmakers, E.-J., van der Maas, H. Li. J., & Grasman, R. (2007).
%   An EZ-diffusion model for response time and accuracy.
%   Psychonomic Bulletin & Review, 14 (1), 3-22.
%  
% See also EZdiffusionfit, dprime2Pcorr.

% Original coding by Alexander Petrov, Ohio State University.
% $Revision: 1.0 $  $Date: 2008-07-20 $
%
% Part of the utils toolbox version 1.2 for MATLAB version 5 and up.
% http://alexpetrov.com/softw/utils/
% Copyright (c) Alexander Petrov 1999-2008, http://alexpetrov.com
% Please read the LICENSE and NO WARRANTY statement in ../utils_license.m

% 1.0 2008-07-20 AP: Wrote it, using Wagenmakers et al (2007)

%-- The scaling parameter s has a ubiquitous default of 0.1
if (nargin<4) ; s = 0.1 ; end
s2 = s.^2 ;

%-- Probability correct
y = -v .* a ./ s2 ;
Pc = 1 ./ (1+exp(y)) ;

%-- Postpone special case v=0
idx = find(v~=0) ;
vi = v(idx) ;
yi = y(idx) ;
if (numel(a)==1) ; ai = a ; else ai = a(idx) ; end    % allow for scalar a

%-- Variance of response times -- Equation 6 in Wagenmakers et al
exp_yi = exp(yi) ;
VRTi = (2.*yi.*exp_yi - exp(2.*yi) + 1) ./ (exp_yi+1).^2 ;
VRTi = VRTi .* (ai.*s2./(2.*vi.^3)) ;
VRT = NaN(size(v)) ;
VRT(idx) = VRTi ;

%-- Mean decision time -- Equation 9
MDTi = (1-exp_yi)./(1+exp_yi) ;
MDTi = MDTi .* (ai./(2.*vi)) ;
MDT = NaN(size(v)) ;
MDT(idx) = MDTi ;

%-- Now handle special case v=0
idx = find(v==0) ;
if (~isempty(idx))
    if (numel(a)==1) ; ai = a ; else ai = a(idx) ; end
    %- Variance of response times
    VRTi = ai.^4 ./ (24*s2.^2);     % Eq. on top of p. 9, right
    VRT(idx) = VRTi ;
    %- Mean decision time
    MDTi = ai.^2 ./ (4*s2) ;
    MDT(idx) = MDTi ;
end

%-- Observable Mean RT = Mean decision time + Ter
MRT = MDT + Ter ;   % Equation 8

%---   Return Pc, VRT, MRT
%%%%%% End of file EZdiffusion.m

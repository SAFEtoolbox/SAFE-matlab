function [ y, V, Si_ex, STi_ex ] = ishigami_homma_function( x )
%
% Implements the Ishigami-Homma function, a standard benchmark function in
% the Sensitivity Analysis literature (see for instance Eq. (4.34) in 
% Saltelli et al. (2008)).
%
% Usage:
% [ y, V, Si_ex, STi_ex ] = ishigami_homma_function(x)
% 
% x = vector of inputs x(1),x(2),x(3)                        - vector (1,3)
%     x(i) ~ Unif(-pi,+pi) for all i
%
%      y = output                                                  - scalar
%      V = output variance (*)                                     - scalar
%  Si_ex = first-order sensitivity indices (*)               - vector (1,3)
% STi_ex = total-order sensitivity indices (*)               - vector (1,3)
%
% (*) = exact value computed analytically
%
% REFERENCES
%
% Saltelli et al. (2008) Global Sensitivity Analysis, The Primer, Wiley.

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/ 

a = 2 ;
b = 1 ;
y = sin(x(1)) + a * sin(x(2))^2 + b*x(3)^4*sin(x(1)) ;

% By model definition, we should get:
% VARy = VAR(Y) = 1/2 + a^2/8 + b*pi^4/5 + b^2*pi^8/18
% V1 = VAR(E(Y|X1)) = 1/2 + b*pi^4/5 + b^2*pi^8/50
% V2 = VAR(E(Y|X2)) = a^2/8
% V3 = VAR(E(Y|X3)) = 0
% V13 = VAR(E(Y|X1,X3)) = b^2*pi^4/18 - b^2*pi^8/50
% V12 = V23 = V123 = 0
% and thus:
% ST1 = S1 + S13
% ST2 = S2
% ST3 = S13
V  = 1/2+a^2/8+b*pi^4/5+b^2*pi^8/18 ;
Si_ex(1) = ( 1/2 + b*pi^4/5 + b^2*pi^8/50 )/V ;
Si_ex(2) = a^2/8 / V ;
Si_ex(3) = 0 ;
STi_ex(1) = Si_ex(1) + ( b^2*pi^8/18 - b^2*pi^8/50 )/V ;
STi_ex(2) = Si_ex(2) ;
STi_ex(3) = ( b^2*pi^8/18 - b^2*pi^8/50 )/V ;

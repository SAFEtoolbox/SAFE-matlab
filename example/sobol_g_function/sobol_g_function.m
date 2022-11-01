function [ y, V, Si ] = sobol_g_function(x,a)
%
% Implements the Sobol' g-function, a standard benchmark function in
% the Sensitivity Analysis literature (see for instance Sec. 3.6 in 
% Saltelli et al. (2008)).
%
% Usage:
% [ y, V, Si_ex ] = sobol_g_function(x,a)
% 
% x = function inputs [x(i)~Unif(0,1) for all i]             - vector (1,M)     
% a = function parameters (fixed)                            - vector (1,M)
%
%      y = output                                                  - scalar
%      V = output variance (*)                                     - scalar
%  Si_ex = first-order sensitivity indices (*)               - vector (1,3)
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

x = x(:)' ; % row vector
a = a(:)' ; % row vector

if length(a)~=length(x);
    error('x and a must be the same length');
end

g = ( abs( 4*x-2 )+a ) ./ ( 1+a ) ;
y = prod(g) ;

% By model definition:
% Vi = VAR(E(Y|Xi)) = 1 / ( 3(1+ai)^2 )
% VARy = VAR(Y) = - 1 + prod( 1 + Vi )
Vi = 1./(3*(1+a).^2)  ; 
V  = - 1 + prod( 1+Vi ) ;
Si = Vi / V ;


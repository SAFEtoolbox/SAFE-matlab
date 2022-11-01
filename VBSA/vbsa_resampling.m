function [ XA, XB, XC ] = vbsa_resampling(X)
%
% This function implements the resampling strategy needed to build the
% approximators of the first-order (main effects) and total order
% sensitivity indices (e.g. Saltelli et al. 2008; 2010).
%
% Usage:
% [ XA, XB, XC ] = vbsa_resampling(X)
%
% Input:
%  X = matrix of NX input samples                           - matrix (NX,M)
%
% Output:
% XA = first N=NX/2 rows of X                                - matrix (N,M)
% XB = last N=NX/2 rows of X                                 - matrix (N,M)
% XC = Block matrix of M 'recombinations of XA and XB      - matrix (N*M,M)
%      XC = [ XC1 ;
%             XC2 ;
%             ...
%             XCM ]
%      Each block XCi is a (N,M) matrix whose columns are all taken from XB
%      exception made for i-th, which is taken from XA.
%
% This function is meant to be used in combination with 'vbsa_indices'. 
% See help of that function for more details and examples.
%
% REFERENCES:
%
% Saltelli et al. (2008), Global Sensitivity Analysis, The Primer, Wiley.
%
% Saltelli et al. (2010), Variance based sensitivity analysis of model 
% output. Design and estimator for the total sensitivity index, Computer 
% Physics Communications, 181, 259-270.

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

[NX,M]=size(X);
if mod(NX,2)>0; NX=NX-1; fprintf('\n WARNING: input matrix X has an odd number of rows, using the first %d rows only.\n',NX); end

N  = NX/2 ;    
XA = X(  1:N  ,:) ;
XB = X(N+1:2*N,:) ;
XC = nan(N*M,M)   ;
for i=1:M
    Ci  = XB ; Ci(:,i)=XA(:,i) ;
    XC((i-1)*N+1:(i-1)*N+N,:) = Ci ;
end


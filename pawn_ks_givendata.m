function [KS,KS_dummy,YF,Fu,Fc,YY,xc,NC,XX,YU,idx_bootstrap,YUd] = pawn_ks_givendata(X,Y,n)
%
% Compute KS statistic between the unconditional and conditional output
% distributions - to be used to approximate the PAWN sensitivity indices
% according to the approximation strategy by Pianosi and Wagener (2018).
% The KS statistic for the i-th input (i=1,...,M) at the k-th conditioning
% interval (k=1,...,n) is defined as:
%
% KS(k,i) = max | F(y) - F(y| x(i) belongs to I(i,k)) |
%            y
%
% where F(y) is the unconditional output distribution, and F(y|...) is
% the conditional output distribution when input x(i) varies within the 
% interval I(i,k).
%
% Usage:
% [KS,KS_dummy,YF,Fu,Fc,YY,xc,NC,XX] = pawn_ks_givendata(X,Y,n)
%
% Input:
%         X = set of inputs samples                            - matrix (N,M)
%         Y = set of output samples                            - matrix (N,1)
%         n = number of conditional intervals                  - scalar
%
% Output:
% KS = KS statistic between unconditional and conditional      - matrix (n,M)
%      output distributions for each input factor at different
%      conditioning intervals. 
% KS_dummy = KS statistic of the dummy parameter               - scalar
% YF = values of y at which the CDFs are given                 - vector (P,1)
% Fu = values of the empirical unconditional distribution F(y) - vector (P,1)
% Fc = cell array whose element (i,k) is a                 - cell array (M,n)
%       vector (P,1) of values of the empirical 
%       conditional distribution F(y|...)
%  YY = output samples for the conditional distributions   - cell array (M,n) 
%       The element Y{i,k} is the k-th subsample 
%       of output evaluations obtained by fixing the i-th input 
%       (or group of inputs) to its k-th conditioning interval,
%       and it is a matrix of size (NC(i,k),1).
% xc = mean value of conditioning intervals                - cell array (M,1)
%      The element xc{i} is the list of 
%      subsample centers (mean points of the conditioning intervals) 
%      used to ccondition input i, and it is a matrix of size (n,1).
% NC = number of conditional output samples                - matrix (M,n)
%      for each input and each conditional interval 
% XX = input samples to derive the conditional samples     - cell array (M,n)
%      The element X{i,k} is the k-th subsample
%      of input values where the i-th input is fixed to
%      the k-th conditioning interval (while other inputs vary freely), 
%      and it is a matrix of size (NC(i,k),M)
%
% REFERENCES
%
% Pianosi, F. and Wagener, T. (2018), Distribution-based sensitivity 
% analysis from a generic input-output sample, Env. Mod. & Soft.

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

% Check inputs
if ~isnumeric(X) ; error('input ''X'' must be numeric'); end
if ~isnumeric(Y) ; error('input ''Y'' must be numeric'); end
[N,M]=size(X) ;
[N2,m]=size(Y) ;
if N~=N2; error('input ''X'' and ''Y'' must have the same number of rows'); end
if m~=1; error('input ''Y'' must be a column vector'); end
if ~isscalar(n) ; error('input ''n'' must be scalar'); end
if abs(n-round(n)); error('''n'' must be integer'); end
if n<1  ; error('input ''n'' must be a positive integer'); end

% Split sample
[YY,xc,NC,XX] = pawn_split_sample(X,Y,n);
bootsize = round(mean(mean(NC)));

% Bootrstrap from unconditional sample:
idx_bootstrap = randperm(N,bootsize) ;
YU = Y(idx_bootstrap)        ; % (bootsize,1)
% Compute empirical CDFs of unconditional and conditional samples:
[ YF, Fu, Fc  ] = pawn_cdfs(YU,YY) ;
% Compute KS among CDFs:
KS  = pawn_ks(YF,Fu,Fc)      ; % (n,M)

% Compute KS statistic for dummy parameter:
% Bootrstrap again from unconditional sample:
YUd = Y(randperm(N,bootsize))  ; % (bootsize,1)
% Compute empirical CDFs of the two unconditional samples:
[ YFd, Fud, Fcd  ] = pawn_cdfs(YU,{YUd}) ;
% Compute KS among CDFs:
KS_dummy = pawn_ks(YFd,Fud,Fcd)       ;



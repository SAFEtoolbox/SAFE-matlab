function [YY,xc,NC,XX] = pawn_split_sample(X,Y,varargin)
%
% Split a generic input-output dataset in order to create
% unconditional and conditional samples for the approximation
% of PAWN sensitivity indices (Pianosi and Wagener, 2018)
%
% [YY,xc,NC,XX] = pawn_split_sample(X,Y)
% [YY,xc,NC,XX] = pawn_split_sample(X,Y,n)
%
% Input:
%         X = set of inputs samples                          - matrix (N,M)
%         Y = set of output samples                          - matrix (N,1)
%         n = number of conditional intervals (default: 10)  - scalar
%
% Output:
%  YY = output samples for the conditional CDFs          - cell array (M,n) 
%       The element Y{i,k} is the k-th subsample 
%       of output evaluations obtained by fixing the i-th input 
%       (or group of inputs) to its k-th conditioning interval,
%       and it is a matrix of size (NC,1).
% xc = cell array of size (M,1). The element xc{i} is the list of 
%      subsample centers (mean points of the conditioning intervals) 
%      used to compute the sensitivity of factor i
%      and it is a matrix of size (n,1).
% XX = input samples to derive the conditional samples   - cell array (M,n)
%      The element X{i,k} is the k-th subsample
%      of input values where the i-th input is fixed to
%      the k-th conditioning interval (while other inputs vary freely), 
%      and it is a matrix of size (NC,M) 
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

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(X) ; error('input ''X'' must be numeric'); end
if ~isnumeric(Y) ; error('input ''Y'' must be numeric'); end
[N,M]=size(X) ;
[n,m]=size(Y) ;
if N~=n; error('input ''X'' and ''Y'' must have the same number of rows'); end
if m~=1; error('input ''Y'' must be a column vector'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
n = 10 ;

% Recover and update optional arguments:
if nargin > 2
    if ~isempty(varargin{1})
        n = varargin{1} ;
        if ~isscalar(n) ; error('input ''n'' must be scalar'); end
        if abs(n-round(n)); error('''n'' must be integer'); end
        if n<1  ; error('input ''n'' must be a positive integer'); end
    end
end

%%%%%%%%%%%%%
% Create sub-samples:
%%%%%%%%%%%%

YY = cell(M,n);
XX = cell(M,n);
xc = cell(1,M);
NC = nan(M,n) ;
for i=1:M
edges   = min(X(:,i)):(max(X(:,i))-min(X(:,i)))/n:max(X(:,i)) ;
centers = nan(n,1) ;
for k=1:n
    centers(k) = mean([edges(k),edges(k+1)]);
    if k<n
        indices = X(:,i)>=edges(k) & X(:,i) < edges(:,k+1);
    else
        indices = X(:,i)>=edges(k) & X(:,i) <= edges(:,k+1);
    end
    YY{i,k} = Y(indices)  ;
    NC(i,k) = sum(indices);
    XX{i,k} = X(indices,:) ;
    
end
xc{i}=centers;
end
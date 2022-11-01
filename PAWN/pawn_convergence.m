function [ T_m, T_lb, T_ub ] = pawn_convergence( Y, YY, stat, NUb, NCb,varargin )
%
% Analyze the convergence of the PAWN sensitivity indices (Pianosi and
% Wagener, 2015) by repeated estimation with an increasing number of 
% samples.
%
% The PAWN index is a statistic of the Kolomogorov-Smirnov statistic
% between unconditional output CDF and conditional one evaluated
% at different conditioning values. For more info, see help of
% pawn_indices.m
%
% Usage:
% [ T_m, T_lb, T_ub ] = pawn_convergence(Y,YY, stat,NUb,NCb)
% [ T_m, T_lb, T_ub ] = pawn_convergence(Y,YY,stat,NUb,NCb,idx)
% [ T_m, T_lb, T_ub ] = pawn_convergence(Y,YY,stat,NUb,NCb,idx,Nboot) 
% [ T_m, T_lb, T_ub ] = pawn_convergence(Y,YY,stat,NUb,NCb,idx,Nboot,alfa) 
%
% Input:
%   Y = output sample for estimating the unconditional CDF  - matrix (NU,P)
%  YY = output samples for the conditional CDFs          - cell array (M,n) 
%       The element Y{i,k} is the k-th subsample 
%       of output evaluations obtained by fixing the i-th input 
%       (or group of inputs) to its k-th conditioning value,
%       and it is a matrix of size (NC,P).
% stat = name of a matlab function (e.g. 'max' or                  - string
%       'median') to compute the statistic across x_ik
% NUb = sample sizes at which unconditional CDF will         - vector (1,R)
%       be estimated (note that max(NUb) must be <= NU)
% NCb = sample sizes at which conditional CDFs will          - vector (1,R)
%       be estimated (note that max(NCb) must be <= NC)
% idx = when dealing with multiple output (P>1),                   - scalar
%       this function will actually estimate the conditional and
%       unconditional CDFs of one output only. 'idx' thus is the
%       index of the output to be analyzed (i.e. the column of
%       matrix Y and YY{i,k})                            (default value: 1)
% Nboot = number of resamples for boostrapping (default: 1)        - scalar
%  alfa = significance level for the confidence intervals estimated
%         by bootstrapping (default: 0.05)                         - scalar
%
% Output:
%  T_m = mean value of the index T(i) over Nboot estimates   - matrix (R,M)
%        at different sample sizes
% T_lb = lower bound of the index T(i) over Nboot estimates  - matrix (R,M)
%        at different sample sizes
% T_ub = upper bound of the index T(i) over Nboot estimates  - matrix (R,M)
%        at different sample sizes
%
% -------------------------------------------------------------------------
% ADVANCED USAGE 
% for Regional-Response Global Sensitivity Analysis:
% -------------------------------------------------------------------------
%
% [T_m,T_lb,T_ub]=pawn_convergence(Y,YY,stat,NUb,NCb,idx,Nboot,alfa,'above',Ythres)
%
% check help of the function 'pawn_ks' for more details about this.
%
% REFERENCES
%
% Pianosi, F. and Wagener, T. (2015), A simple and efficient method 
% for global sensitivity analysis based on cumulative distribution 
% functions, Env. Mod. & Soft., 67, 1-11.

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


if ~isnumeric(Y); error('''Y'' must be a vector'); end
[ NU, P ] = size(Y) ;
if ~iscell(YY); error('''YY'' must be a cell array'); end
[ M, n ] = size(YY)  ;
if ~ischar(stat); error('''stat'' must be a string'); end
[NC ,P2]= size(YY{1}); 
if P~=P2; error('''Y'' and ''YY{i}'' must have the same number of colums'); end


if ~isnumeric(NCb); error('''NCb'' must be numeric'); end
[m,R]= size(NUb) ;
if m>1; error('''NCb'' must be a row vector' ); end
if any(NCb<=0); error('All elements of ''NCb'' must be positive' ); end
if any(abs(NCb-round(NCb))); error('All elements of ''NCb'' must be integer'); end
if any(NCb > NC); error('All elements of ''NCb'' must be lower or equal to NC (%d)',NC); end
    
if ~isnumeric(NUb); error('''NUb'' must be numeric'); end
[m,R2]= size(NUb) ;
if m>1; error('''NUb'' must be a row vector' ); end
if R2~=R; error('''NUb'' and ''NCb'' must have the same number of columns' ); end
if any(NUb<=0); error('All elements of ''NUb'' must be positive' ); end
if any(abs(NUb-round(NUb))); error('All elements of ''NUb'' must be integer'); end
if any(NUb > NU); error('All elements of ''NUb'' must be lower or equal to NU (%d)',NU); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set optional arguments to their default values:
idx=1;
Nboot=1;
alfa=0.05;
output_condition = 'allrange' ;
par = NaN ;

% Recover and update optional arguments:
if nargin > 5
    if ~isempty(varargin{1})
        idx = varargin{1} ;
        if ~isscalar(idx) ; error('input ''idx'' must be scalar'); end
        if abs(idx-round(idx)); error('''idx'' must be integer'); end
        if idx<1  ; error('input ''idx'' must be a positive integer'); end
        if idx>P  ; error('''idx'' exceeds the number of columns in Y'); end
    end
end

if nargin > 6
    if ~isempty(varargin{2})
        Nboot=varargin{2};
        if ~isscalar(Nboot); error('''Nboot'' must be scalar'); end
        if Nboot<0; error('''Nboot'' must be positive' ); end
        if abs(Nboot-round(Nboot)); error('''Nboot'' must be an integer'); end
    end
end
if nargin > 7
    if ~isempty(varargin{3})
        alfa=varargin{3};
        if ~isscalar(alfa); error('''alfa'' must be scalar'); end
        if any([alfa<0,alfa>1]); error('''alfa'' must be in [0,1]' ); end
    end
end

if nargin > 8
    if ~isempty(varargin{4})
        output_condition = varargin{4} ;
        if ~ischar(output_condition); error('''output_condition'' must be a string'); end
    end
    if isempty(varargin{5})
        error('please specify the parameters for the output condition')
    else
        par = varargin{5} ;
    end
end

%%%%%%%%%%%%%%%%%
% Compute indices
%%%%%%%%%%%%%%%%%

T_m  = nan(R,M) ;
T_lb = nan(R,M) ;
T_ub = nan(R,M) ;

for r=1:R
    
    Y_NUr  = Y(1:NUb(r),:);
    YY_NCr = cell(M,n);
    for k=1:n
        for i=1:M
            YYik=YY{i,k};
            YY_NCr{i,k}=YYik(1:NCb(r),:);
        end
    end
    
    if Nboot > 1
        [ T_m(r,:), T_lb(r,:), T_ub(r,:) ] = pawn_indices( Y_NUr,YY_NCr,stat,idx,Nboot,alfa,output_condition,par) ;
    else % no bootstrapping
        T_m(r,:) = pawn_indices(Y_NUr,YY_NCr,stat,idx,Nboot,alfa,output_condition,par) ;        
        T_lb=[];
        T_ub=[];        
    end
    
end




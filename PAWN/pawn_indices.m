function [ T_m, T_lb, T_ub, ks ] = pawn_indices(Y,YY,stat,varargin)
%
% Compute the PAWN sensitivity index T (Pianosi and Wagener, 2015).
% This is achieved by two steps: 
% - compute the Kolmogorov-Smirnov (KS) statistic between the empirical 
% unconditional CDF and the conditional CDFs at different conditioning
% values; 
% - take a statistic (e.g. the median) of the results
%
% T(i) = stat KS(k,i) = stat max | F(y) - F(y|x(i)=x_ik) |
%         k              k    y
%
% i=1,...,M (number of inputs)
% k=1,...,n (number of conditioning values for each input)
%
% Bootstrapping can be used to derive confidence bounds of T(i)
% by repeating computions over 'Nboot' resamples of the original 
% datasets. Resamples can have the same size as the original datasets or be
% reduced.
%
% Usage:
% [ T_m, T_lb, T_ub, ks ] = pawn_indices( Y, YY, stat)
% [ T_m, T_lb, T_ub, ks ] = pawn_indices( Y, YY, stat, idx)
% [ T_m, T_lb, T_ub, ks ] = pawn_indices( Y, YY, stat, idx, Nboot)
% [ T_m, T_lb, T_ub, ks ] = pawn_indices( Y, YY, stat, idx, Nboot, alfa)
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
% idx = when dealing with multiple output (P>1),                   - scalar
%       this function will actually estimate the conditional and
%       unconditional CDFs of one output only. 'idx' thus is the
%       index of the output to be analyzed (i.e. the column of
%       matrix Y and YY{i,k})                  (default: 1)
% Nboot = number of resamples for boostrapping (default: 1)        - scalar
%  alfa = significance level for the confidence intervals estimated
%         by bootstrapping (default: 0.05)                         - scalar
%
% Output:
%  T_m = mean value of the index T(i) over Nboot estimates   - vector (1,M)
% T_lb = lower bound of the index T(i) over Nboot estimates  - vector (1,M)
% T_ub = upper bound of the index T(i) over Nboot estimates  - vector (1,M)
%
% Optional output (when Nboot=1)
%   ks = value of the KS-statistics for the n conditioning    - matrix(n,M)
%        values and for the M input factors
%
% -------------------------------------------------------------------------
% ADVANCED USAGE 
% for Regional-Response Global Sensitivity Analysis:
% -------------------------------------------------------------------------
%
% [ T_m, T_lb, T_ub, ks ] = pawn_indices(Y,YY,stat,idx,Nboot,alfa,'above',Ythres)
%
% check help of the function 'pawn_ks' for more details about this.
%
%
% NOTE: If the column 'idx' of input Y or of any YY{i,k} includes any 
% NaN values, the function will identify them and exclude them from further 
% computation. A Warning message about the number of discarded NaN elements 
% (and hence the actual number of samples used for estimating the indices 
% T) will be displayed.
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
[ ~, P ] = size(Y) ;
if ~iscell(YY); error('''YY'' must be a cell array'); end
[ M, n ] = size(YY)  ;
if ~ischar(stat); error('''stat'' must be a string'); end

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
if nargin > 3
    if ~isempty(varargin{1})
        idx = varargin{1} ;
        if ~isscalar(idx) ; error('input ''idx'' must be scalar'); end
        if abs(idx-round(idx)); error('''idx'' must be integer'); end
        if idx<1  ; error('input ''idx'' must be a positive integer'); end
        if idx>P  ; error('''idx'' exceeds the number of columns in Y'); end
    end
end

if nargin > 4
    if ~isempty(varargin{2})
        Nboot=varargin{2};
        if ~isscalar(Nboot); error('''Nboot'' must be scalar'); end
        if Nboot<0; error('''Nboot'' must be positive' ); end
        if abs(Nboot-round(Nboot)); error('''Nboot'' must be an integer'); end
    end
end
if nargin > 5
    if ~isempty(varargin{3})
        alfa=varargin{3};
        if ~isscalar(alfa); error('''alfa'' must be scalar'); end
        if any([alfa<0,alfa>1]); error('''alfa'' must be in [0,1]' ); end
    end
end

if nargin > 6
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove NaN values and select output 'idx'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Y,YY]=pawn_exclude_NaN(Y,YY,idx);
% Y is a vector without NaN values
% YY is a cell matrix of size (M,n) where each element YY{i,k} is vector of
% size (NC,1) without NaN values.
NU=length(Y); % actually size of unconditional sample after the exclusion 
% of NaNs

%%%%%%%%%%%%%%%%%
% Compute indices
%%%%%%%%%%%%%%%%%

if Nboot > 1
    ks_stat = nan(Nboot,M) ;
    for j=1:Nboot
        bootsize = NU ;
        Yb = Y( floor((rand(bootsize,1)*NU+1)))   ;
        YYb = cell(M,n) ;
        for k=1:n
            for i=1:M
                ytmp = YY{i,k} ; % vector (NCik,1)
                NCik=length(YY{i,k}); % actual size of conditional sample 
                % after the exclusion of NaNs
                
% Option 1: resampling with replacement
% Creates a new sample of the same size as 'ytmp' including 'bootsize'
% components of 'ytmp' (some of which will be resampled several times)
                %bootsize = NCik ; 
                %YYb{i,k} = ytmp( floor((rand(bootsize,1)*NCik+1))) ;
                
% Option 2: resampling without replacement
% Creates a new sample of smaller size than 'ytmp' including 'bootsize'
% components of 'ytmp' (each taken once only).
% We recommend this option because it avoids inducing an 'artificial'
% upper shift of the sample CDF that would instead be induced by
% bootstrapping *with* replacement (because the same value of 'y' appears
% several times in the sample)
                bootsize = round(0.8*NCik) ; % Take 80% of total sample
                YYb{i,k} = ytmp( randperm(NCik,bootsize)) ;
                
            end
        end
        [ YFb, FUb, FCb  ] = pawn_cdfs(Yb,YYb) ;
        ks = pawn_ks(YFb, FUb, FCb,output_condition,par) ; % (n,M)
        ks_stat(j,:) = feval(stat,ks) ; % (1,M)
    end
    T_m  = mean(ks_stat,1) ;
    T_lb = sort(ks_stat) ; T_lb = T_lb(max(1,round(Nboot*alfa/2)),:);
    T_ub = sort(ks_stat) ; T_ub = T_ub(round(Nboot*(1-alfa/2)),:);
    ks=[];

else % No bootstrapping

    [ YF, FU, FC  ] = pawn_cdfs(Y,YY) ;
    ks = pawn_ks(YF, FU, FC,output_condition,par); % (n,M)
    T_m = feval(stat,ks) ; % (1,M)
    T_lb = [] ;
    T_ub = [];
end

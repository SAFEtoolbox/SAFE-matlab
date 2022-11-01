function [stat,idx,Yk,stat_sd,stat_lb,stat_ub] = RSA_indices_groups(X,Y,varargin)
%
% Computation function for Regional Sensitivity Analysis with grouping
% (as first proposed by Wagener et al., 2001).
%
% The function splits the samples in a dataset X into 'ngroups' sub-sets
% corresponding to 'ngroup' equally spaced values of Y.
% Then it assesses the distance between the CDFs of X in the 
% different datasets by a statistic (maximum or median) of the
% maximum vertical distance between the CDFs, i.e.
%
%    stat = max( max( | Fi(x) - Fj(x) | ) )
%                 x
% or
%    stat = median( max( | Fi(x) - Fj(x) | ) )
%                    x
%
% where Fi() is the CDF of X in i-the dataset and Fj() is the CDF in the
% j-th dataset.
%
% See also 'RSA_plot_groups' on how to visualize results.
%
% Basic usage: 
% [stat,idx,Yk] = RSA_indices_groups(X,Y)   
% [stat,idx,Yk] = RSA_indices_groups(X,Y,ngroup)                                                   
% [stat,idx,Yk] = RSA_indices_groups(X,Y,ngroup,flag)
%
% Input:
%         X = set of inputs samples                          - matrix (N,M)
%         Y = set of output samples                          - matrix (N,1)
% ngroup = number of groups considered (default: 10)            - scalar         
%      flag = statistic for the definition of the RSA index        - scalar
%             flag=1: stat = median (default)      
%             flag=2: stat = maximum (see 'spread' below)
%
% Output:
%      stat = distance measure between CDFs for each input   - vector (1,M)
%       idx = respective group of the samples                - vector (N,1)
%             You can easily derive the n_groups datasets {Xi} as:
%             Xi = X(idx==i);
%        Yk = range of Y in each group                - vector (ngroup+1,1)   
%
% Usage with boostrapping:
%
% [stat,idx,Yk,stat_sd,stat_lb,stat_ub] = ...
%                         RSA_indices_groups(X,Y,ngroup,flag,Nboot,alfa)
%
% Optional input:
%     Nboot = number of resamples used for boostrapping            - scalar
%             (if not specified: Nboot=0, i.e. no bootstrapping)
%      alfa = significance level for the confidence                - scalar 
%             intervals estimated by bootstrapping (default: 0.05)
%
% Additional output:
% stat_sd = standard deviation of 'stat' from bootstrapping  - vector (1,M)
% stat_lb = lower bound of 'stat' from bootstrapping         - vector (1,M)
% stat_ub = upper bound of 'stat' from bootstrapping         - vector (1,M)
%
% REFERENCES
%
% Wagener, T., Boyle, D. P., Lees, M. J., Wheater, H. S., Gupta, H. V., 
% and Sorooshian, S. (2001): A framework for development and application of 
% hydrological models, Hydrol. Earth Syst. Sci., 5, 13-26.

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
ngroup=10;
flag = 1 ;
Nboot=0;
alfa=0.05;

% Recover and update optional arguments:

if nargin > 2
    if ~isempty(varargin{1})
        ngroup = varargin{1} ;
        if ~isscalar(ngroup) ; error('input ''ngroup'' must be scalar'); end
        if abs(ngroup-round(ngroup)); error('''ngroup'' must be integer'); end
        if ngroup<=1  ; error('input ''ngroup'' must be a positive integer above 1'); end
    end
end
if nargin > 3
    if ~isempty(varargin{2})
        flag = varargin{2} ;
        if ~isnumeric(flag) ; error('input ''flag'' must be a number among 1 or 2'); end
        if ~any([1 2]==flag); error('input ''flag'' must be a number among 1 or 2'); end
    end
end
if nargin > 4
    if ~isempty(varargin{3})
        Nboot=varargin{3};
        if ~isscalar(Nboot); error('''Nboot'' must be scalar'); end
        if Nboot<0; error('''Nboot'' must be positive' ); end
        if abs(Nboot-round(Nboot)); error('''Nboot'' must be an integer'); end
    end
end
if nargin > 5
    if ~isempty(varargin{4})
        alfa=varargin{4};
        if ~isscalar(alfa); error('''alfa'' must be scalar'); end
        if any([alfa<0,alfa>1]); error('''alfa'' must be in [0,1]' ); end
    end
end

%%%%%%%%%%%%%
% RSA_indices_groups
%%%%%%%%%%%%

% Compute idx and Yk for the entire sample

% Arrange inputs and output according to ascending output values
[Y_sort, ord] = sort(Y) ;

% Define indices for splitting inputs into ngroup:
split=0:floor(N/ngroup):N;

idx=nan(N,1);

Yk=nan(ngroup+1,1);

for i=1:ngroup
   idx(ord(split(i)+1:split(i+1)))= i;
   Yk(i)=Y_sort(split(i)+1);
end

Yk(end)=Y_sort(split(end));

if Nboot>1
    bootsize=N;
    B   = floor((rand(bootsize,Nboot)*N+1));
    stat_n=nan(Nboot,M);
    for n=1:Nboot
        stat_n(n,:)=RSA_groups_compute_stat(X(B(:,n),:),Y(B(:,n)),ngroup,flag);
    end  
    stat = mean(stat_n) ;
    stat_sd = std(stat_n)  ;
    stat_sorted = sort(stat_n) ;
    stat_lb = stat_sorted(max(1,round(Nboot*alfa/2)),:) ; 
    stat_ub = stat_sorted(max(1,round(Nboot*(1-alfa/2))),:) ;
    
    else
        stat=RSA_groups_compute_stat(X,Y,ngroup,flag);
        stat_sd=[];
        stat_lb=[];
        stat_ub=[];
end


function [stat] = RSA_groups_compute_stat(X,Y,ngroup,flag)

[N,M]=size(X) ;

[Y_sort, ord] = sort(Y) ;

% Define indices for splitting inputs into ngroup:
split=0:floor(N/ngroup):N;
idx=nan(N,1);


for i=1:ngroup
   idx(ord(split(i)+1:split(i+1)))= i;
end

mvd=nan(sum(2:ngroup-1),M); % Vector of distance between the CDFs

for i=1:M
    
    % Approximate CDF of the i-th parameter for each group:
    L=length(unique(X(:,i)));
    CDF_= nan(L,ngroup);
    xx  = nan(L,ngroup);    
    
    for j=1:ngroup
            xx(:,j) = unique(sort(X(:,i))) ;
            CDF_(:,j)  = empiricalcdf(X(idx==j,i),xx(:,j))    ;
    %        [xx(:,j),CDF_(:,j)] = approximatecdfpair(X(idx==j,i),X(:,i));  
    end

    % Compute the distance between the diferent CDFs
    compt=1;
    for j=1:ngroup

        for k=1:ngroup   
            if j~=k
                mvd(compt,i)=max(abs(CDF_(:,j)-CDF_(:,k)));
                compt=compt+1;
            end
        end
    end

end

if flag==1
    stat=median(mvd);
elseif flag==2
    stat=max(mvd);
end



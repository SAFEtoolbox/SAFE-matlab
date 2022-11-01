function [stat,stat_lb,stat_ub,stat_all] = RSA_convergence_thres(X,Y,NN,varargin)
%
% This function computes and plots the sensitivity indices obtained by
% using function 'RSA_indices_thres' with an increasing number of output 
% samples
% (see help of RSA_indices_thres for details about RSA and references)
%
% Basic usage:
%
% stat = RSA_convergence_thres(X,Y,NN)
% stat = RSA_convergence_thres(X,Y,NN,threshold,flag)
%
% Usage with bootstrapping:
%
% [stat,stat_lb,stat_ub] = ...
%                 RSA_convergence_thres(X,Y,NN,threshold,flag,Nboot)
% [stat,stat_lb,stat_ub] = ...
%                 RSA_convergence_thres(X,Y,NN,threshold,flag,Nboot,alfa)
% [stat,stat_lb,stat_ub,stat_all] = ...
%                 RSA_convergence_thres(X,Y,NN,threshold,flag,Nboot,alfa)
%
% Input:
%         X = set of inputs samples                          - matrix (N,M)
%         Y = set of output samples                          - matrix (N,P)
%        NN = vector of subsample sizes at which indices     - vector (1,R)
%             will be estimated (max(NN) must not exceed N)
% threshold = threshold for output values                    - vector (1,P) 
%             (if not specified: threshold=median(Y) )   
%      flag = specify the statistic to assess the distance         - scalar
%             between CDFs                                   
%             flag=1: stat = maximum vertical distance (*)         
%             flag=2: stat = area between curves (*)
%             flag=3: stat = input range reduction (*)
%             (if not specified: flag=1 )
%     Nboot = number of resamples used for boostrapping            - scalar
%             (if not specified: Nboot=0, i.e. no bootstrapping)
%      alfa = significance level for the confidence                - scalar 
%             intervals estimated by bootstrapping (default: 0.05)
% Ouput:
%      stat = distance measure between CDFs for each input   - matrix (R,M)
%             at different sampling size (if using boostrapping, 
%             this is the mean distance over bootstrap resamples)
%   stat_lb = lower bound of 'stat' from bootstrapping       - matrix (R,M)
%             at different sampling size              [or empty if Nboot<1]
%   stat_ub = upper bound of 'stat' from bootstrapping       - matrix (R,M)
%             at different sampling size              [or empty if Nboot<1]
% stat_all = collection of 'stat' estimated for each            - cell(1,R)
%            bootstrap resample 
%            stat_all{j} is a (Nboot,M) matrix of the
%            Nboot estimates of 'stat' for the jth sample size.
%
% (*) for more information about these sensitivity measures, see help of
% function 'RSA_indices_thres.m'

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

if ~isnumeric(X) ; error('input ''X'' must be a matrix of size (N,M)'); end
if ~isnumeric(Y) ; error('input ''Y'' must be a vector of size (N,1)'); end
[N,M]=size(X) ;
[n,P]=size(Y) ;
if N~=n; error('input ''X'' and ''Y'' must have the same number of rows'); end

if ~isnumeric(NN); error('input ''NN'' must be a vector of integer numbers'); end
[n,R] = size(NN);
if n~=1; error('''NN'' must be a row vector'); end
if sum(abs(NN-round(NN))); error('all elements of ''NN'' must be integer numbers'); end
if any(NN<=0); error('all elements in ''NN'' must positive'); end
if any(diff(NN)<=0); error('elements in ''NN'' must be sorted in ascending order'); end
if NN(end)>N; error('maximum value in ''NN'' must not exceed N=%g',N); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
threshold = median(Y);
flag  =1   ;
Nboot =0   ;
alfa  =0.05;

% Recover and update optional arguments:
if nargin > 3
    if ~isempty(varargin{1})
        threshold = varargin{1} ;
        if ~isnumeric(threshold) ; error('input ''threshold'' must be a matrix of size (N,P)'); end
        [n,p]=size(threshold);
        if n~=1; error('input ''threshold'' should be a row vector'); end
        if p~=P; error('''threshold'' and ''Y'' must have the same number of columns'); end
    end
end
if nargin > 4
    if ~isempty(varargin{2})
        flag = varargin{2} ;
        if ~isscalar(flag) ; error('input ''flag'' must be scalar'); end
    end
end
if nargin > 5
    if ~isempty(varargin{3})
        Nboot=varargin{3};
        if ~isscalar(Nboot); error('''Nboot'' must be scalar'); end
        if Nboot<0; error('''Nboot'' must be positive' ); end
        if abs(Nboot-round(Nboot)); error('''Nboot'' must be an integer'); end
    end
end
if nargin > 6
    if ~isempty(varargin{4})
        alfa=varargin{4};
        if ~isscalar(alfa); error('''alfa'' must be scalar'); end
        if any([alfa<0,alfa>1]); error('''alfa'' must be in [0,1]' ); end
    end
end

%%%%%%%%%%%%%%%%%
% Compute indices
%%%%%%%%%%%%%%%%%

stat = nan(R,M) ;
stat_lb = nan(R,M) ;
stat_ub = nan(R,M) ;
stat_all=cell(1,R);     
    
R = length(NN) ;
for j=1:R
    %Xj = X(1:NN(j),:) ;
    %Yj = Y(1:NN(j),:) ;
    [Xj,idx_new]=lhcube_shrink(X,NN(j));
    Yj=Y(idx_new,:);
    if Nboot > 1
        [stat(j,:),idxb,stat_lb(j,:),stat_ub(j,:),stat_all{1,j}] = RSA_indices_thres(Xj,Yj,threshold,flag,Nboot,alfa) ;
    else
        stat(j,:) = RSA_indices_thres(Xj,Yj,threshold,flag) ;
        stat_lb = [] ;
        stat_ub = [] ;
        stat_all= [] ;
    end

end

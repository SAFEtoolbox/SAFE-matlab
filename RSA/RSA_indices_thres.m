function [stat,idxb,stat_lb,stat_ub,stat_all] = RSA_indices_thres(X,Y,varargin)
%
% Computation function for Regional Sensitivity Analysis
% (as first proposed by Spear and Hornberger (1980);
% for another application example see e.g. Sieber and Uhlenbrook (2005)).
%
% It splits the samples in a dataset X into two datasets
% ('behavioural' Xb and 'non-behavioural' Xnb)
% depending on whether the associated sample in Y satisfies 
% the condition (behavioural): 
%
%                    Y(i,j)<threshold(j)     for j=1,...,P
%
% or not (non-behavioural).
% Then assess the distance between the CDFs of Xb and Xnb by a suitable
% measure (maximum vertical distance, area between curves, etc.).
%
% See also 'RSA_plot_thres' on how to visualize results.
%
% Basic usage: 
%          [stat,idxb] = RSA_indices_thres(X,Y)
%          [stat,idxb] = RSA_indices_thres(X,Y,threshold)
%          [stat,idxb] = RSA_indices_thres(X,Y,threshold,flag)
%
% Input:
%         X = set of inputs samples                          - matrix (N,M)
%         Y = set of output samples                          - matrix (N,P)
% threshold = threshold for output values                    - vector (1,P)
%             (if not specified: threshold=median(Y) )   
%      flag = specify the statistic to assess the distance         - scalar
%             between CDFs                                   
%             flag=1: stat = maximum vertical distance (see 'mvd' below)         
%             flag=2: stat = area between curves (see 'spread' below)
%             flag=3: stat = input range reduction (see 'irr' below)
%             (if not specified: flag=1 )
% Output:
%
%    stat = distance measure between CDFs for each input     - vector (1,M)
%    idxb = indices of samples statisfying the condition     - vector (N,1)
%           You can easily derive the two datasets Xb and Xnb as:
%             Xb  = X(idxb,:)  ;
%             Xnb = X(~idxb,:) ;
%
% Usage with boostrapping:
%
%          [stat,idxb,stat_lb,stat_ub,stat_all] = ...
%                        RSA_indices_thres(X,Y,threshold,flag,Nboot)
%          [stat,idxb,stat_lb,stat_ub,stat_all] = ...
%                        RSA_indices_thres(X,Y,threshold,flag,Nboot,alfa)
%
% Optional input:
%
%     Nboot = number of resamples used for boostrapping            - scalar
%             (if not specified: Nboot=0, i.e. no bootstrapping)
%      alfa = significance level for the confidence                - scalar 
%             intervals estimated by bootstrapping (default: 0.05)%
% Output:
%
%  stat_lb = lower bound of 'stat' from bootstrapping        - vector (1,M)
%  stat_ub = upper bound of 'stat' from bootstrapping        - vector (1,M)
% stat_all = collection of 'stat' estimated for each      - matrix(Nboot,M)
%            bootstrap resample - matrix
% 
% Available measures of distance between CDFs:
%
% 1. mvd = maximum vertical difference between the two CDFs
%        (behavioural and no-behavioural).        
%        Again, the larger the mvd of an input factor,
%        the higher the sensitivity to that factor.
%        However, in contrast to ''spread'', this is an
%        absolute measure, i.e. it has meaningful value per
%        se, regardless of the units of measures of X and Y
%        In fact, by definition: 
%        - it varies between 0 and 1         
%        - if equal to 0 then the two CDFs are exactely the same  
%        - if equal to 1 then the two CDFs are 'mutually exclusive'
%          (the same value is given probability 0 by one CDF and 1 
%           by the other)
%        - the higher it is, the higher the sensitivity
%
% 2. spread = area between the inputs CDF estimated over the 
%           behavioural set, and the CDF of the non-behavioural set.
%           The larger the spread of an input factor,
%           the higher the sensitivity to that factor. 
%
% 3. irr = input range reduction, that is, the relative 
%       reduction of the input range wrt the original range
%       when considering the behavioural set only  
%       Again an absolute measure of sensitivity, in fact: 
%       - it varies between 0 and 1
%       - if equal to 0 then the range has been reduced to one point 
%       - if equal to 1 then the range has not been reduced at all   
%       - the lower it is, the higher the sensitivity to that input
%
% REFERENCES
% 
% Spear, R.C. and Hornberger, G.M. (1980). Eutrophication in peel inlet,
% II, identification of critical uncertianties via generalized sensitivity
% analysis, Water Resour. Res., 14, 43-49.
%
% Sieber, A., Uhlenbrook, S. (2005). Sensitivity analysis of a distributed
% catchment model to verify the model structure, Journal of Hydrology, 
% 310(1-4), 216-235.

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
if ~isnumeric(Y) ; error('input ''Y'' must be a matrix of size (N,P)'); end
[N,M]=size(X) ;
[n,P]=size(Y) ;
if N~=n; error('input ''X'' and ''Y'' must have the same number of rows'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
threshold = median(Y);
flag  = 1   ;
Nboot = 0   ;
alfa  = 0.05;

% Recover and update optional arguments:
if nargin > 2
    if ~isempty(varargin{1})
        threshold = varargin{1} ;
        if ~isnumeric(threshold) ; error('input ''threshold'' must be a matrix of size (N,P)'); end
        [n,p]=size(threshold);
        if n~=1; error('input ''threshold'' should be a row vector'); end
        if p~=P; error('''threshold'' and ''Y'' must have the same number of columns'); end
    end
end
if nargin > 3
    if ~isempty(varargin{2})
        flag = varargin{2} ;
        if ~isnumeric(flag) ; error('input ''flag'' must be a number among 1,2 or 3'); end
        if ~any([1 2 3]==flag); error('input ''flag'' must be a number among 1,2 or 3'); end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RSA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Nboot > 1
    
        bootsize=N ;
        B = floor((rand(bootsize,Nboot)*N+1));
        stat_all = nan(Nboot,M) ;
        idxb   = false(N,Nboot) ;
        
        for n=1:Nboot
            Xi = X(B(:,n),:) ;
            Yi = Y(B(:,n),:) ;
            %figure; scatter_plots(Xi,Yi(:,1))
            %figure; scatter_plots(Xi,Yi(:,2))
            stat_all(n,:)= compute_indices(Xi,Yi,threshold,flag) ;

        end
        % Notice that RSA may return a vector of NaNs if the bootstrap
        % resample Yi is such that the threshold condition is never
        % satisfied (or always satisfied).
        % Therefore, in the following calculations, we will make sure that
        % any row of NaNs in 'stat_all' be excluded (and display a warning
        % message).
        stat      = nanmean(stat_all) ;
        idx = ~isnan(stat_all(:,1)) ;       
        if sum(idx)<Nboot
            warning('Statistics were computed using %d bootstrap resamples instead of %d',sum(idx),Nboot)
        end
        stat_lb_sorted = sort(stat_all(idx,:)) ;
        stat_lb = stat_lb_sorted(max(1,round( sum(idx)*alfa/2)),:) ; 
        stat_ub = stat_lb_sorted(max(1,round( sum(idx)*(1-alfa/2))),:) ;
        
        % Last, let's call the function once more to obtain the vector
        % 'idxb' of the indices of behavioural input samples (needed to be
        % returned among the output arguments):
        [~,idxb] = compute_indices(X,Y,threshold,flag);
    
else
    stat_all = [] ;
    stat_lb = [] ;
    stat_ub = [] ;
    [stat,idxb] = compute_indices(X,Y,threshold,flag) ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%% built-in function %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stat,idxb] = compute_indices(X,Y,threshold,flag)

[N,M]=size(X) ;
[N,P]=size(Y) ;

idxb = sum(Y<repmat(threshold,N,1),2)==P ;
stat = nan(1,M) ;

% Define below and above subsamples:
Xb  = X(idxb,:)   ;
Nb  = size(Xb,1)  ; % number of behavioural parameterizations
Xnb = X(~idxb,:)  ;
Nnb = size(Xnb,1) ; % number of non-behavioural parameterizations

if Nb<=0
    warning('Cannot find any output value below the threshold! Try increasing the threshold value')

elseif Nnb<=0
    warning('Cannot find any output value above the threshold! Try reducing the threshold value')

else % Perform RSA
    for i=1:M
        % Empirical CDF of behavioural and non-behavioural inputs:     
        xx = unique(sort(X(:,i))) ;
        CDFb  = empiricalcdf(Xb(:,i),xx)    ;
        CDFnb = empiricalcdf(Xnb(:,i),xx)   ;       
        % Compute distance between CDFs:
        if flag==1
            stat(i) = max(abs(CDFb-CDFnb)) ;
        elseif flag==2
            stat(i) = trapz( xx, max([CDFb CDFnb],[],2) ) -trapz(xx,min([CDFb CDFnb],[],2));
        end
    end
    
    if flag==3
        % Ranges of input factors that produce ''behavioural'' output:
        xmin = min(Xb)  ;
        xmax = max(Xnb) ;
        % Compute the relative reduction wrt the original ranges of input factors
        % (as they appear in the sample ''X''):
        stat = 1 - ( xmax - xmin ) ./ ( max(X) - min(X) ) ;
    end
    
end




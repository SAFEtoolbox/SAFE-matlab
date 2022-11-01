function [S_m, S_lb, S_ub] = aggregate_boot(S, varargin)
%
% This function computes the mean and confidence intervals of the
% sensitivity indices across bootstrap resamples.
% 
% Usage:
% S_m, S_lb, S_ub = aggregate_bootstrap(S)
% S_m, S_lb, S_ub = aggregate_bootstrap(S, alfa)
% 
% Input:
%    S = array of sensitivity indices estimated for each  - matrix(Nboot,M)
%        bootstrap resample
%        or cell array of sensitivity indices          or - cell(1,R)
%        estimated for each bootstrap resample (array  
%        of R matrix(Nboot, M) where S{j} are the  
%        estimates of the sensitivity indices at the  
%        jth sample size)
%     
% Optional input:
%  alfa = significance level for the confidence           - scalar
%        intervals estimated by bootstrapping
%        (default: 0.05)
%     
% Output:
%  S_m = mean sensitivity indices across bootstrap        - matrix(R,M)
%        resamples at the different sample sizes 
% S_lb = lower bound of sensitivity indices across        - matrix(R,M)
%        bootstrap resamples at the different
%       sample sizes
% S_lb = upper bound of sensitivity indices across        - matrix(R,M)
%        bootstrap resamples at the different
%        sample sizes
%     
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
if isnumeric(S)
    [Nboot,M] = size(S);
    R = 1;  % number of sample sizes
    tmp = cell(1,1); tmp{1}=S; S=tmp ; % create cell array to simplify computations
elseif iscell(S)
    R = length(S);
    [Nboot,M] = size(S{1});
else
    error('Wrong data type for input ''S''')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if ~isempty(varargin{4})
        alfa=varargin{4};
        if ~isscalar(alfa); error('''alfa'' must be scalar'); end
        if any([alfa<0,alfa>1]); error('''alfa'' must be in [0,1]' ); end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute statistics across bootstrap resamples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise variables
S_m = ones(R,M);
S_lb = ones(R,M);
S_ub = ones(R,M);

for j=1:R % loop over sample sizes
    
    if Nboot>1
        S_m(j,:) = nanmean(S{j}) ; % bootstrap mean
        idx = ~isnan(S{j}(:,1)) ;       
        if sum(idx)<Nboot
            warning('Statistics were computed using %d bootstrap resamples instead of %d',sum(idx),Nboot)
        end
        S_lb_sorted = sort(S{j}(idx,:)) ;
        S_lb(j,:) = S_lb_sorted(max(1,round( sum(idx)*alfa/2)),:) ; % lower bound
        S_ub(j,:) = S_lb_sorted(max(1,round( sum(idx)*(1-alfa/2))),:) ; % upper bound
        
    else
        S_m(j,:) = S{j};
        warning('Nboot <= 1: lower and upper bounds could not be computed')
    end
    
end

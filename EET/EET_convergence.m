function [m_r,s_r,m_lb_r,m_ub_r,s_lb_r,s_ub_r,m_sd_r,s_sd_r,m_all,s_all] = EET_convergence(EE,rr,varargin)
%
% Compute mean and standard deviation of Elementary Effects 
% using sub-samples of the original sample 'EE' of increasing size
% (see help of EET_indices for more details about the EET and references)
%
% Basic usage:
% [ m_r, s_r ] = EET_convergence(EE,rr)
%
% Usage with optional input/output argument:
% [ m_r, s_r ] = EET_convergence(EE,rr,Nboot)
% [ m_r, s_r ] = EET_convergence(EE,rr,Nboot,alfa)
% [ m_r, s_r, m_lb_r, m_ub_r, s_lb_r, s_ub_r, m_sd_r, s_sd_r ] = ...
%                EET_convergence(EE,rr,Nboot)
% [ m_r, s_r, m_lb_r, m_ub_r, s_lb_r, s_ub_r, m_sd_r, s_sd_r, ...
%   m_all,s_all] = EET_convergence(EE,rr,Nboot)
%                
%
% Input:
% EE = matrix of 'r' elementary effects (EEs)                - matrix (r,M)
% rr = reduced number of EEs at which indices will be        - vector (1,R)
%      estimated (must be a vector of integer positive 
%      values not exceeding 'r')
% Nboot = number of resamples used for boostrapping (default:0)    - scalar 
%  alfa = significance level for the confidence intervals estimated
%         by bootstrapping (default: 0.05)                         - scalar
%
% Output:
%  m_r = mean EEs at different sampling size                 - matrix (R,M)
%  s_r = std of EEs at different sampling size               - matrix (R,M)
% m_lb_r = lower bound of mean EEs at different sampling size- matrix (R,M)
% m_ub_r = upper bound of mean EEs  "   "           "      " - matrix (R,M)
% s_lb_r = lower bound of std EEs at different sampling size - matrix (R,M)
% s_ub_r = upper bound of std EEs  "   "           "      "  - matrix (R,M)
% m_sd_r = std. of mean EEs at different sampling size       - matrix (R,M)
% s_sd_r = std. of std of EEs  "   "        "      "         - matrix (R,M)
% m_all  = collection of estimates of 'mi' for each             - cell(1,R)
%          bootstrap resample 
%          m_all{j} is a (Nboot,M) matrix of the Nboot 
%          estimates of 'mi' for the jth sample size.
% s_all  = collection of estimates of 'sigma' for each          - cell(1,R)
%          bootstrap resample 
%          s_all{j} is a (Nboot,M) matrix of the Nboot 
%          estimates of 'sigma' for the jth sample size.
%
% NOTE: If the vector EE includes any NaN values, the function will 
% identify them and exclude them from further computation. A Warning message 
% about the number of discarded NaN elements (and hence the actual number 
% of samples used for estimating mi and sigma) will be displayed.

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

if ~isnumeric(EE); error('input ''EE'' must be numerical'); end
[r,M]  = size(EE) ;

if ~isnumeric(rr); error('input ''rr'' must be a vector of integer numbers'); end
[n,R] = size(rr);
if n~=1; error('''rr'' must be a row vector'); end
if sum(abs(rr-round(rr))); error('all elements of ''rr'' must be integer numbers'); end
if any(rr<=0); error('all elements in ''rr'' must positive'); end
if any(diff(rr)<=0); error('elements in ''rr'' must be sorted in ascending order'); end
if rr(end)>r; error('maximum value in ''r'' must not exceed r=%g',r); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
Nboot=0;
alfa =0.05;
%alfa=0.1;

if nargin>2
    if ~isempty(varargin{1})
        Nboot=varargin{1};
        if ~isscalar(Nboot); error('''Nboot'' must be scalar'); end
        if Nboot<0; error('''Nboot'' must be positive' ); end
        if abs(Nboot-round(Nboot)); error('''Nboot'' must be an integer'); end
    end
end
if nargin>3
    if ~isempty(varargin{2})
        alfa=varargin{2};
        if ~isscalar(alfa); error('''alfa'' must be scalar'); end
        if any([alfa<0,alfa>1]); error('''alfa'' must be in [0,1]' ); end
    end
end

%%%%%%%%%%%%%%%%%
% Compute indices
%%%%%%%%%%%%%%%%%

m_r    = nan(R,M) ;
s_r    = nan(R,M) ;
m_sd_r = nan(R,M) ;
m_lb_r = nan(R,M) ;
m_ub_r = nan(R,M) ;
s_sd_r = nan(R,M) ;
s_lb_r = nan(R,M) ;
s_ub_r = nan(R,M) ;
m_all  = cell(1,R);
s_all  = cell(1,R);

for j=1:R
    
    idx = randperm(r,rr(j)); % randomly select a subset of rr(j) samples 
   % by using 'randperm', resampling is performed without replacement 
   % (which is our recommendation in particular for small sample size)

%    idx = randi(r,rr(j),1); % use 'randi' in place of 'randperm' for 
%    % resampling with replacement 

    EEj = EE(idx,:) ; % extract the selected samples
    
    if Nboot > 1
        
        bootsize=rr(j) ;
        B = floor((rand(bootsize,Nboot)*rr(j)+1));
        
        mi_j  = nan(Nboot,M) ;
        sigma_j = nan(Nboot,M) ;
        idx_EEj= nan(Nboot,M) ;
        for n=1:Nboot  
            [mi_j(n,:),sigma_j(n,:),idx_EEj(n,:)]=compute_indices(EEj(B(:,n),:));
        end
        
        m_r(j,:)    = mean(mi_j) ;
        m_sd_r(j,:) = std(mi_j)  ;
        m_sorted    = sort(mi_j) ;
        m_lb_r(j,:) = m_sorted(max(1,round(Nboot*alfa/2)),:)  ;
        m_ub_r(j,:) = m_sorted(round(Nboot*(1-alfa/2)),:)     ;
        m_all{1,j}  = mi_j ;
        
        s_r(j,:)    = mean(sigma_j) ;
        s_sd_r(j,:) = std(sigma_j)  ;
        s_sorted    = sort(sigma_j) ;
        s_lb_r(j,:) = s_sorted(max(1,round(Nboot*alfa/2)),:)  ;
        s_ub_r(j,:) = s_sorted(round(Nboot*(1-alfa/2)),:)     ;
        s_all{1,j}  = sigma_j ;
        
        if sum(sum(isnan(EEj)))
            fprintf('\n WARNING:')
            fprintf('\n%d NaNs were found in EE',sum(sum(isnan(EEj))))
            fprintf('\nAverage number of samples that could be used to evaluate mean ')
            fprintf('and standard deviation of elementary effects is:')
            fprintf('\nX%d:  %1.0f',[1:M;mean(idx_EEj)]);
            fprintf('\n')
        end
        % Print to screen the sensitivity indices
        fprintf('\n\t mean(EE) std(EE)\n');
        fprintf('X%d:\t %2.3f\t %2.3f\n',[ 1:M; m_r(j,:); s_r(j,:)]);
    
    else
        [m_r(j,:),s_r(j,:),idx_EEj]=compute_indices(EEj);
     
        m_sd_r = [] ;
        m_lb_r = [] ;
        m_ub_r = [] ;
        s_sd_r = [] ;
        s_lb_r = [] ;
        s_ub_r = [] ;
        m_all  = [] ;
        s_all  = [] ;
        
        if sum(sum(isnan(EEj)))
            fprintf('\n WARNING:')
            fprintf('\n%d NaNs were found in EE',sum(sum(isnan(EEj))))
            fprintf('\nAverage number of samples that could be used to evaluate mean ')
            fprintf('and standard deviation of elementary effects is:')
            fprintf('\nX%d:  %1.0f',[1:M;idx_EEj]);
            fprintf('\n')
        end
        % Print to screen the sensitivity indices
        fprintf('\n\t mean(EE) std(EE)\n');
        fprintf('X%d:\t %2.3f\t %2.3f\n',[ 1:M; m_r(j,:); s_r(j,:)]);
        
    end
end

%%%% built-in function

function [mi,sigma,idx_EE ] = compute_indices(EE)
% EE matrix (r,M)
[~,M]=size(EE);

idx_EE= nan(1,M);
mi    = nan(1,M);
sigma = nan(1,M);
    
for j=1:M
    nan_EE=isnan(EE(:,j));% find any NaN in EE
    idx_EE(j)=sum(~nan_EE);% total number of NaNs in EE
    mi(j)   = mean(abs(EE(~nan_EE,j)));% mean absolute value of EE (excluding NaNs)
    sigma(j)= std(EE(~nan_EE,j));% std of EE (excluding NaNs)     
end

% Old version:

% if Nboot > 1
%     m_sd_r = nan(R,M) ;
%     m_lb_r = nan(R,M) ;
%     m_ub_r = nan(R,M) ;
%     s_sd_r = nan(R,M) ;
%     s_lb_r = nan(R,M) ;
%     s_ub_r = nan(R,M) ;
%     m_all=cell(1,R);
%     s_all=cell(1,R);
%     for j=1:R
%         bootsize=rr(j) ;
%         B = floor((rand(bootsize,Nboot)*rr(j)+1)); 
%         mi_j = nan(Nboot,M) ;
%         sigma_j = nan(Nboot,M) ;
%         
%         for n=1:Nboot
%             mi_j(n,:)= mean(EE(B(:,n),:));
%             sigma_j(n,:)=std(EE(B(:,n),:));
%         end
%         
%         m_r(j,:) = mean(mi_j) ;
%         m_sd_r(j,:) = std(mi_j)  ;
%         m_sorted = sort(mi_j) ;
%         m_lb_r(j,:) = m_sorted(max(1,round(Nboot*alfa/2)),:)  ;
%         m_ub_r(j,:) = m_sorted(round(Nboot*(1-alfa/2)),:)     ;
%         m_all{1,j}  = mi_j;
%         
%         s_r(j,:) = mean(sigma_j) ;
%         s_sd_r(j,:) = std(sigma_j)  ;
%         s_sorted = sort(sigma_j) ;
%         s_lb_r(j,:) = s_sorted(max(1,round(Nboot*alfa/2)),:)  ;
%         s_ub_r(j,:) = s_sorted(round(Nboot*(1-alfa/2)),:)     ;
%         s_all{1,j}=sigma_j;
%     end
% else
%     m_sd_r = [] ;
%     m_lb_r = [] ;
%     m_ub_r = [] ;
%     s_sd_r = [] ;
%     s_lb_r = [] ;
%     s_ub_r = [] ;
%     m_all  = [] ;
%     s_all  = [] ;
%     for j=1:R
%         idx = randi(r,rr(j),1); % randomly select a subset of rr(j) samples
%         EEi = EE(idx,:) ; % extract the selected samples
%         % EEi = EE(1:rr(j),:) ;
%         m_r(j,:) = mean(EEi);
%         s_r(j,:) = std(EEi) ;
%     end
% end


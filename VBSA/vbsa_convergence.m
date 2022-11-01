function [ Si, STi, Si_sd, STi_sd, Si_lb, STi_lb, Si_ub, STi_ub,Si_all,STi_all ] = vbsa_convergence(Y,M,NN,varargin)
%
% This function computes the variance-based first-order indices
% (or 'main effects') and total-order ('total effects') indices 
% (Homma and Saltelli, 1996) using sub-samples of the original sample 'Y'
% of increasing size.
%
% Usage:
% [ Si, STi ] = vbsa_convergence(Y,M,NN)
% [ Si, STi ] = vbsa_convergence(Y,M,NN,Nboot)
% [ Si, STi ] = vbsa_convergence(Y,M,NN,Nboot,alfa)
% [ Si, STi, Si_sd, STi_sd, Si_lb, STi_lb, Si_ub, STi_ub ] = ...
%                              vbsa_convergence(Y,M,NN,Nboot,alfa)
% [ Si, STi, Si_sd, STi_sd, Si_lb, STi_lb, Si_ub, STi_ub, ...
%            Si_all,STi_all] = vbsa_convergence(Y,M,NN,Nboot,alfa)
%
% Input:
%  Y = [ YA; YB; YC(:) ] = vector of N*(M+2) model output samples
%  M = number of inputs                                            - scalar 
% NN = sub-sample sizes at which indices will be estimated   - vector (1,R)
%     (must be a vector of integer positive values not exceeding 'N')
% Nboot = number of resamples used for boostrapping (default:0)    - scalar 
%  alfa = significance level for the confidence intervals estimated
%         by bootstrapping (default: 0.05)                         - scalar
%
%  ** [leave empty any optional argument to get default values] ** 
%
% Output:
%    Si  = estimates of the main effects at different sampling size
%    STi = estimates of the total effects at different sampling size
% Si_sd  = standard deviation of main effects at different sampling size
% STi_sd = standard deviation of total effects at different sampling size
% Si_lb  = lower bound of main effects at different sampling size
% STi_lb = lower bound of total effects at different sampling size
% Si_ub  = upper bound of main effects at different sampling size
% STi_ub = upper bound of total effects at different sampling size
%                      - all output arguments are matrices of size (R,M)
%                        If Nboot<=1, then bootstrapping is not applied and 
%                        Si_sd,STi_sd,Si_lb,STi_lb,Si_ub,STi_ub are all NaN
% Si_all = collection of 'Si' estimated for each                - cell(1,R)
%            bootstrap resample 
%            Si_all{j} is a (Nboot,M) matrix of the Nboot 
%            estimates of 'Si' for the jth sample size.
% STi_all = collection of 'STi' estimated for each              - cell(1,R)
%            bootstrap resample 
%            STi_all{j} is a (Nboot,M) matrix of the Nboot 
%            estimates of 'STi' for the jth sample size.
% SEE ALSO vbsa_indices and vbsa_resampling
%
% NOTE: If the vector Y includes any NaN values, the function  
%       will identify them and exclude them from further computation. 
%       A Warning message about the number of discarded NaN elements (and 
%       hence the actual number of samples used for estimating Si and STi) 
%       will be displayed.
%
% EXAMPLE:
%
% fun_test  = 'ishigami_homma_function' ;
% M = 3 ;
% distrfun = 'unif' ;
% distrpar = [ -pi pi ];
% N = 1000;
% sampstrat = 'lhs' ;
% X = AAT_sampling(sampstrat,M,distrfun,distrpar,2*N);
% [ XA, XB, XC ] = vbsa_resampling(X) ;
% YA = model_evaluation(fun_test,XA);
% YB = model_evaluation(fun_test,XB);
% YC = model_evaluation(fun_test,XC);
% Y = [YA;YB;YC] ;
% NN = [ N/10:N/10:N ];
% [ Si, STi ] = vbsa_convergence(Y,M,NN);
% figure
% subplot(211); plot_convergence(Si,NN);set(gca,'YLim',[-.1,1.1])
% subplot(212); plot_convergence(STi,NN);set(gca,'YLim',[-.1,1.1])
%
% REFERENCES:
%
% Homma, T. and A., Saltelli (1996). Importance measures in global 
% sensitivity analysis of nonlinear models. 
% Reliability Engineering & System Safety, 52(1), 1-17.

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

if ~isscalar(M); error('input ''M'' must be scalar'); end
if abs(M-round(M)); error('input ''M'' must be integer'); end
if M<=0; error('input ''M'' must be positive'); end

[N,n]=size(Y) ;
N = N/(M+2) ;
if n~=1; error('input ''Y'' must be a column vector'); end
if abs(N-round(N)); error('number of rows in input ''Y'' must be a multiple of (M+2)'); end

if ~isnumeric(NN); error('input ''NN'' must be a vector of integer numbers'); end
[n,R] = size(NN);
if n~=1; error('''NN'' must be a row vector'); end
if sum(abs(NN-round(NN))); error('all elements of ''NN'' must be integer numbers'); end
if any(NN<=0); error('all elements in ''NN'' must positive'); end
if any(diff(NN)<=0); error('elements in ''NN'' must be sorted in ascending order'); end
if NN(end)>N; error('maximum value in ''NN'' must not exceed N/(M+2)=%g',N); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
Nboot=0;
alfa =0.05;

if nargin>3
    if ~isempty(varargin{1})
        Nboot=varargin{1};
        if ~isscalar(Nboot); error('''Nboot'' must be scalar'); end
        if Nboot<0; error('''Nboot'' must be positive' ); end
        if abs(Nboot-round(Nboot)); error('''Nboot'' must be an integer'); end
    end
end
if nargin>4
    if ~isempty(varargin{2})
        alfa=varargin{2};
        if ~isscalar(alfa); error('''alfa'' must be scalar'); end
        if any([alfa<0,alfa>1]); error('''alfa'' must be in [0,1]' ); end
    end
end

%%%%%%%%%%%%%%%%%
% Compute indices
%%%%%%%%%%%%%%%%%

YA = Y(1:N);
YB = Y(N+1:N+N);
YC = Y(N+N+1:end);
YC = reshape(YC,N,M)  ;
Si = nan(length(NN),M) ;
STi= nan(length(NN),M) ;
Si_sd  = nan(length(NN),M) ;
STi_sd = nan(length(NN),M) ;
Si_lb  = nan(length(NN),M) ;
Si_ub  = nan(length(NN),M) ;
STi_lb = nan(length(NN),M) ;
STi_ub = nan(length(NN),M) ;
Si_all=cell(1,length(NN));
STi_all=cell(1,length(NN));

for j = 1 : length(NN)
%    YAj=YA(1:NN(j))   ;
%    YBj=YB(1:NN(j))   ;
%    YCj=YC(1:NN(j),:) ;

    idx = randperm(N,NN(j)); % randomly select a subset of NN(j) samples 
   % by using 'randperm', resampling is performed without replacement 
   % (which is our recommendation in particular for small sample size)

%    idx = randi(N,NN(j),1); % use 'randi' in place of 'randperm' for 
%    % resampling with replacement 

    YAj=YA(idx)   ; % extract the selected samples
    YBj=YB(idx)   ;
    YCj=YC(idx,:) ;
    if Nboot>1  % and compute indices:
        [Si(j,:),STi(j,:),Si_sd(j,:),STi_sd(j,:),Si_lb(j,:),STi_lb(j,:),...
            Si_ub(j,:),STi_ub(j,:),Si_all{1,j},STi_all{1,j}] = ...
            vbsa_indices(YAj,YBj,YCj(:),Nboot,alfa) ;
    else
        [ Si(j,:),STi(j,:) ]=vbsa_indices(YAj,YBj,YCj(:));
        Si_all=[];
        STi_all=[];
    end
end


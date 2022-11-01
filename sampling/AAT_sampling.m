function X = AAT_sampling(samp_strat,M,distr_fun,distr_par,N)
%
% Generates a sample X composed of N random samples of M uncorrelated
% variables. 
%
% Usage:
% X = AAT_sampling(samp_strat,M,distr_fun,distr_par,N)
%
% Input:
% samp_strat = sampling strategy                                   - string                  
%              Options: 'rsu': random uniform                     
%                       'lhs': latin hypercube                    
%          M = number of variables                                 - scalar 
%  distr_fun = probability distribution function of each variable 
%                  - string (eg: 'unif') if all variables have the same pdf
%                  - cell array of M strings (eg:{'unif','norm'}) otherwise
%  distr_par = parameters of the probability distribution function
%                         - row vector if all input variables have the same
%                         - cell array of M vectors otherwise
%          N = number of samples                                   - scalar 
%
% Output:
%          X = matrix of samples                             - matrix (N,M)    
%              Each row is a point in the input space.  
%              In contrast to OAT_sampling, rows are not sorted in any 
%              specific order, and all elements in a row usually 
%              usually differ from the elements in the following row. 
%
% Supported PDF types:
%
%      'beta'  or 'Beta',
%      'bino'  or 'Binomial',
%      'chi2'  or 'Chisquare',
%      'exp'   or 'Exponential',
%      'ev'    or 'Extreme Value',
%      'f'     or 'F',
%      'gam'   or 'Gamma',
%      'gev'   or 'Generalized Extreme Value',
%      'gp'    or 'Generalized Pareto',
%      'geo'   or 'Geometric',
%      'hyge'  or 'Hypergeometric',
%      'logn'  or 'Lognormal',
%      'nbin'  or 'Negative Binomial',
%      'ncf'   or 'Noncentral F',
%      'nct'   or 'Noncentral t',
%      'ncx2'  or 'Noncentral Chi-square',
%      'norm'  or 'Normal',
%      'poiss' or 'Poisson',
%      'rayl'  or 'Rayleigh',
%      't'     or 'T',
%      'unif'  or 'Uniform',
%      'unid'  or 'Discrete Uniform',
%      'wbl'   or 'Weibull'.
%
% Examples:
%
% % Example 1: 2 inputs, both from Unif[0,3]
% N = 1000 ;
% M = 2 ;
% distr_fun = 'unif' ;
% distr_par = [ 0 3 ];
% samp_strat = 'lhs';
% X = AAT_sampling(samp_strat,M,distr_fun,distr_par,N);
% % Plot results:
% figure; plot(X(:,1),X(:,2),'.k')
% xlabel('x_1'); ylabel('x_2')
% 
% % Example 2: 2 inputs, one from Unif[0,3], one from Unif[1,5]
% distr_fun = 'unif' ;
% distr_par = {[ 0 3 ],[1 5]};
% X = AAT_sampling(samp_strat,M,distr_fun,distr_par,N);
% % (use above code to plot results)
% 
% % Example 3: 2 inputs, one from Unif[0,3], one from discrete, uniform in [1,5]
% distr_fun = {'unif','unid'} ;
% distr_par = {[ 0 3 ],5};
% X = AAT_sampling(samp_strat,M,distr_fun,distr_par,N);
%
% % Example 4: investigate the difference between 'rsu' and 'lhs':
% N = 100 ;
% X1 = AAT_sampling('rsu',2,'unif',[0 1],N);
% X2 = AAT_sampling('lhs',2,'unif',[0 1],N);
% figure
% subplot(121); plot(X1(:,1),X1(:,2),'ok'); title('Random Uniform')
% subplot(122); plot(X2(:,1),X2(:,2),'ok'); title('Latin Hypercube')

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

if ~ischar(samp_strat)
    error('''samp_strat'' must be a string')
end

if ~isscalar(M); error('''M'' must be scalar'); end
if M<=0; error('''M'' must be positive' ); end
if abs(M-round(M)); error('''M'' must be integer'); end

if ischar(distr_fun)
    tmp=cell(1,M);for i=1:M; tmp{i}=distr_fun; end; distr_fun = tmp;
elseif iscell(distr_fun)
    if length(distr_fun)~=M
        error('If ''distr_fun'' is a cell array, it must have M components')
    end
else
    error('Wrong data type for input ''distr_fun''')
end
if isnumeric(distr_par)
    [n1,n2] = size(distr_par) ;
    if n1==1
        tmp=cell(1,M);for i=1:M; tmp{i}=distr_par; end; distr_par = tmp;
    else
        error('If ''distr_par'' is a vector, it must be a row vector')
    end
elseif iscell(distr_par)
    if length(distr_par)~=M
        error('If ''distr_par'' is a cell array, it must have M components')
    end
else
    error('Wrong data type for input ''distr_par''')
end

if ~isscalar(N); error('''N'' must be a scalar'); end
if N<=0; error('''N'' must be positive 0' ); end
if abs(N-round(N)); error('''N'' must be integer'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniformly sample the unit square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(samp_strat,'rsu')      % Uniform sampling
    X   = rand(N,M) ; % (N,M)
elseif strcmp(samp_strat,'lhs') % Latin Hypercube sampling
    X = lhcube(N,M);
else
    error('Sampling_strategy should be either ''rsu'' or ''lhs''')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map back into the specified distribution
% by inverting the CDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:M
    pars = distr_par{i};
    np   = length(pars) ;
    name = distr_fun{i}    ;
    if any(strcmpi(name,{'chi2','exp','geo','poiss','rayl','t','unid'}))
        if np~=1; error('Input %d: Number of PDF parameters not consistent with PDF type',i);
        else
            X(:,i) = feval([name 'inv'],X(:,i),pars);
        end
        
    elseif any(strcmpi(name,{'beta','bino','ev','f','gam','logn','nbin','nct','ncx2','norm','unif','wbl'}))
        if np~=2; error('Input %d: Number of PDF parameters not consistent with PDF type',i);
        else
            X(:,i) = feval([name 'inv'],X(:,i),pars(1),pars(2));
        end
        
    elseif any(strcmpi(name,{'gev','gp','hyge','ncf'}))
        if np~=3; error('Input %d: Number of PDF parameters not consistent with PDF type',i);
        else
            X(:,i) = feval([name 'inv'],X(:,i),pars(1),pars(2),pars(3));
        end
    else
        error('Input %d: Unknown PDF type',i)
    end
end


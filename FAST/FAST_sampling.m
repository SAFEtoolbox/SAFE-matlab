function [X,s] = FAST_sampling(distr_fun,distr_par,M,varargin)
%
% Implements sampling for the Fourier Amplitude Sensitivity Test (FAST;
% Cukier et al., 1978) and returns a matrix 'X' of N input samples.
% (See also FAST_sampling_unif.m for details and references about FAST 
% sampling strategy)
%
% Usage:
% [X,s] = FAST_sampling(distr_fun,distr_par,M)
% [X,s] = FAST_sampling(distr_fun,distr_par,M,N)
% [X,s] = FAST_sampling(distr_fun,distr_par,M,N,Nharm)
% [X,s] = FAST_sampling(distr_fun,distr_par,M,N,Nharm,omega)
%
% Generates a sample X composed of N random samples of M uncorrelated
% variables. 
%
% Input:
%  distr_fun = probability distribution function of each variable 
%                  - string (eg: 'unif') if all variables have the same pdf
%                  - cell array of M strings (eg:{'unif','norm'}) otherwise
%  distr_par = parameters of the probability distribution function
%                         - row vector if all input variables have the same
%                         - cell array of M vectors otherwise
%         M = number of inputs                                     - scalar
%         N = number of samples (default is 2*Nharm*max(omega)+1   - scalar
%             which is the minimum sampling size according to         (odd)
%             Cukier et al. (1978))
%     Nharm = interference factor, i.e.the number of higher        - scalar
%             harmonics to be considered (default is 4)
%     omega = angular frequencies associated to inputs       - vector (1,M)
%         (default values computed by function 'generate_FAST_frequency.m')
%
% Output:
%     X = matrix of input samples                            - matrix (N,M)
%     s = vector of sampled points over the search curve     - vector (N,1)

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

if ~isscalar(M); error('''M'' must be a scalar'); end
if (M - round(M))~=0 ; error('''M'' must be an integer'); end
if M<=0 ; error('''M'' must be a positive integer'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
Nharm = 4;  % taken from Saltelli et al. (1999; page 42)
omega = generate_FAST_frequency(M);
N     = [];

% Recover and update optional arguments:
if nargin > 3
    if ~isempty(varargin{1})
        N=varargin{1} ;
        if ~isscalar(N); error('''N'' must be scalar'); end
        if (N - round(N))~=0 ; error('''N'' must be integer'); end
        if N<=0 ; error('''N'' must be positive'); end
        if mod(N,2)==0 ; error('''N'' must be odd'); end
    end
end
if nargin > 4
    if ~isempty(varargin{2})
        Nharm = varargin{2};
        if ~isscalar(Nharm); error('''Nharm'' must be a scalar'); end
        if (Nharm - round(Nharm))~=0 ; error('''Nharm'' must be an integer'); end
        if Nharm<=0 ; error('''Nharm'' must be a positive integer'); end
    end
end
if nargin > 5
    if ~isempty(varargin{3})
        omega = varargin{3};
        if ~isnumeric(omega) ; error('input ''omega'' must be a vector of size (1,M)'); end
        [tmp,M2] = size(omega);
        if tmp>1 ; error('input ''omega'' must be a vector of size (1,M)'); end
        if M2~=M ; error('input ''omega'' must be a vector of size (1,M)'); end
        if any(omega<0); error('all components of ''omega'' should be positive'); end
        if any(omega-round(omega)~=0); error('all components of ''omega'' should be integer'); end
    end
end

if isempty(N) ; % If user did not specify the sample size
    N = 2*Nharm*max(omega)+1 ; % ... set it to the minimum sample size
end
if N < 2*Nharm*max(omega)+1 % and finally check that is is consistent with omega
    Nuser = N ;             % (in case omega was specified by user)
    N = 2*Nharm*max(omega)+1 ;
    warning('Sample size specified by user (%d) is smaller than minimum sample size. Using the latter (%d) instead',Nuser,N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniformly sample the unit square using FAST sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,s] = FAST_sampling_unif(M,N,Nharm,omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map back into the specified distribution
% by inverting the CDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:M
    pars = distr_par{i};
    np   = length(pars) ;
    name = distr_fun{i}    ;
    if any(strcmpi(name,{'chi2','exp','geo','poiss','rayl','tinv','unid'}))
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
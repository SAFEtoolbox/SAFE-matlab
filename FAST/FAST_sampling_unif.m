function [X,s] = FAST_sampling_unif(M,varargin)
%
% Implements sampling for the Fourier Amplitude Sensitivity Test (FAST;
% Cukier et al., 1978) and returns a matrix 'X' of N input samples.
% Inputs are assumed to be uniformly distributed in the unit hypercube
% [0,1]^M. 
% Samples are taken along the search curve defined by transformations
% 
%      x_i(s) = G_i( sin( omega_i*s ) )      i=1,...,M         (*)
%
% where 's' is a scalar variable that varies in (-pi/2,pi/2)
%
% Usage:
% [X,s] = FAST_sampling_unif(M)
% [X,s] = FAST_sampling_unif(M,N)
% [X,s] = FAST_sampling_unif(M,N,Nharm)
% [X,s] = FAST_sampling_unif(M,N,Nharm,omega)
%
% Input:
%     M = number of inputs                                         - scalar
%     N = number of samples (default is 2*Nharm*max(omega)+1       - scalar
%         which is the minimum sampling size according to             (odd)
%         Cukier et al. (1978))
% Nharm = interference factor, i.e.the number of higher            - scalar
%         harmonics to be considered (default is 4)
% omega = angular frequencies associated to inputs           - vector (1,M)
%         (default values computed by function 'generate_FAST_frequency.m')
%
% Output:
%     X = matrix of input samples                            - matrix (N,M)
%     s = vector of sampled points over the search curve     - vector (N,1)
%
% Notes:
% (*)  Here we use the curve proposed by Saltelli et al. (1999):
%          x_i = 1/2 + 1/pi * arcsin( sin( omega_i*s ) )
%
% References:
%
% Cukier, R.I., Levine, H.B., and Shuler, K.E. (1978), Nonlinear
% Sensitivity Analyis of Multiparameter Model SYstems, Journal of
% Computational Physics, 16, 1-42.
%
% Saltelli, A., Tarantola, S. and Chan, K.P.S. (1999), A Quantitative
% Model-Independent Method ofr Global Sensitivty Analysis of Model Output,
% Technometrics, 41(1), 39-56.

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

if ~isscalar(M); error('''M'' must be a scalar'); end
if (M - round(M))~=0 ; error('''M'' must be an integer'); end
if M<=0 ; error('''M'' must be a positive integer'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
Nharm = 4;  % taken from Saltelli et al. (1999; page 42)
omega = generate_FAST_frequency(M);
N     = [] ; % minimum sample size

% Recover and update optional arguments:
if nargin > 1
    if ~isempty(varargin{1})
        N = varargin{1};
        if ~isscalar(N); error('''N'' must be scalar'); end
        if (N - round(N))~=0 ; error('''N'' must be integer'); end
        if N<=0 ; error('''N'' must be positive'); end
        if mod(N,2)==0 ; error('''N'' must be odd'); end
    end
end
if nargin > 2
    if ~isempty(varargin{2})
        Nharm = varargin{2};
        if ~isscalar(Nharm); error('''Nharm'' must be a scalar'); end
        if (Nharm - round(Nharm))~=0 ; error('''Nharm'' must be an integer'); end
        if Nharm<=0 ; error('''Nharm'' must be a positive integer'); end
    end
end
if nargin > 3
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform sampling over the search curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = pi/2*(2*[1:N]'-N-1)/N ;
% s = [ pi/2*(-N+1)/N
%       pi/2*(-N+3)/N
%       ...
%       pi/2*( -1 )/N
%       pi/2*( +1 )/N
%       ...
%       pi/2*(+N-3)/N
%       pi/2*(+N-1)/N  ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map back sampled points in the input space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = nan(N,M) ;
for n=1:N
    % X(n,:) = feval(searchcurve,s(n),omega) ;
    X(n,:) = 1/2+1/pi*asin( sin(omega*s(n)) ) ;
end


%function x = Saltelli1999_searchcurve(s,omega)
%
%x = 1/2+1/pi*arcsin( sin(omega*s) ) ;



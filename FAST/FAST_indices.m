function [Si,V,A,B,Vi] = FAST_indices(Y,M,varargin)
%
% Computes main effect (first-order) sensitivity index
% according to the Fourier Amplitude Sensitivity Test (FAST)
% (Cukier et al., 1978; Saltelli et al., 1999)
%
% Usage:
% [Si,V,A,B,Vi] = FAST_indices(Y,M)
% [Si,V,A,B,Vi] = FAST_indices(Y,M,Nharm)
% [Si,V,A,B,Vi] = FAST_indices(Y,M,Nharm,omega)
%
% Input:
%     Y = set of model output samples                       - vector (N,1)
%     M = number of inputs                                  - scalar
% Nharm = interference factor, i.e.the number of higher     - scalar
%         harmonics to be considered (default is 4)
% omega = angular frequencies associated to inputs          - vector (1,M)
%         (default values computed by function 'generate_FAST_frequency.m')
%
% Output:
%    Si = main effect (first-order) sensitivity indices     - vector (1,M)
%     V = total output variance                             - scalar
%     A = Fourier coefficients                              - vector (1,N)
%     B = Fourier coefficients                              - vector (1,N)
%    Vi = output variances from each input                  - vector (1,M)
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

if ~isnumeric(Y) ; error('input ''Y'' must be a vector of size (N,1)'); end
[N,tmp] = size(Y) ;
if tmp>1 ; error('input ''Y'' must be a vector of size (N,1)'); end
if N<1   ; error('input ''Y'' must be a vector of size (N,1)'); end

if ~isscalar(M); error('''M'' must be a scalar'); end
if (M - round(M))~=0 ; error('''M'' must be an integer'); end
if M<=0 ; error('''M'' must be a positive integer'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
Nharm = 4;  % taken from Saltelli et al. (1999; page 42)
omega = generate_FAST_frequency(M);

% Recover and update optional arguments:
if nargin > 2
    if ~isempty(varargin{1})
        Nharm = varargin{1};
        if ~isscalar(Nharm); error('''Nharm'' must be a scalar'); end
        if (Nharm - round(Nharm))~=0 ; error('''Nharm'' must be an integer'); end
        if Nharm<=0 ; error('''Nharm'' must be a positive integer'); end
    end
end
if nargin > 3
    if ~isempty(varargin{2})
        omega = varargin{2};
        if ~isnumeric(omega) ; error('input ''omega'' must be a vector of size (1,M)'); end
        [tmp,M2] = size(omega);
        if tmp>1 ; error('input ''omega'' must be a vector of size (1,M)'); end
        if M2~=M ; error('input ''omega'' must be a vector of size (1,M)'); end
        if any(omega<0); error('all components of ''omega'' should be positive'); end
        if any(omega-round(omega)~=0); error('all components of ''omega'' should be integer'); end
    end
end
   
if N < 2*Nharm*max(omega)+1 % and finally check that is is consistent with omega
    error('Sample size (i.e. the length of vector Y) is %d, which is lower than the minimum sample size (i.e. 2*Nharm*max(omega)+1=2*%d*%d+1=%d',N,Nharm,max(omega),2*Nharm*max(omega)+1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Fourier coefficients (vectors A and B) from Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The code below implements the equations given in Appendix C
% of Saltelli et al. (1999).

A = zeros(1,N) ;
B = zeros(1,N) ;
baseplus  = sum( reshape(Y(2:end),2,(N-1)/2) )'   ; % ((N-1)/2,1)
baseminus = -diff( reshape(Y(2:end),2,(N-1)/2) )' ; % ((N-1)/2,1)    
for j=1:N
    if mod(j,2)==0 % j is even
        sp = Y(1) ;
        for k=1:(N-1)/2
            sp = sp + baseplus(k)*cos(j*k*pi/N) ;
        end
        A(j) = sp/N ;
    else           % j is odd
        sp = 0 ;
        for k=1:(N-1)/2
            sp = sp + baseminus(k)*sin(j*k*pi/N) ;
        end
        B(j) = sp/N ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute main effect from A and B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The code below implements the equations given in Appendix B
% of Saltelli et al. (1999) (here we use 'V' and 'Vi' for the output
% variances while in that paper they are called 'D' and 'Di')

V = 2*sum(A.^2+B.^2); % total output variance
Vi = nan(1,M) ; % output variances from the i-th input
for i=1:M
     idx = [1:Nharm]*omega(i) ;
     Vi(i) = 2*sum(A(idx).^2+B(idx).^2) ;
end
Si = Vi / V ; 

fprintf('\n \t main \n');
fprintf(' X%d:\t %1.4f\t \n',[ 1:M; Si ]);
fprintf('\n sum:\t %1.4f\t \n\n',sum(Si));

    
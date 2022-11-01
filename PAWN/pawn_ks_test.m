function [rej, ks_crit] = pawn_ks_test(ks,NC,NU,alfa)
%
% Apply the two-sample Kolmogorov-Smirnov test for inputs screening,
% according to the PAWN approach (Pianosi adn Wagener, 2015).
%
% [rej, ks_crit] = pawn_ks_test(ks,NC,NU,alfa,xc,dimplot)
%
% Input:
%   KS = Kolmogorov-Smirnov statistic                        - matrix (n,M)
%        between empirical unconditional CDF of the output
%        and conditional CDF when fixing the i-th input (i=1,...,M)
%        at its k-th (k=1,...,n) conditioning value 
%   NC = number of output samples used to compute the              - scalar
%        conditional CDFs
%   NU = number of output samples used to compute the              - scalar 
%        unconditional CDF
% alfa = confidence level of the KS test                           - scalar
%        must be one among {0.1,0.05,0.025,0.01,0.005,0.001}
%
% Output:
%     rej = frequency (over n conditioning values)           - vector (M,1)
%           of rejecting the null hypothesis that the
%           unconditional and conditional CDFs be equal,
%           i.e. that the i-th input be non-influential.
%           In other words, rej(i)=0 means that the
%           i-th input was always found non-influential
%           while rej(i)=1 means that it was never
%           found non-influential
% ks_crit = critical value of the KS statistic used to             - scalar
%           reject the hypothesis
%           (computed according to Table AIX in Wall (1996))
%
% REFERENCES
%
% Pianosi, F. and Wagener, T. (2015), A simple and efficient method 
% for global sensitivity analysis based on cumulative distribution 
% functions, Env. Mod. & Soft., 67, 1-11.
%
% Wall, J. (1996). Practical statistics for astronomers - ii. correlation, 
% data-modelling and sample comparison. Quarterly Journal of the 
% Royal Astronomical Society 37, 519?563.

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

if ~isnumeric(ks); error('''ks'' must be a matrix'); end
[n,M] = size(ks) ;
ks = ks' ;

if ~isscalar(NC); error('''NC'' must be a scalar'); end
if NC<=0; error('''NC'' must be positive' ); end
if abs(NC-round(NC)); error('''NC'' must be integer'); end

if ~isscalar(NU); error('''NU'' must be a scalar'); end
if NU<=0; error('''NU'' must be positive' ); end
if abs(NU-round(NU)); error('''NU'' must be integer'); end

if     alfa==0.10 ; c=1.22;
elseif alfa==0.05 ; c=1.36;
elseif alfa==0.025; c=1.48;
elseif alfa==0.01 ; c=1.63;
elseif alfa==0.005; c=1.73;
elseif alfa==0.001; c=1.95;
else
    error('''alfa'' must take a value in {0.1,0.05,0.025,0.01,0.005,0.001}');
end
% c    = [ 1.22 1.36 1.48  1.63 1.73  1.95 ] ; 
%   alfa = 0.10 0.05 0.025 0.01 0.005 0.001


%%%%%%%%%%%%%%%
% Apply KS-test
%%%%%%%%%%%%%%%

ks_crit = c*sqrt((NC+NU)/(NC*NU)) ;
reject = ks > ks_crit ;
rej    = sum(reject,2)/n ;
% frequency (across fixed points) of passing the KS test at confidence level alfa, 
% i.e. the empirical unconditional and conditional CDFs are different
% or, the factor is influential.

% Print results to screen:
str='X%d:'; for k=1:n; str = [ str ' \t %2.3f' ]; end; str = [ str ' \n'];
fprintf('X%d: \t %2.3f\n',[ 1:M; rej' ]);

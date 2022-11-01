function Fi = empiricalcdf(x,xi)
%
% Fi = empiricalcdf(x,xi)
%
% Compute the empirical CDF of the sample 'x' 
% and evaluate it at datapoints 'xi'.
% 
%  x = samples to build the empirical CDF: F(x) - vector (N,1)
% xi = values where to evaluate the CDF        - vector (Ni,1)
% Fi = CDF values at 'xi'                      - vector (Ni,1)
%
% Use: 
% F = approximatecdf(x,x)
% to obtain the CDF values at the same datapoints used for its
% construction.
%
% Example:
%
% x = rand(50,1) ;
% F = empiricalcdf(x,x) ;
% xi = [min(x):0.001:max(x)] ;
% Fi = empiricalcdf(x,xi) ;
% figure; plot(xi,Fi,'k',x,F,'or')

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

if ~isnumeric(x)
    error('input ''x'' must be numeric')
end
if ~isnumeric(xi)
    error('input ''xi'' must be numeric')
end

% Estimate empirical CDF values at 'x':
x = x(:)      ;
N = length(x) ;
x = sort(x)   ;
F = [1:N]/N   ;
% Remove any multiple occurance of 'x'
% and set F(x) to the upper value (recall that F(x) is the percentage of
% samples whose value is lower than *or equal to* x!)
[x,iu] = unique(x,'last');
F      = F(iu);
N      = length(F) ;

% Interpolate the empirical CDF at 'xi':
xi = xi(:) ;
Fi = ones(size(xi)) ;
for j=N:-1:1
    Fi(xi<=x(j)) = F(j) ;
end
Fi(xi<x(1)) = 0 ;

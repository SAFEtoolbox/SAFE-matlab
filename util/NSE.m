function nse = NSE(y_sim,y_obs)
%
% Computes the Nash-Sutcliffe Efficiency (NSE) coefficient
%
% nse = NSE(Y_sim,y_obs)
% 
% Y_sim = time series of modelled variable     - matrix (N,T)
%         (N>1 different time series can be evaluated at once)
% y_obs = time series of observed variable     - vector (1,T)
%
% nse   = vector of NSE coefficients           - vector (N,1)
%
% References:
%
% Nash, J. E. and J. V. Sutcliffe (1970), 
% River flow forecasting through conceptual models part I 
% A discussion of principles, Journal of Hydrology, 10 (3), 282-290.

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

[N,T] = size(y_sim) ;
[M,D] = size(y_obs) ;

if T~=D
    error('input y_sim and y_obs must have the same number of columns')
end
if M>1
    error('input y_obs must be a row vector')
end

Err  = y_sim - repmat(y_obs,N,1)   ;
Err0 = y_obs - mean(y_obs) ;
nse = 1 - sum(Err.^2,2) / sum(Err0.^2) ;
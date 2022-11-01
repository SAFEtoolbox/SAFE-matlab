function Y = tv_objfun(TSsim,TSobs,obj_fun,W)
%
% This function returns a time-varying metric of model performance:
%
% Y(t) = PERF[ sim(t-w, ..., t, ..., t+w) ,
%              obs(t-w, ..., t, ..., t+w)  ]     for t=1,...,T
%
% Usage:
% Y =  tv_objfun(TSsim,TSobs,obj_fun,W)
% 
% Input:
%  TSsim = set of N time series of simulated output          - matrix (N,T)
%  TSobs = set of N time series of observed output           - matrix (N,T)
%          (or 1 time series if all rows in TSsim           or vector (1,T)
%           should be evaluated against the same
%           time series of observations)
% obj_fun = performance metric to be computed                - string
%           options: 'MAE','RMSE','BIAS' (*)
% W = semi-length of the window for averaging performances   - scalar
%
% Output:
%      Y = set of N time series of time-varying objective    - matrix (N,T)
%          function
%
% (*) Comment: 
% Performance metrics are defined as:
% MAE(t)  = mean( | sim(t-w,...,t+w) - obs(t-w,..., t+w) | )
% RMSE(t) = mean( ( sim(t-w,...,t+w) - obs(t-w,..., t+w) )^2 )
% BIAS(t) = | mean(sim(t-w,...,t+w)) - mean(obs(t-w,..., t+w)) |

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

if ~isnumeric(TSsim) ; error('input ''TSsim'' must be a matrix of real numbers'); end
if ~isnumeric(TSobs) ; error('input ''TSobs'' must be a matrix of real numbers'); end
[N,T] = size(TSsim) ;
[n,t] = size(TSobs) ;
if t~=T
    error('inputs ''TSsim'' and ''TSobs'' should have the same number of columns');
end
if n~=N
    if n==1
        TSobs = repmat(TSobs,N,1);
    else
        error('inputs ''TSsim'' and ''TSobs'' should have the same number of rows');
    end
end

if ~ischar(obj_fun) ; error('input ''obj_fun'' must be a string'); end

if ~isnumeric(W) ; error('input ''W'' must be a scalar number'); end
if ~isscalar(W); error('''W'' must be scalar'); end
if W<1; error('''W'' must be >= 1' ); end
if abs(W-round(W)); error('''W'' must be an integer'); end

Err = TSsim - TSobs ; % (N,T)
Y   = nan(N,T)      ;

for t=1:T    
    if strcmp(obj_fun,'MAE')
        
        err = Err(:,max(1,t-W):min(T,t+W)) ; % (N,W+1)
        Y(:,t) = mean(abs(err),2) ;          % (N,1)
    
    elseif strcmp(obj_fun,'RMSE')
        
        err = Err(:,max(1,t-W):min(T,t+W)) ; % (N,W+1)
        Y(:,t) = sqrt(mean(err.^2,2)) ;      % (N,1)
    
    elseif strcmp(obj_fun,'BIAS')
        
        qobs = TSobs(:,max(1,t-W):min(T,t+W)) ; % (N,W+1)
        qsim = TSsim(:,max(1,t-W):min(T,t+W)) ; % (N,W+1)       
        Y(:,t) = abs( mean(qsim,2) - mean(qobs,2) ) ; % (N,1)
        
    else
        fprintf('\n ERROR! obj_fun %s not available\n',obj_fun)
    end
    
end

function [ Y, varargout ] = model_execution(fun_test,X,varargin)
%
% Evaluates the model coded in the matlab function 'fun_test' against
% each sample input vector in matrix X, and returns the associated output
% vectors (in matrix Y).
%
% Usage:
%
% Y = model_execution(fun_test,X)
%
% fun_test = name of the function implementing the model (string)
%        X = matrix (NxM) of N input samples
%           (each input sample has M components)
%
%        Y = vector (NxP) of associated output samples
%            (P being the number of model outputs associated to
%            each sampled input combination)
%
% NOTE: If the 'fun_test' function requires other arguments besides 'X',
% or produce other output besides 'Y', they can be passed as optional
% arguments after 'X' and recovered as optional output after 'Y'.
% Example:
%
% fun_test = 'sobol_g_function' ; 
% X = rand(3,4) ;
% a = ones(1,4) ;
% Y = model_execution(fun_test,X,a) ; % Or:
% [Y,V] = model_execution(fun_test,X,a) ; 

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

if ~ischar(fun_test); error('''fun_test'' must be a string'); end
if exist(fun_test,'file')~=2; error('function ''%s'' does not exist',fun_test); end

%%%%%%%%%%%%%%%%%%
% Model evaluation
%%%%%%%%%%%%%%%%%%

NumExtraArgOut = max(nargout,1) - 1;
NumExtraArgIn  = length(varargin)  ;

if nargin > 2   
    tmp = feval(fun_test,X(1,:),varargin{1:NumExtraArgIn}) ;
else
    tmp = feval(fun_test,X(1,:)) ;
end
P = length(tmp) ; % number of model output

N = size(X,1);
Y = nan(N,P) ;

if (NumExtraArgIn>0)&&(NumExtraArgOut>0)
    for j=1:N
        [Y(j,:),varargout{1:NumExtraArgOut}] = feval(fun_test,X(j,:),varargin{1:NumExtraArgIn}) ;
    end
elseif (NumExtraArgIn>0)&&(NumExtraArgOut<=0)
    for j=1:N
        Y(j,:) = feval(fun_test,X(j,:),varargin{1:NumExtraArgIn}) ;
    end
elseif (NumExtraArgIn<=0)&&(NumExtraArgOut>0)
    for j=1:N
        [Y(j,:),varargout{1:NumExtraArgOut}] = feval(fun_test,X(j,:)) ;
    end
else
    for j=1:N
        Y(j,:) = feval(fun_test,X(j,:)) ;
    end
end



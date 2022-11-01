function [YY, comp_time] = pawn_model_execution(fun_test,XX,varargin)
%
% This function evaluate a given model against the random samples generated
% by the PAWN_SAMPLING function. 
%
% Usage:
%               YY = pawn_model_execution(fun_test,XX)
% [YY ,compt_time] = pawn_model_execution(fun_test,XX)
%
% Input:
% XX = cell array of size (M,n). The element X{i,k} is the subsamples
%      of input datapoints to compute the sensitivity of factor i in
%      subsample k, and it is a matrix of size (NC,M).
%      XX can be generated using the PAWN_SAMPLING function
%      (see help pawn_sampling to learn more).
% fun_test = name of the function implementing the model (string)

%
% Output:
% YY = cell array of size (M,n).  The element X{i,k} is the subsamples
%      of output values to compute the sensitivity of factor i in
%      subsample k, and it is a vector of size (NC,P)
% comp_time = total computing time for model evaluation (sec)
%
% NOTE: If the 'fun_test' function requires other arguments besides 'XX',
% they can be passed as optional arguments after 'X' (see function 
% model_evaluation).

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

if nargin > 2
    NumExtraArg = length(varargin);
end

[M,n] = size(XX) ;
YY    = cell(M,n);
tic
for i=1:M
    for k=1:n
        YY{i,k} = model_evaluation(fun_test,XX{i,k},varargin{1:NumExtraArg});
    end
end
comp_time=toc;
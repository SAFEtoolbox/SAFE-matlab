function [Y_noNaN,YY_noNaN]=pawn_exclude_NaN(Y,YY,varargin)

% Exclude NaN values in output samples generated to estimate unconditional
% and conditional CDFs.
% 
% Usage:
% [Y_noNaN,YY_noNaN]=pawn_exclude_nan(Y,YY)
% [Y_noNaN,YY_noNaN]=pawn_exclude_nan(Y,YY,idx)
%
% Input:
% Y = output sample for estimating the unconditional CDF    - matrix (NU,P)
% YY = output samples for the conditional CDFs           - cell array (M,n) 
%       The element Y{i,k} is the k-th subsample 
%       of output evaluations obtained by fixing the i-th input 
%       (or group of inputs) to its k-th conditioning value,
%       and it is a matrix of size (NC,P).
% idx = when dealing with multiple output (P>1),                   - scalar
%       this function will actually estimate the conditional and
%       unconditional CDFs of one output only. 'idx' thus is the
%       index of the output to be analyzed (i.e. the column of
%       matrix Y and YY{i,k})                  (default: 1)
%
% Output
% Y_noNaN = output sample for estimating the            - matrix (NU_new,P)
%           unconditional  CDF where NaN values 
%           were removed
% YY = output samples for the conditional CDFs           - cell array (M,n) 
%      all NaN values were removed for each 
%      element Y{i,k}
% 
% NOTE: This function is called in pawn_indices and pawn_cdfs so that NaNs
% are excluded from the calculation of the indices and CDFs.
%
% REFERENCES
%
% Pianosi, F. and Wagener, T. (2015), A simple and efficient method 
% for global sensitivity analysis based on cumulative distribution 
% functions, Env. Mod. & Soft., 67, 1-11.

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

if ~isnumeric(Y); error('''Y'' must be a matrix'); end
[ ~, P ] = size(Y) ;
if ~iscell(YY); error('''YY'' must be a cell array'); end
[ M, n ] = size(YY) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
idx=1;

% Recover and update optional arguments:
if nargin > 2
    if ~isempty(varargin{1})
        idx = varargin{1} ;
        if ~isscalar(idx) ; error('input ''idx'' must be scalar'); end
        if abs(idx-round(idx)); error('''idx'' must be integer'); end
        if idx<1  ; error('input ''idx'' must be a positive integer'); end
        if idx>P  ; error('''idx'' exceeds the number of columns in Y'); end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Exclude NaN value in Y
%%%%%%%%%%%%%%%%%%%%%%%%%
if all(isnan(Y(:,idx))); error('all data in ''Y'' are NaN'); end
nan_Y=isnan(Y(:,idx)); % find any NaN in Y
idx_Y=sum(nan_Y); % total number of NaN in Y
Y_noNaN=Y(~nan_Y,idx);

% Print warning message if NaN values are found
if idx_Y
    fprintf('\n WARNING: %d NaN were found in Y',idx_Y)
    fprintf('\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%
% Exclude NaN value in YY
%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate variable
YY_noNaN=cell(size(YY));

for i=1:M % for each input factor   
    for k=1:n %for each conditioning value
        
        nan_YYik=isnan(YY{i,k}(:,idx)); % find any NaN in YY{i,k}
        idx_YYik=sum(nan_YYik); % total number of NaN in YY{i,k}
        YY_noNaN{i,k}=YY{i,k}(~nan_YYik,idx);
        
        % Print warning message if NaN values are found
        if idx_YYik
            fprintf('\n WARNING: %d NaN were found in YY{%d,%d}',idx_YYik,i,k)
            fprintf('\n')
            
        end
    end
end


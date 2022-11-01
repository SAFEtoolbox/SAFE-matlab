function KS = pawn_ks(YF,FU,FC,varargin)
%
% Compute the Kolmogorov-Smirnov (KS) statistic between the empirical
% unconditional CDF and the conditional CDFs at different conditioning
% values:
%
% KS(k,i) = max | F(y) - F(y|x(i)=x_ik) |
%            y
%
% i=1,...,M (number of inputs)
% k=1,...,n (number of conditioning values for each input)
%
% Usage:
% KS = pawn_ks(YF,FU,FC)
%
% Input:
%  YF = values of y at which the CDFs are given              - vector (P,1)
%  FU = values of the empirical unconditional CDF F(y)       - vector (P,1)
%  FC = cell array whose element (i,k) is a              - cell array (M,n)
%       vector (P,1) of values of the empirical
%       conditional CDF F(y|x(i)=x_ik)
%
% Output:
% KS = matrix (n,M) of KS statistic
%
% Example:
%
% NU = 150 ;
% NC = 100 ;
% n  = 10  ;
% M  = 3   ;
% XU = AAT_sampling('lhs',M,'unif',[-pi,pi],NU);
% YU = model_evaluation('ishigami_homma_function',XU) ;
% XX = pawn_sampling('lhs',M,'unif',[-pi,pi],n,NC);
% YY = pawn_model_evaluation('ishigami_homma_function',XX) ;
% [ YF, FU, FC  ] = pawn_cdfs(YU,YY) ;
% KS = pawn_ks(YF,FU,FC);
%
% -------------------------------------------------------------------------
% ADVANCED USAGE
% for Regional-Response Global Sensitivity Analysis:
% -------------------------------------------------------------------------
%
% KS = pawn_ks(YF,FU,FC,'above',Ythreshold)
%
% Compute KS considering only output values satisfying the condition:
% Y >= Ythreshold
% Use 'below' instead of 'above' for the condition Y<=Ythreshold.
% For more sophisticate conditions, the user can define its own
% function and use it as follows:
%
% KS = pawn_ks(YF,FU,FC,'output_condition',param)
%
% Where 'output_condition' is a function with the following structure:
%
% idx = output_condition(Y,param)
%
% where   Y = vector (N,1) of output samples
%     param = vector (any size) of parameters to check the condition
%       idx = vector (N,1) of logical values (1 if condition is satisfied,
%                                             0 otherwise)
%
% NOTE: The function discard empty conditional outputs (the function
% pawn_cdfs returns an empty FC{i,k} when all values are NaNs in the
% corresponding sample YY{i,k}.

% REFERENCES:
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

if ~isnumeric(YF); error('input ''YF'' must be numeric'); end
[P,m]=size(YF) ;
if m~=1; error('input ''YF'' must be a column vector'); end

if ~isnumeric(FU); error('input ''FU'' must be numeric'); end
[PU,m]=size(FU) ;
if m~=1; error('input ''FU'' must be a column vector'); end
if P~=PU; error('input ''UF'' and ''FU'' must have the same number of rows'); end

if ~iscell(FC); error('input ''FC'' must be a cell array'); end
[M,n]  = size(FC)  ;
for i=1:M;
    for k=1:n;
        if ~isempty(FC{i,k})
            if ~isnumeric(FC{i,k}); error('element (%d,%d) of ''FC'' is not numeric',i,k); end
            [PC,m]=size(FC{i,k}) ;
            if m~=1; error('element (%d,%d) of ''FC'' must be a column vector',i,k); end
            if P~=PC; error('element (%d,%d) of ''FC'' must have the same number of rows as ''YF''',i,k); end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
output_condition = 'allrange' ;
par = NaN ;
% Recover and update optional arguments:
if nargin > 3
    if ~isempty(varargin{1})
        output_condition = varargin{1} ;
        if ~ischar(output_condition); error('''output_condition'' must be a string'); end
    end
    if isempty(varargin{2})
        error('please specify the parameters for the output condition')
    else
        par = varargin{2} ;
    end
end

KS = NaN(n,M)  ;

% Find subset of output values specifying a given condition
idx = feval(output_condition,YF,par);

for i=1:M % for each input:
    
    for k=1:n % for each conditioning value:
        if ~isempty(FC{i,k} )
            FCik = FC{i,k} ;
            % Compute KS:
            KS(k,i)=max(abs(FU(idx)-FCik(idx))) ;
        end
    end
    
end

% % Print results to screen:
% str='X%d:'; for k=1:n; str = [ str ' \t %2.3f' ]; end; str = [ str ' \n'];
% fprintf(str,[ 1:M; KS ]);

% Built-in functions:

function idx = above(y,par)

idx = y >= par ;

function idx = below(y,par)

idx = y <= par ;

function idx = allrange(y,par)

idx = true(size(y)) ;
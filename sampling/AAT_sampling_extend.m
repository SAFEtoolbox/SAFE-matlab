function Xext = AAT_sampling_extend(X,distr_fun,distr_par,Next,varargin)
%
% This function create an expanded sample 'Xext' starting from
% a sample 'X' and using latin hypercube and the maximin criterion.
% 
% Usage:
% Xext = AAT_sampling_extend(X,distr_fun,distr_par,Next)
% Xext = AAT_sampling_extend(X,distr_fun,distr_par,Next,N_rep)
% 
% Input:
%          X = initial sample                                 - matrix(N,M)
%  distr_fun = probability distribution function of each variable 
%                  - string (eg: 'unif') if all variables have the same pdf
%                  - cell array of M strings (eg:{'unif','norm'}) otherwise
%  distr_par = parameters of the probability distribution function
%                         - row vector if all input variables have the same
%                         - cell array of M vectors otherwise
%       Next = new dimension of the sample (must be >N)            - scalar
%       nrep = number of replicate to select the maximin hypercube - scalar
%             (default value: 10)
% 
% Output:
%       Xext = expanded sample                             - matrix(Next,M)
% 
% See also AAT_sampling.
% 
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

if ~isnumeric(X) ; error('input ''X'' must be a matrix of size (N,M)'); end
[N,M]=size(X) ;

if ischar(distr_fun)
    tmp=cell(1,M);for i=1:M; tmp{i}=distr_fun; end; distr_fun = tmp;
elseif iscell(distr_fun)
    if length(distr_fun)~=M
        error('If ''distr_fun'' is a cell array, it must have M components')
    end
else
    error('Wrong data type for input ''distr_fun''')
end
if isnumeric(distr_par)
    [n1,n2] = size(distr_par) ;
    if n1==1
        tmp=cell(1,M);for i=1:M; tmp{i}=distr_par; end; distr_par = tmp;
    else
        error('If ''distr_par'' is a vector, it must be a row vector')
    end
elseif iscell(distr_par)
    if length(distr_par)~=M
        error('If ''distr_par'' is a cell array, it must have M components')
    end
else
    error('Wrong data type for input ''distr_par''')
end

if ~isscalar(Next); error('''Next'' must be a scalar'); end
if (Next - round(Next))~=0 ; error('''Next'' must be an integer'); end
if Next<=0 ; error('''Next'' must be a positive integer'); end
if Next < N ; error('''Next'' must be larger than ''N''');end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
nrep=10;

% Recover and update optional arguments:
if nargin > 4
    if ~isempty(varargin{1})
        nrep = varargin{1};
        if ~isscalar(nrep); error('''nrep'' must be a scalar'); end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map back the original sample into
% the uniform unit square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = nan(size(X)) ;

for i=1:M
    pars = distr_par{i} ;
    np   = length(pars) ;
    name = distr_fun{i} ;
    if any(strcmpi(name,{'chi2','exp','geo','poiss','rayl','t','unid'}))
        if np~=1; error('Input %d: Number of PDF parameters not consistent with PDF type',i);
        else
            U(:,i) = feval([name 'cdf'],X(:,i),pars);
        end
        
    elseif any(strcmpi(name,{'beta','bino','ev','f','gam','logn','nbin','nct','ncx2','norm','unif','wbl'}))
        if np~=2; error('Input %d: Number of PDF parameters not consistent with PDF type',i);
        else
            U(:,i) = feval([name 'cdf'],X(:,i),pars(1),pars(2));
        end
        
    elseif any(strcmpi(name,{'gev','gp','hyge','ncf'}))
        if np~=3; error('Input %d: Number of PDF parameters not consistent with PDF type',i);
        else
            U(:,i) = feval([name 'cdf'],X(:,i),pars(1),pars(2),pars(3));
        end
    else
        error('Input %d: Unknown PDF type',i)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add samples in the unit square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uext = lhcube_extend(U,Next,nrep) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map back into the specified distribution
% by inverting the CDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xext = nan(size(Uext)) ;

for i=1:M
    pars = distr_par{i};
    np   = length(pars) ;
    name = distr_fun{i}    ;
    if any(strcmpi(name,{'chi2','exp','geo','poiss','rayl','t','unid'}))
        if np~=1; error('Input %d: Number of PDF parameters not consistent with PDF type',i);
        else
            Xext(:,i) = feval([name 'inv'],Uext(:,i),pars);
        end
        
    elseif any(strcmpi(name,{'beta','bino','ev','f','gam','logn','nbin','nct','ncx2','norm','unif','wbl'}))
        if np~=2; error('Input %d: Number of PDF parameters not consistent with PDF type',i);
        else
            Xext(:,i) = feval([name 'inv'],Uext(:,i),pars(1),pars(2));
        end
        
    elseif any(strcmpi(name,{'gev','gp','hyge','ncf'}))
        if np~=3; error('Input %d: Number of PDF parameters not consistent with PDF type',i);
        else
            Xext(:,i) = feval([name 'inv'],Uext(:,i),pars(1),pars(2),pars(3));
        end
    else
        error('Input %d: Unknown PDF type',i)
    end
end



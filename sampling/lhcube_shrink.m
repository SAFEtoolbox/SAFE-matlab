function [X_new,idx_new]=lhcube_shrink(X,N_new,varargin)
%
% This function drop rows from a latin hypercube using the maximin 
% criterion.
% 
% Usage:
% [X_new,idx_new]=lhcube_shrink(X,N_new)
% [X_new,idx_new]=lhcube_shrink(X,N_new,nrep)
% 
% 
% Input:
%       X = initial latin hypercube                           - matrix(N,M)
%   N_new = new dimension for the latin hypercube             - scalar
%    nrep = number of replicate to select the best hypercube  - scalar
%          (default value: 10)
% 
% Output:
%   X_new = new latin hypercube                          - matrix(N_new,M)
% idx_new = indices of the rows selected from X          - vector(N_new,1)
%            [ i.e. Xnew = X(idx_new,:) ]
% 
% Example:
% 
% N = 30 ;
% M =  2 ;
% X = lhcube(N,M); % create LHS
% figure(1)
% plot(X(:,1),X(:,2),'x')
% set(gca,'XTick',[0:1/N:1],'XTickLabel',{},'YTick',[0:1/N:1],'YTickLabel',{})
% grid on
% 
% N_new = 20 ;
% X_new = lhcube_shrink(X,N_new);
% figure(1); hold on
% plot(X_new(:,1),X_new(:,2),'or')
% 
% References:
% 
% this is the matlab/octave version of the code in:
%         J.C. Rougier, calibrate package 
%         http://www.maths.bris.ac.uk/~mazjcr/#software
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
if ~isscalar(N_new); error('''N_new'' must be a scalar'); end
if (N_new - round(N_new))~=0 ; error('''N_new'' must be an integer'); end
if N_new<=0 ; error('''N_new'' must be a positive integer'); end
if N_new>N  ; error('''N_new'' must be smaller than ''N''');end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
nrep=10;

% Recover and update optional arguments:
if nargin > 2
    if ~isempty(varargin{1})
        nrep = varargin{1};
        if ~isscalar(nrep); error('''nrep'' must be a scalar'); end
    end
end

%%%%%%%%%%%%
% Drop rows
%%%%%%%%%%%%

X_new   = nan ;
ddbest  = 0   ;
idx_new = nan ;

% Generate nrep subsamples of the initial latin hypercube sample X and keep 
% the one that maximises the minimum inter-point Euclidean distance between 
% any two sampled points:
for i=1:nrep 
    idx=randperm(N,N_new);
    Xi=X(idx,:);
    dd=min(pdist(Xi)); % use this if you have the Matlab Statistical Toolbox
    % dd=(ipdm(Xi,'metric',2)); dd=min(dd(dd>0));
    if dd>ddbest
        X_new=Xi;
        idx_new=idx;
        ddbest=dd;
    end
end
        
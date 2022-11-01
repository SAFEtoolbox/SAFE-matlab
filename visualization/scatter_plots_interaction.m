function [] = scatter_plots_interaction(X,Y,varargin)
%
% This function produces scatter plots of X(i) against X(j),
% for all possible combinations of (i,j). In each plot,
% the marker colour is proportional to the value of Y.
%
% Usage: 
% scatter_plots_interaction(X,Y)
% scatter_plots_interaction(X,Y,ms)
% scatter_plots_interaction(X,Y,ms,X_Labels)
%
% Input:
%        X = set of input samples                      - matrix (N,M)
%        Y = set of output samples                     - vector (N,1)
%       ms = marker size (default: 10)                 - scalar
% X_Labels = cell array of model inputs names          - cell array (1,M)
%            (default: {'#1','#2',...,'#M'})
%
% Example:
%
% X = rand(100,3) ; 
% Y = sin(X(:,1)) + 2 * sin(X(:,2)).^2  + X(:,3).^4.*sin(X(:,1));
% scatter_plots_interaction(X,Y)
% scatter_plots_interaction(X,Y,16)
% scatter_plots_interaction(X,Y,[],{'x1','x2','x3'})

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isnumeric(X); error('''X'' must be a matrix of size (N,M)'); end
if ~isnumeric(Y); error('''Y'' must be a vector of size (N,1)'); end
[N,M]=size(X) ;
[n,m]=size(Y) ;
if N~=n; error('input ''X'' and ''Y'' must have the same number of rows'); end
if m~=1; error('input ''Y'' must be a column vector'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
ms = 10;
X_Labels = cell(1,M); for i=1:M; X_Labels{i}=['#' num2str(i)]; end

% Recover and update optional arguments:
if nargin > 2
    if ~isempty(varargin{1})
        ms = varargin{1}; % marker size
        if ~isscalar(ms); error('''ms'' must be a scalar'); end
        if ms<=0; error('''ms'' must be positive' ); end
        if abs(ms-round(ms)); error('''ms'' must be integer'); end
    end
end
if nargin > 3
    if ~isempty(varargin{2})
        X_Labels = varargin{2};
        if ~iscell(X_Labels); error('''X_Labels'' must be a cell array'); end
        if length(X_Labels)~=M; error('''X_Labels'' must have M=%d components',M); end
        for i=1:M; if ~ischar(X_Labels{i}); error('all components of ''X_Labels'' must be string'); end; end
    end
end

figure
k=1;
for i=1:M-1
    for j=i+1:M
        subplot(M-1,M-1,k)
        scatter(X(:,i),X(:,j),ms,Y,'fill'); set(gca,'XTick',[],'YTick',[])
        axis([min(X(:,i)),max(X(:,i)),min(X(:,j)),max(X(:,j))])
        title(['( ' X_Labels{i} ' vs ' X_Labels{j} ' )'])
        k=k+1;
    end
    k=k+i;
end


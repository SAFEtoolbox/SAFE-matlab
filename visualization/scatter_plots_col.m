function scatter_plots_col(X,Y,i1,i2,varargin)
%
% This function produces scatter plots of X(i1) against X(i2).
% The marker colour is proportional to the value of Y.
%
% Usage: 
% scatter_plots_col(X,Y,i1,i2)
% scatter_plots_col(X,Y,i1,i2,marker_size)
% scatter_plots_col(X,Y,i1,i2,marker_size,X_Labels)
%
% Input:
%           X = set of input samples                   - matrix (N,M)
%           Y = set of output samples                  - vector (N,1)
%          i1 = index of input on the horizontal axis  - scalar
%          i2 = index of input on the vertical axis    - scalar
% marker_size = size of marker (default: 7)            - scalar
%    X_Labels = cell array of model inputs names       - cell array (1,M)
%               (default: {'#1','#2',...,'#M'})
%
% Example:
%
% X = rand(100,3) ; 
% Y = sin(X(:,1)) + 2 * sin(X(:,2)).^2  + X(:,3).^4.*sin(X(:,1));
% scatter_plots_col(X,Y,1,2)
% scatter_plots_col(X,Y,1,2,16)
% scatter_plots_col(X,Y,1,2,[],{'x1','x2','x3'})

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

% Options for the graphic:
fn = 'Helvetica' ; % font type of axes, labels, etc.
%fn = 'Courier' ;
fs = 20 ; % font size of axes, labels, etc.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isnumeric(X); error('''X'' must be a matrix of size (N,M)'); end
if ~isnumeric(Y); error('''Y'' must be a vector of size (N,1)'); end
[N,M]=size(X) ;
[n,m]=size(Y) ;
if N~=n; error('input ''X'' and ''Y'' must have the same number of rows'); end
if m~=1; error('input ''Y'' must be a column vector'); end

if ~isscalar(i1); error('''i1'' must be a scalar'); end
if i1<=0; error('''i1'' must be positive' ); end
if abs(i1-round(i1)); error('''i1'' must be integer'); end
if i1>M; error('input indes %d exceeds the number of columns in X',i1); end

if ~isscalar(i2); error('''i2'' must be a scalar'); end
if i2<=0; error('''i2'' must be positive' ); end
if abs(i2-round(i2)); error('''i2'' must be integer'); end
if i2>M; error('input indes %d exceeds the number of columns in X',i2); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
marker_size=7;
X_Labels = cell(1,M); for i=1:M; X_Labels{i}=['#' num2str(i)]; end

% Recover and update optional arguments:
if nargin > 4
    if ~isempty(varargin{1})
        marker_size = varargin{1};
        if ~isscalar(marker_size); error('''marker_size'' must be a scalar'); end
        if marker_size<=0; error('''marker_size'' must be positive' ); end
        if abs(marker_size-round(marker_size)); error('''marker_size'' must be integer'); end
    end
end

if nargin > 5
    if ~isempty(varargin{2})
        X_Labels = varargin{2};
        if ~iscell(X_Labels); error('''X_Labels'' must be a cell array'); end
        if length(X_Labels)~=M; error('''X_Labels'' must have M=%d components',M); end
        for i=1:M; if ~ischar(X_Labels{i}); error('all components of ''X_Labels'' must be string'); end; end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%figure
scatter(X(:,i1),X(:,i2),marker_size,Y,'fill'); box on
xlabel(X_Labels{i1},'FontSize',fs,'FontName',fn)
ylabel(X_Labels{i2},'FontSize',fs,'FontName',fn)
colorbar
axis([min(X(:,i1)),max(X(:,i1)),min(X(:,i2)),max(X(:,i2))])
set(gca,'FontSize',fs,'FontName',fn)


function scatter_plots(X,Y,varargin)
%
% This function produces M scatter plots, each one plotting the output
% sample (Y) against one component of the input vector sample
% (i.e. one column of X).
%
% Usage: 
% scatter_plots(X,Y)
% scatter_plots(X,Y,n_col)
% scatter_plots(X,Y,n_col,Y_Label)
% scatter_plots(X,Y,n_col,Y_Label,X_Labels)
% scatter_plots(X,Y,n_col,Y_Label,X_Labels,idx)
%
% Input:
%           X = set of model inputs samples             - matrix (N,M)
%           Y = set of model output samples             - vector (N,1)
%       n_col = number of panels per row                - scalar
%               (default: min(5,M))
%     Y_Label = label of vertical axis                  - string
%               (default: 'Y')
%    X_Labels = cell array of model inputs names        - cell array (1,M)
%               (default: {'X1','X2',...,'XM'})
%         idx = indices of datapoints to be highlighted - vector (N,1)
%
% Example:
%
% X = rand(100,3) ; 
% Y = sin(X(:,1)) + 2 * sin(X(:,2)).^2  + X(:,3).^4.*sin(X(:,1));
% scatter_plots(X,Y)
% scatter_plots(X,Y,2)
% scatter_plots(X,Y,[],'y')
% scatter_plots(X,Y,[],'y',{'x1','x2','x3'})
% scatter_plots(X,Y,[],'y',{'x1','x2','x3'},Y>2)

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
mt  = 'o' ; % marker type
mth = 'o' ; % marker type for highlighted datapoints

% Options for the colours:
% You can produce a coloured plot or a black and white one
% (printer-friendly).
% Option 1 - coloured plot: uncomment the following 4 lines:
me  = 'b' ; % marker edge colour
meh = 'r' ; % marker edge colour for highlighted datapoints
mc  = 'b' ; % marker face colour
mch = 'r' ; % marker face colour for highlighted datapoints
% Option 2 - B&W plot: uncomment the following 4 lines:
% me  = [100 100 100]/256 ; % marker edge colour
% meh = 'k' ; % marker edge colour for highlighted datapoints
% mc  = 'w' ; % marker face colour
% mch = 'k' ; % marker face colour for highlighted datapoints

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
n_col    = 5 ;
Y_Label  = 'Y' ;
X_Labels = cell(1,M); for i=1:M; X_Labels{i}=['X' num2str(i)]; end
idx = [];

% Recover and update optional arguments:
if nargin > 2
    if ~isempty(varargin{1})
        n_col = varargin{1}; % max number of panels per row
        if ~isscalar(n_col); error('''n_col'' must be a scalar'); end
        if n_col<=0; error('''n_col'' must be positive' ); end
        if abs(n_col-round(n_col)); error('''n_col'' must be integer'); end
    end
end

if nargin > 3
    if ~isempty(varargin{2})
        Y_Label = varargin{2};
        if ~ischar(Y_Label); error('''Y_Label'' must be a string'); end
    end
end

if nargin > 4
    if ~isempty(varargin{3})
        X_Labels = varargin{3};
        if ~iscell(X_Labels); error('''X_Labels'' must be a cell array'); end
        if length(X_Labels)~=M; error('''X_Labels'' must have M=%d components',M); end
        for i=1:M; if ~ischar(X_Labels{i}); error('all components of ''X_Labels'' must be string'); end; end
    end
end

if nargin > 5
    if ~isempty(varargin{4})
        idx = varargin{4}; 
        [n,m]=size(idx) ;
        if N~=n; error('input ''idx'' must have the same number of rows as ''X'' and ''Y'' '); end
        if m~=1; error('input ''idx'' must be a column vector'); end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_col = min(n_col,M)  ;
n_row = ceil(M/n_col) ;
%figure
for i=1:M
    
    subplot(n_row,n_col,i)
    plot( X(:,i), Y, mt , 'MarkerFaceColor', mc, 'MarkerEdgeColor', me )
    xlabel(X_Labels{i},'FontSize',fs,'FontName',fn)
    
    % add axis, labels, etc.
    if mod(i,n_col)==1; ylabel(Y_Label,'FontSize',fs,'FontName',fn); end
    axis([min(X(:,i)),max(X(:,i)),min(Y)-std(Y),max(Y)+std(Y)])
    set(gca,'FontSize',fs,'FontName',fn)
    
    if ~isempty(idx);
       hold on; plot( X(idx,i), Y(idx), mth, 'MarkerFaceColor', mch, 'MarkerEdgeColor', meh )
    end

end


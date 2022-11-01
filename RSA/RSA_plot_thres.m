function [] = RSA_plot_thres(X,idxb,varargin)
%
% Plotting function for Regional Sensitivity Analysis.
% Plot the CDF of the samples in dataset X that satisfy a given condition
% ('behavioural'), and the CDF of the samples that do not satisfy 
% the condition ('non-bahvioural').
% (see help of RSA_indices_thres for details about RSA and references)
%
% Usage: 
%           RSA_plot_thres(X,idxb) 
%           RSA_plot_thres(X,idxb,n_col)  
%           RSA_plot_thres(X,idxb,n_col,labels) 
%           RSA_plot_thres(X,idxb,n_col,labels,str_legend)
%
% Input:
%         X = set of input samples                           - matrix (N,M)
%      idxb = indices of samples statisfying the condition   - vector (N,1)
%     n_col = number of panels per row in the plot                 - scalar
%             (default: min(5,M))
%    labels = labels for the horizontal axis             - cell array (1,M)
%             (default: {' X1','X2',...})
% str_legend = text for legend                           - cell array (1,2)

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

lwn = 2 ; % line width for non-behavioural records
lwb = 2 ; % line width for behavioural records

% Options for the colours:
% You can produce a coloured plot or a black and white one
% (printer-friendly).
% Option 1 - coloured plot: uncomment the following 2 lines:
lcn = 'b' ; % line colour for non-behavioural
lcb = 'r' ; % line colour for behavioural
% Option 2 - B&W plot: uncomment the following 2 lines:
%lcn =  [150 150 150]/256 ; % line colour
%lcb = 'k' ; % line colour for behavioural


%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(X) ; error('input ''X'' must be a matrix of size (N,M)'); end
if ~islogical(idxb) ; error('input ''idxb'' must be a vector of N logical'); end
[N,M]=size(X) ;
[n,m]=size(idxb) ;
if N~=n; error('input ''X'' and ''idxb'' must have the same number of rows'); end
if m~=1; error('input ''idxb'' must be a column vector'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
n_col = 5 ;
X_Labels = cell(1,M); for i=1:M; X_Labels{i}=['X' num2str(i)]; end
str_legend={'behavioural','non-behavioural'};

% Recover and update optional arguments:
if nargin > 2
    if ~isempty(varargin{1})
        n_col = varargin{1};
        if ~isscalar(n_col); error('''n_col'' must be a scalar'); end
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
if nargin > 4
    if ~isempty(varargin{3})
        str_legend = varargin{3};
        if ~iscell(str_legend); error('''str_legend'' must be a cell array'); end
        if length(str_legend)~=2; error('''str_legend'' must have 2 components'); end
        for i=1:2; if ~ischar(str_legend{i}); error('all components of ''str_legend'' must be string'); end; end
    end
end


%%%%%%%%
% Plot
%%%%%%%%

n_col = min(floor(n_col),M)  ;
n_row = ceil(M/n_col) ;

Xb  = X(idxb,:) ;
Xnb = X(~idxb,:) ;

figure
for i=1:M       
    % Approximate behavioural and non-behavioural CDFs:
    xx = unique(sort(X(:,i))) ;
    CDFb  = empiricalcdf(Xb(:,i),xx)    ;
    CDFnb = empiricalcdf(Xnb(:,i),xx)   ;       
    % Plot CDFs:
    subplot(n_row,n_col,i); hold on
    plot(xx,CDFb,'Color',lcb,'LineWidth',lwb)
    plot(xx,CDFnb,'Color',lcn,'LineWidth',lwn)
    hold on
    xlabel(X_Labels{i},'FontName',fn,'FontSize',fs)
    set(gca,'FontName',fn,'FontSize',fs)
    if mod(i,n_col) == 1 % first column
        ylabel('cdf','FontName',fn,'FontSize',fs)
    end
    axis([min(X(1:N,i)),max(X(1:N,i)),0,1])
    if i==1;
        if ~isempty(str_legend); legend(str_legend); end
    end
    box on
    
 %   [ks,imax]=max(abs(CDFb-CDFnb)) ;
 %   hold on
 %   plot([xx(imax) xx(imax)],[CDFb(imax) CDFnb(imax)],'k')
    
end

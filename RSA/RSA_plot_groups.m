function []= RSA_plot_groups(X,idx,Yk,varargin)
%
% Plotting function for Regional Sensitivity Analysis with grouping.
% Plot 'Ng' CDFs of the samples in X with different colours.
% (see help of RSA_indices_groups for details about RSA and references)
%
% Usage: 
% RSA_plot_groups(X,idx,Yk)     
% RSA_plot_groups(X,idx,Yk,n_col) 
% RSA_plot_groups(X,idx,Yk,n_col,X_Labels,legend_title)
%
% Input:
%         X = set of input samples                           - matrix (N,M)
%       idx = index of group to which input samples belong   - vector (N,1) 
%        Yk = range of Y in each group                    - vector (Ng+1,1)   
%     n_col = number of panels per row (default: min(5,M))         - scalar
%  X_Labels = labels of horizontal axis                  - cell array (1,M)
%             (default: {'X1','X2',...})
% legend_title = label for legend (default: 'Y')            - string

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

fn = 'Helvetica' ; % font type of axes, labels, etc.
fs = 20 ; % font size of axes, labels, etc.
lw = 2 ; % line width
colorscale = 'jet' ; % use this for coloured plot
%colorscale = 'gray' ; % use this for black and white plot

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(X) ; error('input ''X'' must be a matrix of size (N,M)'); end
if ~isnumeric(Yk) ; error('input ''Yk'' must be a vector of size (1,nb_groups+1)'); end
if ~isnumeric(idx) ; error('input ''idx'' must be a vector of size (L,nb_groups)'); end
[N,M]=size(X) ;
[n,m]=size(Yk) ;
[N2,L]=size(idx) ;
if L~=1; error('input ''idx'' must be a vector of size (1,N)'); end
if N~=N2; error('input ''X'' and ''idx'' must have a consistent number of rows'); end
if m~=1; error('input ''Yk'' must be a column vector'); end
Ng=n-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:

n_col = 5 ;
X_Labels = cell(1,M); for i=1:M; X_Labels{i}=['X' num2str(i)]; end
legend_title = 'Y';

% Recover and update optional arguments:
if nargin > 3
    if ~isempty(varargin{1})
        n_col = varargin{1};
        if ~isscalar(n_col); error('''n_col'' must be a scalar'); end
    end
end
if nargin > 4
    if ~isempty(varargin{2})
        X_Labels = varargin{2};
        if ~iscell(X_Labels); error('''X_Labels'' must be a cell array'); end
        if length(X_Labels)~=M; error('''X_Labels'' must have M=%d components',M); end
        for i=1:M; if ~ischar(X_Labels{i}); error('all components of ''X_Labels'' must be string'); end; end
    end
end

if nargin > 5
    if ~isempty(varargin{2})
        legend_title = varargin{3};
        if ~ischar(legend_title); error('''legend_title'' must be string'); end; 
    end
end


%%%%%%
% RSA
%%%%%%

% set number of columns and rows in the plot:
n_col = min(floor(n_col),M)  ;
n_row = ceil(M/n_col) ;

figure

% set colour scale:
clrs = feval(colorscale,Ng);

for i=1:M % for each input

    xx  = unique(sort(X(:,i))) ;
    
    for j=1:Ng % for each group

        % Compute empirical CDF:
        CDFj = empiricalcdf(X(idx==j,i),xx) ;
        
        % Plot the CFD:
        subplot(n_row,n_col,i); hold on
        plot(xx,CDFj,'color',clrs(j,:),'LineWidth',lw)
    end 
    box on
    
    % Costumize axes:
    axis([min(X(1:N,i)),max(X(1:N,i)),0,1])
    set(gca,'FontName',fn,'FontSize',fs)
    
    % x-axis label:
    xlabel(X_Labels{i},'FontName',fn,'FontSize',fs)
    
    % y-axis label:
    if mod(i,n_col) == 1 % first column
        ylabel('cdf','FontName',fn,'FontSize',fs)
        %set(gca,'YTick',[0 1])
    else
        set(gca,'YTick',[])
    end
    
    
end

% Add colorbar:
axes('Position', [0.10 0.05 0.9 0.9], 'Visible', 'off');
for i=1:Ng+1; ctick{i}=sprintf('%1.2g',Yk(i)); end
colormap(feval(colorscale,Ng))
c=colorbar('YTickLabel',ctick,'FontName',fn,'Fontsize',fs);
ylabel(c,legend_title,'FontName',fn,'Fontsize',fs)

 
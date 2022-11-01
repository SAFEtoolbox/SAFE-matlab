function [] = pawn_plot_cdf(YF, FU, FC, varargin )
%
% Plot the unconditional output CDF (i.e. when all inputs vary)
% and the conditional CDFs (when one or more inputs are fixed).
%
% Usage:
% pawn_plot_cdf(YF, FU, FCi)
% pawn_plot_cdf(YF, FU, FCi, xci)
% pawn_plot_cdf(YF, FU, FCi, xci, y_label)
% pawn_plot_cdf(YF, FU, FCi, xci, y_label, x_label)
%
% Input:
%  YF = values of y at which the CDFs are given              - vector (P,1)
%  FU = values of the empirical unconditional CDF            - vector (P,1)
% FCi = cell array whose k-th element is a               - cell array (1,n)
%       vector (P,1) of values of the empirical 
%       conditional CDF at the k-th conditioning value
% xci = conditioning values for the input (used for legend)  - vector (n,1)
%       (default: not used, no legend)
% y_label = legend for the horizontal axis                         - string
%           (default: 'output y')
% x_label = legend for the colour bar                              - string
%           (default: empty)
%
% Example:
%
% NOTE: The function discard empty conditional outputs (the function
% pawn_cdfs returns an empty FC{i,k} when all values are NaNs in the
% corresponding sample YY{i,k}.
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
% figure
% for i=1:M
% subplot(1,M,i)
% pawn_plot_cdf(YF, FU, FC(i,:))
% end

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

fs = 20 ; % Fontsize

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
if M>1; error('input ''FC'' must be a cell array of size (1,n)'); end
for k=1:n
    if ~isempty(FC{k})
        if ~isnumeric(FC{k}); error('element %d of ''FC'' is not numeric',k); end
        [PC,m]=size(FC{k}) ;
        if m~=1; error('element %d of ''FC'' must be a column vector',k); end
        if P~=PC; error('element %d of ''FC'' must have the same number of rows as ''YF''',k); end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag = 0 ;
xc = 1:n ;
y_label = 'output y' ;
if nargin>3
    if ~isempty(varargin{1})
        xc = varargin{1} ;
        if ~isnumeric(xc) ; error('xc must be a column vector of %d elements',n); end
        [n1,m1] = size(xc) ;
        if m1~=1 ; error('xc must be a column vector'); end
        if n1~=n ; error('xc must be a column vector of %d elements',n); end
        flag    = 1  ; % colorbar with be added to the plot
        x_label = '' ;
    end
end
if nargin>4
    if ~isempty(varargin{2})
        y_label = varargin{2} ;
        if ~ischar(y_label); error('''y_label'' must be a string'); end
    end
end
if nargin>5
    if ~isempty(varargin{3})
        x_label = varargin{3} ;
        if ~ischar(x_label); error('''x_label'' must be a string'); end
        flag = 1 ;
    end
end

%%%%%%%%%%%%%%
% Create plot
%%%%%%%%%%%%%%

% Identify and discard empty conditional CDFs
idx_empty=false(1,n); %indices of the empty conditional CDFS
for i=1:n
    if isempty(FC{i})
        idx_empty(i)=true;
    end
end
FC=FC(~idx_empty); % remove empty CDFs
n=length(FC); % compute actual number of conditional CDFs
xc=xc(~idx_empty); % remove conditioning values when the corresponding CDF 
                   % is empty

% Prepare color and legend
col=gray(n+1);% Color for colormap, (n+1) so that the last line is not white
map = colormap(col(1:n,:));%Remove white color
[ccc,iii]=sort(xc);

% Plot a horizontal dashed line at F=1:
plot(YF,ones(size(YF)),'--k')
hold on

% Plot conditional CDFs in gray scale:
for j=1:n
    plot(YF,FC{iii(j)},'Color',map(j,:),'LineWidth',2)    
    set(gca,'FontSize',fs)%,'YLim',[0 1.1],'YTick',[0 0.5 1])
    xlabel(y_label,'FontSize',fs);
end
axis([min(YF),max(YF),-0.02,1.02])

if flag % Add colorbar on the left 
    ColorBarLabels = cell(1,n);
    for j=1:n; ColorBarLabels{j} = sprintf('%5.4g',ccc(j)); end
    h=colorbar('location','eastoutside','FontSize',12,...
        'YTick',1.5:1.5+n,'YTickLabel',ColorBarLabels);
    title(h,x_label,'FontSize',fs);
end

% Plot unconditional CDF in red:
plot(YF,FU,'r','LineWidth',3)

box on

%    mylegend=cell(1,ncluster+1);
%    for j=1:ncluster; mylegend{j}=['F(y|' x_labels{idx} '=' num2str(-ccc(j)) ')']; end
%    mylegend{ncluster+1}='F(y)';
%    legend(mylegend)
%    axis([Ym,YM,0,1.1])
%    if i==1; ylabel('CDF','FontSize',fs); end





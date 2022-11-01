function [] = parcoor(X,varargin)
%
% Create Parallel Coordinate Plot.
% 
% Usage:
%       parcoor(X)
%       parcoor(X,XLabel)
%       parcoor(X,XLabel,i_axis)
%       parcoor(X,XLabel,i_axis,idx)
% 
% Input:
%      X = set of samples                                    - matrix (N,M)
% XLabel = labels for the horizontal axis                - cell array (1,M)
%          (default: {' X1','X2',...})
% i_axis = index of input to be used for assigning units           - scalar
%          of measurement to the vertical axis if all the 
%          inputs have the same range.
%          (if empty or not specified the ranges of all the 
%          inputs are displayed)
%    idx = indices of samples statisfying some condition     - vector (N,1)
%          which will be highlighted in different color
% 
% Example:
% X = rand(100,3) ; 
% X(:,2)=X(:,2)+2 ; 
% parcoor(X,{'a','b','c'})
% parcoor(X,{'a','b','c'},2)
% idx = X(:,1)<0.3 ;
% parcoor(X,{'a','b','c'},2,idx)
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

% Options for the graphic:
fn = 'Helvetica' ; % font type of axes, labels, etc.
%fn = 'Courier' ;
fs = 20 ; % font size of axes, labels, etc.
% Text formating for parameter ranges
yticklabels_form = '%3.2f'; % float with 1 character after decimal point
% yticklabels_form = '%3.0f'; % integer
% lw  = 1 ; % line width
% lwh = 2 ; % line width for highlighted records
% 
% % Options for the colours:
% % You can produce a coloured plot or a black and white one
% % (printer-friendly).
% % Option 1 - coloured plot: uncomment the following 2 lines:
% lc  = 'b' ; % line colour
% lch = 'r' ; % line colour for highlighted records
% % Option 2 - B&W plot: uncomment the following 2 lines:
% %lc  =  [150 150 150]/256 ; % line colour
% %lch = 'k' ; % line colour for highlighted records

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(X) ; error('input ''X'' must be a matrix of size (N,M)'); end
[N,M] = size(X) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
XLabel = cell(1,M); for i=1:M; XLabel{i}=['X' num2str(i)]; end
i_axis = [] ;
idx = [] ;
if nargin > 1
    if ~isempty(varargin{1})
        XLabel = varargin{1};
        if ~iscell(XLabel); error('''XLabel'' must be a cell array'); end
        if length(XLabel)~=M; error('''XLabel'' must have M=%d components',M); end
        for i=1:M; if ~ischar(XLabel{i}); error('all components of ''XLabel'' must be string'); end; end
    end
end
if nargin > 2
    if ~isempty(varargin{2})
        i_axis = varargin{2};
        if ~isscalar(i_axis); error('''i_axis'' must be scalar'); end
        if i_axis<=0; error('''i_axis'' must be positive' ); end
        if abs(i_axis-round(i_axis)); error('''i_axis'' must be integer'); end
        if i_axis>M; error('''i_axis'' cannot be higher than M' ); end
    end
end
if nargin > 3
    if ~isempty(varargin{3})
        idx = varargin{3};
        [n,m]=size(idx) ;
        if N~=n; error('input ''X'' and ''idx'' must have the same number of rows'); end
        %if m~=1; error('input ''idx'' must be a column vector'); end
   end
end

xmin  = min(X)  ;
xmax  = max(X)  ;
X2 = (X-repmat(xmin,N,1))./(repmat(xmax-xmin,N,1)) ; % (N,M) 
X2 = X2' ; % (M,N)

% Create PCP:
if isempty(idx)
    m=0;
%    idx = true(N,1);
    lc = [0 0 0] ;
    lw = 1 ;
else
    if m==1
       lc  = [150 150 150]/256 ;
       lch = [0 0 0] ;
       lw  = 1 ;
       lwh = 2 ;
    else
       lc  = [150 150 150]/256 ;
       lw  = 1 ;
       lwh = 2 ;
%       lch = autumn(m);
       lch = lines(m) ;
    end
end
 
%figure
plot(X2,'Color',lc,'LineWidth',lw)
for j =1:m
    hold on
    plot(X2(:,idx(:,j)),'Color',lch(j,:),'LineWidth',lwh)
end

% Costumize horizontal axis:
set(gca,'XTick',1:M,'XTickLabel',XLabel,'XLim',[0.7 M+0.3])

% Add vertical grid:
set(gca,'XGrid','On','GridLineStyle','-')

% Costumize font
set(gca,'FontSize',fs,'FontName',fn)

% % Highlight in black some of the lines (if required):
% if ~isempty(idx)
%         hold on 
%         plot(X2(:,idx),'Color',lch,'LineWidth',lwh)
% end

% Label ticks on vertical axis:

if ~isempty(i_axis) % Show only the input range specified by user
    str=cell(1,10);
    [tmp,YTick] = hist(X(:,i_axis),10);
    for j=1:10; str{j}=num2str(YTick(j),yticklabels_form); end
    [tmp,YTick] = hist(X2(i_axis,:),10);
    set(gca,'YTick',YTick,'YTickLabel',str,'YLim',[-0.1,1.1])
    ylabel(XLabel{i_axis})

else % Show ranges of all inputs
    for i=1:M
        text(i+0.05,0,num2str(xmin(i), yticklabels_form),'HorizontalAlignment','Center','VerticalAlignment','Top','Fontsize',fs,'Fontname',fn)
        text(i+0.05,1,num2str(xmax(i), yticklabels_form),'HorizontalAlignment','Center','VerticalAlignment','Bottom','Fontsize',fs,'Fontname',fn)
        set(gca,'YTick',[],'YTickLabel',[],'YLim',[-0.1,1.1])
    end
end



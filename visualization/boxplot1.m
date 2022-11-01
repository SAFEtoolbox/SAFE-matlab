function [] = boxplot1(S,varargin)
%
% Plot a set of indices as coloured boxes. The height of the box is fixed
% (and thin) if the indices are given as deterministic values, and it is 
% equal to the uncertainty interval if the indices are associated with a 
% lower and upper bound.
%
%      boxplot1(S)
%      boxplot1(S,labels,Y_Label)
%      boxplot1(S,labels,Y_Label,S_lb,S_ub)
%
%       S = vector of indices                           - vector (1,M)
%  labels = strings for the x-axis labels               - cell array (1,M)
%           default: {'X1','X2',...,XM'}
% Y_Label = y-axis label  (default: 'Sensitivity')      - string
%    S_lb = lower bound of 'S' (default: empty)         - vector (1,M)
%    S_ub = upper bound of 'S' (default: empty)         - vector (1,M)
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

dh = 0.40 ; % semi-width of the box
dv = 0.01 ; % semi-height of the box (for deterministic value)

% Options for the colours:
ec  = 'k' ; % color of edges
% You can produce a coloured plot or a black and white one
% (printer-friendly). Furthermore, you can use matlab colourmaps or
% repeat 5 'easy-to-distinguish' colours (see http://colorbrewer2.org/).
% Option 1a - coloured using colorbrewer: uncomment the following line:
col = [[228,26,28];[55,126,184];[77,175,74];[152,78,163];[255,127,0]]/256;
% Option 1b - coloured using matlab colormap: uncomment the following line:
%col=hsv(length(S));
% Option 2a - B&W using matlab colorbrewer: uncomment the following line:
%col = [[37 37 37];[90 90 90];[150 150 150];[189 189 189];[217 217 217]]/256;
% Option 2b - B&W using matlab colormap: uncomment the following line:
%col=gray(length(S));


%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(S); error('''S'' must be a vector of size (1,M)'); end
[N,M] = size(S)  ;
if N~=1; error('''S'' must be a row vector'); end
if M<1 ; error('''S'' must have at least one component'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Set optional arguments to their default values:
labelinput = cell(1,M); for i=1:M; labelinput{i}=['X' num2str(i)]; end
Y_Label = 'Sensitivity';
S_lb    = nan(1,M) ; 
S_ub    = nan(1,M) ; 

% Recover and update optional arguments:
if nargin > 1
    if ~isempty(varargin{1})
        labelinput = varargin{1};
        if ~iscell(labelinput); error('''labelinput'' must be a cell array'); end
        if length(labelinput)~=M; error('''labelinput'' must have M=%d components',M); end
        for i=1:M; if ~ischar(labelinput{i}); error('all components of ''labelinput'' must be string'); end; end
    end
end
if nargin > 2
    if ~isempty(varargin{2})
        Y_Label = varargin{2};
        if ~ischar(Y_Label); error('''Y_Label'' must be a string'); end
    end
end
if nargin == 4
    error('If any of S_lb or S_ub is specified, than both must be specified')
end
if nargin > 3  
    if isempty(varargin{3})
        error('If any of S_lb or S_ub is specified, than both must be specified')
    else
        S_lb = varargin{3} ;
        if ~isnumeric(S_lb); error('''S_lb'' must be a vector of size (1,M)'); end
        if any(isnan(S_lb)); error('''S_lb'' has NaN elements'); end
        [n,m] = size(S_lb) ;
        if n~=1; error('''S_lb'' must be a row vector'); end
        if M~=m; error('''S'' and''S_lb'' must have the same number of elements'); end
        if any(S_lb>S); error('''S_lb'' must be lower or equal to ''S'''); end
    end
    if isempty(varargin{4})
        error('If any of S_lb or S_ub is specified, than both must be specified')
    else
        S_ub = varargin{4} ;
        if ~isnumeric(S_ub); error('''S_ub'' must be a vector of size (1,M)'); end
        if any(isnan(S_ub)); error('''S_ub'' has NaN elements'); end
        [n,m] = size(S_ub) ;
        if n~=1; error('''S_ub'' must be a row vector'); end
        if M~=m; error('''S'' and''S_ub'' must have the same number of elements'); end
        if any(S_ub<S); error('''S_ub'' must be larger or equal to ''S'''); end
    end
end
    
%%%%%%%%%%%%%%%%%%%
% Produce plot
%%%%%%%%%%%%%%%%%%%

[A,B]=size(col);
L=ceil(M/A);
fc=repmat(col,L,1);

%figure
for j=1:M
        if sum(isnan(S_lb)) % no confidence bounds
            % Plot the value as a tick line:
            patch([j-dh j-dh j+dh j+dh ],[S(j)-dv S(j)+dv S(j)+dv S(j)-dv],fc(j,:),'EdgeColor',ec,'LineWidth',1)
            hold on
        else
            % Plot the confidence interval as a rectangle:
            patch([j-dh j-dh j+dh j+dh ],[S_lb(j) S_ub(j) S_ub(j) S_lb(j)],fc(j,:),'EdgeColor',ec,'LineWidth',1.5)
            hold on
            % Plot the mean as a tick line:
            dv = 0.005 ;
            patch([j-dh j-dh j+dh j+dh ],[S(j)-dv S(j)+dv S(j)+dv S(j)-dv],ec,'EdgeColor',ec,'LineWidth',1)
        end
end

x1 = 0   ;
x2 = M+1 ;

plot(x1:x2,zeros(size(x1:x2)),':k')

set(gca,'FontSize',fs,'XTick',1:M,'XTickLabel',labelinput,'XGrid','On','FontName',fn)
ylabel(Y_Label,'FontSize',fs,'FontName',fn)
if sum(isnan(S_lb))
    y1 = min(-0.1, min(S));
else
    y1 = min(-0.1, min(S_lb));
end
if sum(isnan(S_ub))
    y2 = max(1.1, max(S));
else
    y2 = max(1.1, max(S_ub));
end
% axis([x1,x2,-0.1,1.1])
xlim([x1,x2])
ylim([y1,y2])
box on
hold off


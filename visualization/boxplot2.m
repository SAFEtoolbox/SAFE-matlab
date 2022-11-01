function [] = boxplot2(S,varargin)
%
% Plot a set of indices as coloured boxes. The height of the box is fixed
% (and thin) if the indices are given as deterministic values, and it is 
% equal to the uncertainty interval if the indices are associated with a 
% lower and upper bound.
%
%      boxplot1(S)
%      boxplot1(S,labels)
%      boxplot1(S,labels,S_lb,S_ub)
%
%      S = vector of indices                    - matrix (2,M)
% labels = strings for the x-axis labels        - cell array (1,M)
%          default: {'X1','X2',...,XM'}
%   S_lb = lower bound of 'S'                   - matrix (2,M)
%   S_ub = upper bound of 'S'                   - matrix (2,M)
%          default: none

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

% Options for the colours:
ec1  = 'k' ; % color of edges
ec2  = 'k' ; % color of edges
% You can produce a coloured plot or a black and white one
% (printer-friendly). 
% Option 1 - coloured using colorbrewer: uncomment the following lines:
fc1 = [239,138,98]/256;
fc2 = [103,169,207]/256;
% Option 2 - B&W using matlab colorbrewer: uncomment the following lines:
%fc1 = [150 150 150]/256;
%fc2 = [37 37 37]/256;

dh = 0.40 ; % semi-width of the box
dv = 0.01 ; % semi-height of the box (for deterministic value)


%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(S); error('''S'' must be a matrix of size (2,M)'); end
[N,M] = size(S)  ;
if N~=2; error('''S'' must have 2 rows'); end
if M<1 ; error('''S'' must have at least one column'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Set optional arguments to their default values:
labelinput = cell(1,M); for i=1:M; labelinput{i}=['X' num2str(i)]; end
S_lb    = nan(2,M) ; 
S_ub    = nan(2,M) ; 

% Recover and update optional arguments:
if nargin > 1
    if ~isempty(varargin{1})
        labelinput = varargin{1};
        if ~iscell(labelinput); error('''labelinput'' must be a cell array'); end
        if length(labelinput)~=M; error('''labelinput'' must have M=%d components',M); end
        for i=1:M; if ~ischar(labelinput{i}); error('all components of ''labelinput'' must be string'); end; end
    end
end
if nargin == 3
    error('If any of S_lb or S_ub is specified, than both must be specified')
end
if nargin > 2  
    if isempty(varargin{2})
        error('If any of S_lb or S_ub is specified, than both must be specified')
    else
        S_lb = varargin{2} ;
        if ~isnumeric(S_lb); error('''S_lb'' must be a matrix of size (2,M)'); end
        if any(isnan(S_lb(:))); error('''S_lb'' has NaN elements'); end
        [n,m] = size(S_lb) ;
        if n~=2; error('''S_lb'' have 2 rows'); end
        if M~=m; error('''S'' and''S_lb'' must have the same number of columns'); end
        if any(S_lb(1,:)>S(1,:)); error('''S_lb'' must be lower or equal to ''S'''); end
        if any(S_lb(2,:)>S(2,:)); error('''S_lb'' must be lower or equal to ''S'''); end
    end
    if isempty(varargin{3})
        error('If any of S_lb or S_ub is specified, than both must be specified')
    else
        S_ub = varargin{3} ;
        if ~isnumeric(S_ub); error('''S_ub'' must be a matrix of size (2,M)'); end
        if any(isnan(S_ub(:))); error('''S_ub'' has NaN elements'); end
        [n,m] = size(S_ub) ;
        if n~=2; error('''S_ub'' have 2 rows'); end
        if M~=m; error('''S'' and''S_ub'' must have the same number of columns'); end
        if any(S_ub(1,:)<S(1,:)); error('''S_ub'' must be larger or equal to ''S'''); end
        if any(S_ub(2,:)<S(2,:)); error('''S_ub'' must be larger or equal to ''S'''); end
    end
end
    
%%%%%%%%%%%%%%%%%%%
% Produce plot
%%%%%%%%%%%%%%%%%%%

%figure
for j=1:M
        if sum(isnan(S_lb)) % no confidence bounds
            % Plot the mean as a tick line:
            patch([j-dh j-dh j-0.05 j-0.05 ],[S(1,j)-dv S(1,j)+dv S(1,j)+dv S(1,j)-dv ],fc1,'EdgeColor',ec1)
            hold on
            patch([j+0.05 j+0.05 j+dh j+dh ],[S(2,j)-dv S(2,j)+dv S(2,j)+dv S(2,j)-dv ],fc2,'EdgeColor',ec2)
        else
            % Plot the confidence interval as a rectangle:
            patch([j-dh j-dh j-0.05 j-0.05 ],[S_lb(1,j) S_ub(1,j) S_ub(1,j) S_lb(1,j)],fc1,'EdgeColor',ec1,'LineWidth',1.5)
            hold on
            patch([j+0.05 j+0.05 j+dh j+dh ],[S_lb(2,j) S_ub(2,j) S_ub(2,j) S_lb(2,j)],fc2,'EdgeColor',ec2,'LineWidth',1.5)
            % Plot the mean as a tick line:
            dv = 0.005 ;
            patch([j-dh j-dh j-0.05 j-0.05 ],[S(1,j)-dv S(1,j)+dv S(1,j)+dv S(1,j)-dv ],ec1,'EdgeColor',ec1)
            patch([j+0.05 j+0.05 j+dh j+dh ],[S(2,j)-dv S(2,j)+dv S(2,j)+dv S(2,j)-dv ],ec2,'EdgeColor',ec2)
        end
end

x1 = 0   ;
x2 = M+1 ;

plot(x1:x2,zeros(size(x1:x2)),':k') % plot a horizontal line at zero

for i=0:M; hold on; plot([i+0.5 i+0.5],[-1,2],'k'); end % plot vertical lines between one input and another


set(gca,'FontSize',fs,'XTick',1:M,'XTickLabel',labelinput,'XGrid','On','FontName',fn)
ylabel('Sensitivity','FontSize',fs,'FontName',fn)
% axis([x1,x2,-0.1,1.1])
if sum(isnan(S_lb))
    y1 = min(-0.1, min(S(:)));
else
    y1 = min(-0.1, min(S_lb(:)));
end
if sum(isnan(S_ub))
    y2 = max(1.1, max(S(:)));
else
    y2 = max(1.1, max(S_ub(:)));
end

xlim([x1,x2])
ylim([y1,y2])
ylim
box on
hold off


function [] = plot_convergence(S,NN,varargin)
%
% Plots sequence of one or more indices estimated using samples
% of different sizes
%
% Usage: 
% plot_convergence(S,NN)
% plot_convergence(S,NN,S_lb,Sd_ub) 
% plot_convergence(S,NN,S_lb,Sd_ub,SExact) 
% plot_convergence(S,NN,S_lb,Sd_ub,SExact,X_Label,Y_Label,labelinput)
%
% Inputs:
%           S = sequence of estimates of M indices           - matrix (R,M)
%          NN = vector of sample sizes at which indices were estimated
%                               - vector (1,R) of increasing integer values
% Optional inputs
%        S_lb = lower bound of (uncertain) index estimates   - matrix (R,M)
%        S_ub = upper bound of (uncertain) index estimates   - matrix (R,M)
%      SExact = exact value of the indices (if known)        - vector (1,M)
%     X_Label = x-axis label                                 - string
%     Y_Label = y-axis label                                 - string
%  labelinput = legend labels                            - cell array (1,M)
%
% (leave empty any optional argument that you don't want to specify) 

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

% Options for the legend:
sorting   = 1  ; % If 1, inputs will be displayed in the legend 
% according to their influence, i.e. from most sensitive to least sensitive
% (if 0 they will be displayed according to their original order)
nb_legend = 5  ; % number of input names that will be displayed in the legend
end_length=0.3; %adjust the space left for the legend

% Options for the colours:
% You can produce a coloured plot or a black and white one
% (printer-friendly). Furthermore, you can use matlab colourmaps or
% repeat 5 'easy-to-distinguish' colours (see http://colorbrewer2.org/).
% Option 1a - coloured using colorbrewer: uncomment the following line:
col = [[228,26,28];[55,126,184];[77,175,74];[152,78,163];[255,127,0]]/256;
% Option 1b - coloured using matlab colormap: uncomment the following line:
%col=hsv(size(S,2));
% Option 1a - B&W using matlab colorbrewer: uncomment the following line:
%col = [[37 37 37];[90 90 90];[150 150 150];[189 189 189];[217 217 217]]/256;
% Option 1b - B&W using matlab colormap: uncomment the following line:
%col=gray(size(S,2)); 

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(S); error('''S'' must be a matrix of size (R,M)'); end
[Ri,M] = size(S)  ;
if ~isnumeric(NN); error('input ''NN'' must be a vector of integer numbers'); end
if ~isnumeric(NN); error('input ''NN'' must be a vector of integer numbers'); end
[n,R] = size(NN);
if n~=1; error('''NN'' must be a row vector'); end
if sum(abs(NN-round(NN))); error('all elements of ''NN'' must be integer numbers'); end
if any(NN<=0); error('all elements in ''NN'' must positive'); end
if any(diff(NN)<=0); error('elements in ''NN'' must be sorted in ascending order'); end
if R~=Ri; error('number of columns in ''S'' must be equal to the number of elements in ''NN'''); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
S_lb  = zeros(R,M);
S_ub  = zeros(R,M);
SExact = [];
% X_Label = 'sample size';
% Y_Label = 'sensitivity' ;
X_Label = '';
Y_Label = '' ;
labelinput = cell(1,M); for i=1:M; labelinput{i}=['#' num2str(i)]; end

% Recover and update optional arguments:
if nargin > 2
    if ~isempty(varargin{1})
        S_lb = varargin{1} ;
        if ~isnumeric(S_lb); error('''S_lb'' must be a vector of size (1,M)'); end
        [n,m] = size(S_lb) ;
        if n~=R; error('''S'' and''S_lb'' must have the same number of rows'); end
        if M~=m; error('''S'' and''S_lb'' must have the same number of columns'); end
    end
end
if nargin > 3
     if ~isempty(varargin{2})
        S_ub = varargin{2} ;
        if ~isnumeric(S_ub); error('''S_ub'' must be a vector of size (1,M)'); end
        [n,m] = size(S_ub) ;
        if n~=R; error('''S'' and''S_ub'' must have the same number of rows'); end
        if M~=m; error('''S'' and''S_ub'' must have the same number of columns'); end
    end
end

if nargin > 4
    if ~isempty(varargin{3})
        SExact = varargin{3};
        [n,Mi] = size(SExact);
        if n~=1; error('''SExact'' must be a row vector'); end
        if Mi~=M; error('''S'' and''SExact'' must have the same number of columns'); end
    end
end
if nargin > 5
    if ~isempty(varargin{4})
        X_Label = varargin{4};
        if ~ischar(X_Label); error('''X_Label'' must be a string'); end
    end
end
if nargin > 6
    if ~isempty(varargin{5})
        Y_Label = varargin{5};
        if ~ischar(Y_Label); error('''Y_Label'' must be a string'); end
    end
end
if nargin > 7
    if ~isempty(varargin{6})
        labelinput = varargin{6};
        if ~iscell(labelinput); error('''labelinput'' must be a cell array'); end
        if length(labelinput)~=M; error('''labelinput'' must have M=%d components',M); end
        for i=1:M; if ~ischar(labelinput{i}); error('all components of ''labelinput'' must be string'); end; end
    end
end

%%%%%%%%%%%%%%
% plot results
%%%%%%%%%%%%%%

[A,B]=size(col);
L=ceil(M/A);
clrs=repmat(col,L,1);

% Set horizontal and vertical limits:
if NN(1) - mean(diff(NN))>0; H1 = NN(1)- mean(diff(NN)) ; else ; H1=0 ; end
H2 = NN(end)+end_length*(NN(end)-NN(1)) ;

if  any(any(S_lb))~=0
    V1 = min( 0, min(min(S_lb)) ) ;
else
    V1 = min( 0, min(min(S)) ) ;
end
if  any(any(S_ub))~=0
    V2 = max( 1, max(max(S_ub)) ) ;
else
    V2 = max( 1, max(max(S)) ) ;
end


labelinput_new=cell(1,M);

if sorting
    [tmp,Sidx]=sort(S(end,:),'descend') ;
    S = S(:,Sidx) ;
    S_ub = S_ub(:,Sidx) ;
    S_lb = S_lb(:,Sidx) ;
    for i=1:M; labelinput_new{i} = labelinput{Sidx(i)} ;end;
    if ~isempty(SExact)
       SExact = SExact(Sidx) ;
    end
end
%
if nb_legend<M
    labelinput_new=labelinput_new(1:nb_legend);
    labelinput_new{end}=[labelinput_new{end},'...'];
end


% For each index, plot final estimated value:
for i=1:M
    plot(NN(end),S(end,i),'o','MarkerSize',10,'MarkerFaceColor',clrs(i,:),'MarkerEdgeColor','k')
    hold on
end

% Draw an horizontal line at zero:
plot([H1 H2],[0 0],'k')


for i=1:M
    % Plot trajectory with increasing number of samples:
    plot(NN,S(:,i),'Color',clrs(i,:),'Linewidth',2.5);
    box 'on'    
     %Plot exact index value (if known):
    if ~isempty(SExact)
     
      %  plot([NN(end) NN(end)+mean(diff(NN))],[S(end,i) SExact(i)],'o:','MarkerSize',10,'MarkerFaceColor',tmp(clrs(i),:),'Color',tmp(clrs(i),:))
        plot([H1 H2],[SExact(i) SExact(i)],'--','Color',clrs(i,:),'LineWidth',2)
    end
    
end


% Plot confidence bounds

if  any(any(S_lb))~=0
    for i=1:M
        plot(NN,S_lb(:,i),'Color',clrs(i,:),'LineStyle','--','Linewidth',1.2)
    end
end
if  any(any(S_ub))~=0
    for i=1:M
        plot(NN,S_ub(:,i),'Color',clrs(i,:),'LineStyle','--','Linewidth',1.2)
    end
end
% Axes labels:

xlabel(X_Label,'FontName',fn,'FontSize',fs)
ylabel(Y_Label,'FontName',fn,'FontSize',fs)


legend(labelinput_new, 'Location','NorthWest')



%Tick labels for horizontal axis:
str=cell(size(NN)); for k=1:R; str{k}=num2str(NN(k)); end
axis([H1,H2,V1,V2])
set(gca,'XTick',NN,'XTickLabel',str,'FontName',fn,'FontSize',fs)
grid on


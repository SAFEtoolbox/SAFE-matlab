function [] = stackedbar(S,varargin)
%
% Plot and compare N different set of indices using stacked bar plot. 
%
%      stackedbar(S)
%      stackedbar(S,labelinput,Y_Label)
%      stackedbar(S,labelinput,Y_Label,horiz_tick)
%      stackedbar(S,labelinput,Y_Label,horiz_tick,horiz_tick_label)
%
% S = set of indices (N sets, each composed of M indices)    - matrix (N,M)
% labelinput = names of inputs to appear in the legend   - cell array (1,M)
%               (default: no legend)
% Y_Label = label for vertical axis  (default: 'Sensitivity')      - string
% horiz_tick = ticks for horizontal axis                     - vector (1,N)     
%              (default: [1,2,...,N] )
% horiz_tick = labels for ticks of horizontal axis       - cell array (1,N)     
%              (default: {'1','2',...,'N'} )
%
% Example:
% S = rand(5,4);
% figure; stackedbar(S)
% figure; stackedbar(S,{'a','b','c','d'})
% figure; stackedbar(S,{'a','b','c','d'},'total')
% figure; stackedbar(S,[],'total',[1,2,9,11,12])
% figure; stackedbar(S,[],'total',[1,2,9,11,12],{'d1','d2','d3','d4','d5'})

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

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(S); error('''S'' must be a vector of size (1,M)'); end
[N,M] = size(S)  ;
%if N~=1; error('''S'' must be a row vector'); end
if M<1 ; error('''S'' must have at least one component'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Set optional arguments to their default values:
labelinput = [];
Y_Label = 'Sensitivity';
horiz_tick = 1:N ;
horiz_tick_label=[];

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
if nargin > 3
     if ~isempty(varargin{3})
        horiz_tick = varargin{3};
        if ~isnumeric(horiz_tick); error('''horiz_tick'' must be a vector of length N=%d',N); end
        if length(horiz_tick(:))~=N ; error('''horiz_tick'' must be a vector of length N=%d',N); end
    end
end

if nargin > 4
    if ~isempty(varargin{4})
        horiz_tick_label = varargin{4};
        if ~iscell(horiz_tick_label); error('''horiz_tick_label'' must be a cell array'); end
        if length(horiz_tick_label)~=N; error('''horiz_tick_label'' must have N=%d components',N); end
        for i=1:N; if ~ischar(horiz_tick_label{i}); error('all components of ''horiz_tick_label'' must be string'); end; end
    end
end

%%%%%%%%%%%%%%%%%%%
% Produce plot
%%%%%%%%%%%%%%%%%%%

% figure
bar(horiz_tick,S,'stacked')
set(gca,'FontSize',fs,'XTick',horiz_tick,'XGrid','On','FontName',fn)
if ~isempty(horiz_tick_label); set(gca,'XTickLabel',horiz_tick_label); end
if ~isempty(labelinput); legend(labelinput); end
ylabel(Y_Label,'FontSize',fs,'FontName',fn)
box on


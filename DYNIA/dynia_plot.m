function hfig = dynia_plot(xi,fi,varargin)
%
% Plot time pattern of posterior distribution of a model input
% (posterior distribution is approximated by histogram)
%
% Usage:
% dynia_plot(xi,fi)
% dynia_plot(xi,fi,X_Label)
% dynia_plot(xi,fi,X_Label,y)
%
% Input:
% xi = bins centers                                       - vector (nbin,1)
% fi = frequency associated to each bin at each time step - matrix (nbin,T)
% X_Label = label for horizontal axis                              - string
%           default: 'input range'
% y = time series of model output                            - vector (T,1)
%     (to be plotted on the top of the input posterior probability)
%
% See also dynia_histograms for further information and interpretation.

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
fs = 24 ; % font size of axes, labels, etc.
nb_fig=2; %number of digits after decimal points for vertical axis

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(xi) ; error('input ''xi'' must be a vector of real numbers'); end
if ~isnumeric(fi) ; error('input ''fi'' must be a matrix of real numbers'); end
[N,M]=size(xi) ;
[n,T]=size(fi) ;
if M~=1; error('input ''xi'' must be a column vector'); end
if N~=n; error('input ''xi'' and ''fi'' must have the same number of rows'); end
if T<=1; error('input ''fi'' must be a matrix with more than one column'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:

X_Label = 'input range' ;
y = [] ;

% Recover and update optional arguments:

if nargin > 2
    if ~isempty(varargin{1})
        X_Label = varargin{1};
        if ~ischar(X_Label); error('''X_Label'' must be a string'); end
    end
end
if nargin > 3
    if ~isempty(varargin{2})
        y = varargin{2};
        y = y(:) ; 
        if length(y)~=T; error('''y'' must be a vector with the same number of components as the rows in ''fi'''); end
    end
end

%%%%%%%%%%%%%%%%
% Generate plots
%%%%%%%%%%%%%%%%

nbins = length(xi) ;

% Pre-processing: invert the order of the rows in 'xi' and 'fi' so that
% data relevant to higher parameter values be on the first rows (and
% thus will be displaced at the top of the vertical axis by the
% function 'imagesc' used later)
[xi_plot,idx] = sort(xi,'descend') ;
fi_plot = fi(idx,:) ;
%xi_plot = xi(end:-1:1)   ;
%fi_plot = fi(end:-1:1,:) ;

% add two fake columns filled in with ones and zeros so that 'imagesc'
% will rescale to the range [0,1] (0 and 1 having now become the
% minimum and maximum of 'fi_plot'):
fi_plot = [ fi_plot ones(nbins,1) zeros(nbins,1) ] ;

% Plot results:
hfig = figure;
%clrs = gray ;
clrs = autumn ;
% invert ordering of row sin 'clrs' so that black colour is associated
% to value of 1 and white to the value of 0:
clrs = clrs(end:-1:1,:);
colormap(clrs)
% plot the frequencies:
imagesc( fi_plot )
% set axes, labels, etc.:
axis([1 T 0.5 nbins+0.5])
xlabel('time','FontName',fn,'FontSize',fs)
ylabel(X_Label,'FontName',fn,'FontSize',fs)
str=cell(1,nbins);
%for j=1:nbins; str{j}=num2str(xi_plot(j)); end
for j=1:nbins; str{j}=num2str(round(10^nb_fig*xi_plot(j))/(10^nb_fig)); end

set(gca,'YTick',1:nbins,'YTickLabel',str)
c = colorbar;
set(get(c,'title'),'String','freq','FontName',fn,'FontSize',fs)

% Plot (rescaled) output values (y) on the top:
hold on
y_norm = y/max(y) ; % now 'y' varies between zero and one
plot(nbins+0.5-y_norm*(nbins-1),'k','LineWidth',2)

set(gca,'FontName',fn,'FontSize',fs)

function [fi,yi] = plot_pdf(y,varargin)
%
% This function plot the empirical Probability Distribution Function (PDF)
% of a given sample (y). The empirical PDF is approximated by the histogram
% 
% plot_pdf(y)
% plot_pdf(y,Y_Label)
% 
% Input:
%       y = sample variable values         - matrix/vector %
% Y_label = label for horizontal axis      - string        %
%    nbin = number of bins for histogram   - scalar        %
%           (default: min(100,length(y)/10))               %
% Output:
%      fi = frequency of each bin        - vector (1,nbin) %
%      yi = center of each bin           - vector (1,nbin) %
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

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(y); error('''y'' must be a matrix or a vector'); end

if all(isnan(y)); error('all data in ''y'' are NaN'); end
if all(isinf(y)); error('all data in ''y'' are Inf'); end
  
if any(any(isnan(y))); disp('Warning: some data in ''y'' are NaN'); end
if any(any(isinf(y))); disp('Warning: some data in ''y'' are Inf'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Set optional arguments to their default values:
Y_Label = 'y' ;
y = y(:)     ;
n = length(y);
nbin = max(5,min(100,n/10));

% Recover and update optional arguments:
if nargin > 1
    if ~isempty(varargin{1})
        Y_Label = varargin{1} ;
        if ~ischar(Y_Label); error('''Y_Label'' must be a string'); end;
    end
end
if nargin > 2
    if ~isempty(varargin{2})
        nbin = varargin{2} ;
        if ~isscalar(nbin); error('''nbin'' must be scalar'); end
        if nbin<+1; error('''nbin'' must be >1' ); end
        if abs(nbin-round(nbin)); error('''nbin'' must be an integer'); end
    end
end

%%%%%%%%%%%%%%%%%%%
% Produce plot
%%%%%%%%%%%%%%%%%%%

[ni,yi] = hist(y,nbin);
fi = ni/(n*(max(y)-min(y))/nbin) ;

% Limit for horizontal axis:
ym = min(y);
yM = max(y);
if ym==yM % (i.e., if all data have the same value)
    ym = ym - ym/10 ;
    yM = yM + yM/10 ;
end
% 
% Limit for vertical axis:    
fm = 0 ;
fM = max(fi);
% if mean(fi)/std(fi)>2; % if frequency is rather constant (''flat''),
%     % expand the vertical axis so to highlight this:
%     fM = min(1,max(fi)+10*mean(fi)) ; 
% else
%     fM = min(1,max(fi)+std(fi))     ;
% end

% Plot:
plot(yi,fi,'.-','LineWidth',2) 
set(gca,'FontName',fn,'FontSize',fs)
axis([ym,yM,fm,fM])
xlabel(Y_Label,'FontName',fn,'FontSize',fs)
ylabel('PDF','FontName',fn,'FontSize',fs)
box on

    
    
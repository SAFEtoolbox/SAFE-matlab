function [] = plot_cdf(y,varargin)
%
% This function plot the empirical Cumulative Distribution Function (CDF)
% of a given sample (y).
%
% plot_cdf(y)
% plot_cdf(y,Y_Label)
%
%       y = sample variable values    - matrix/vector
% Y_label = label for horizontal axis - string
%
% NOTE: If Y includes any NaN values, the function  will identify them from
% the CDF.

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

% Recover and update optional arguments:
if nargin > 1
    if ~isempty(varargin{1})
        Y_Label = varargin{1} ;
        if ~ischar(Y_Label); error('''Y_Label'' must be a string'); end;
    end
end

%%%%%%%%%%%%%%%%%%%
% Produce plot
%%%%%%%%%%%%%%%%%%%

y  = y(:)      ;
y = y(~isnan(y));% remove NaNs

%n  = length(y) ;
%F  = [1:n]/n   ;
%ys = sort(y)   ;
%plot(ys,F,'.-','LineWidth',2)

Nmax = 5000 ;
if length(y)>Nmax
    yi = sort(y) ;
    F = empiricalcdf(y,yi);
    plot(yi,F,'.')
else
    yi=min(y):(max(y)-min(y))/Nmax:max(y);
    F = empiricalcdf(y,yi);
    plot(yi,F)
end

% Customize plot:
set(gca,'FontName',fn,'FontSize',fs)
xlabel(Y_Label,'FontName',fn,'FontSize',fs)
ylabel('CDF','FontName',fn,'FontSize',fs)
box on

% Limit for horizontal axis:
ym = min(y);
yM = max(y);
if ym==yM % (i.e., if all data have the same value)
    ym = ym - ym/10 ;
    yM = yM + yM/10 ;
end
axis([ym,yM,0,1])



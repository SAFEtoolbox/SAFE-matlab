function [] = pawn_plot_kstest(KSi,NC,NU,varargin)
%
% Visual version of the two-sample Kolmogorov-Smirnov test for inputs 
% screening, according to the PAWN approach (Pianosi and Wagener, 2015).
% Plot the Kolomogorov-Smirnov statistic between unconditional output CDF 
% (i.e. when all inputs vary) and the conditional CDFs (when one or more
% inputs are fixed) at different conditioning values for the fixed
% input(s), and compare it with the KS critical value at given level of
% confidence.
%
% Usage:
% pawn_plot_kstest(KSi,NC,NU,alfa)
% pawn_plot_kstest(KSi,NC,NU,alfa,xci)
% pawn_plot_kstest(KSi,NC,NU,alfa,xci,x_label)
%
% Input:
%  KSi = Kolmogorov-Smirnov statistic                        - vector (n,1)
%        between empirical unconditional CDF of the output
%        and conditional CDF when fixing the i-th input
%        at its k-th (k=1,...,n) conditioning value 
%   NC = number of output samples used to compute the              - scalar
%        conditional CDFs
%   NU = number of output samples used to compute the              - scalar 
%        unconditional CDF
% alfa = significance level of the KS test                         - scalar
%        must take values in among {0.1,0.05,0.025,0.01,0.005,0.001}
%        If alfa is empty, then all critical values will be plotted 
%        (6 lines)
% xci = conditioning values for the fixed input              - vector (n,1)
%       (default: 1:n)
% x_label = legend for the horizontal axis                         - string
%           (default: 'conditioning points')
%
% REFERENCES:
%
% Pianosi, F. and Wagener, T. (2015), A simple and efficient method 
% for global sensitivity analysis based on cumulative distribution 
% functions, Env. Mod. & Soft., 67, 1-11.

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

fn='Helvetica';
fs = 13 ;
ms=8; % Marker size

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(KSi); error('''KSi'' must be a vector'); end
[n,M] = size(KSi) ;
if M>1; error('''KSi'' must be a column vector'); end
    
if ~isscalar(NC); error('''NC'' must be a scalar'); end
if NC<=0; error('''NC'' must be positive' ); end
if abs(NC-round(NC)); error('''NC'' must be integer'); end

if ~isscalar(NU); error('''NU'' must be a scalar'); end
if NU<=0; error('''NU'' must be positive' ); end
if abs(NU-round(NU)); error('''NU'' must be integer'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xci = 1:n ;
x_label = 'conditioning points' ;
alfa =[ 0.10 0.05 0.025 0.01 0.005 0.001];
c =[ 1.22 1.36 1.48  1.63 1.73  1.95 ]; 
n_alfa=6; %number of significance levels

if nargin>3
    if ~isempty(varargin{1})
        alfa=varargin{1} ;
        if ~isnumeric(alfa); error('''alfa'' must be a vector'); end
        [m,n_alfa]=size(alfa);
        if m~=1;error('''alfa'' must be a scalar');end
        if n_alfa~=1;error('''alfa'' must be a scalar');end
        if     alfa==0.10 ; c=1.22;
        elseif alfa==0.05 ; c=1.36;
        elseif alfa==0.025; c=1.48;
        elseif alfa==0.01 ; c=1.63;
        elseif alfa==0.005; c=1.73;
        elseif alfa==0.001; c=1.95;
        else
            error('''alfa'' must take values in {0.1,0.05,0.025,0.01,0.005,0.001}');
        end
    end
end

if nargin>4
    if ~isempty(varargin{2})
        xci = varargin{2} ;
        if ~isnumeric(xci) ; error('xci must be a column vector of %d elements',n); end
        [n1,m1] = size(xci) ;
        if m1~=1 ; error('xci must be a column vector'); end
        if n1~=n ; error('xci must be a column vector of %d elements',n); end
    end
end
if nargin>5
    if ~isempty(varargin{3})
        x_label = varargin{3} ;
        if ~ischar(x_label); error('''x_label'' must be a string'); end
    end
end

%%%%%%%%%%%%%%
% Create plot
%%%%%%%%%%%%%%

[xi,ii] = sort(xci);
ksi     = KSi(ii);

%figure

% plot horizontal line at the KS critical value
col=hsv(n_alfa);
ks_crit=nan(1,n_alfa);
for i=1:n_alfa
    ks_crit(i) = c(i)*sqrt((NC+NU)/(NC*NU)) ;
    plot(xi,ks_crit(i)*ones(n,1),'--k','LineWidth',2,'Color',col(i,:))
    hold on
end
%str=cell(1,n_alfa);for i=1:n_alfa; str{i}=num2str(alfa(i));end

if n_alfa >1
    [LEGH,OBJH,OUTH,OUTM]=legend([num2str(alfa(1))]);
    for i=2:n_alfa
        [LEGH,OBJH,OUTH,OUTM]=legend(OUTM{:},[num2str(alfa(i))]);
    end
    title = get(LEGH,'Title');
    set(title,'String','Significance level','FontName',fn,'Fontsize',fs)
end

plot(xi,ksi,'-k','LineWidth',1)
hold on

% plot KS values as coloured circles on a gray scale:
idx=gray(n);
for k=1:n   
    plot(xi(k),ksi(k),'ok','MarkerFaceColor',idx(k,:),'MarkerSize',ms)
    %plot(xi(k),ksi(k),'xk','MarkerFaceColor',idx(k,:),'MarkerSize',12)
end

%axis([xi(1),xi(end),max(min(min(ksi),min(ks_crit))-0.1,0),min(max(max(ksi),max(ks_crit))+0.1,1)])
axis([xi(1),xi(end),0,1])
set(gca,'FontSize',fs)

ylabel('KS','FontSize',fs)
xlabel(x_label,'FontSize',fs)

box on
        

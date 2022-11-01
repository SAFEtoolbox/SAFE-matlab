% This script provides an application example of VBSA where the original
% uncertain model inputs are grouped into a lower number of uncertain
% groups.
% As an example, we use the HBV rainfall-runoff model combined with a snow
% melting routine, and apply SA to investigate propagation of parameter
% uncertainty. The total number of parameters of this model is 13, however
% we here group them into 3 groups: parameters referring to the snowmelt
% routine, parameters of the soil moisture accounting component and
% parameters of the flow routing component.
%
% See also 'workflow_vbsa_hymod.m' about VBSA
% and 'workflow_eet_hbv.m' about HBV model and study area.

% This script prepared by Francesca Pianosi
% University of Bristol, 2015
% mail to: francesca.pianosi@bristol.ac.uk


%% Step 1 - add paths

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath([ my_dir '/sampling'])
addpath([ my_dir '/util'])
addpath([ my_dir '/visualization'])
addpath([ my_dir '/example/hbv'])
addpath([ my_dir '/VBSA'])

%% Step 2 - setup the HBV model and define model output

% Define option for model structure:
% Case = 1 > interflow is dominant
% Case = 2 > percolation is dominant
Case  =  1; 

% Setup warmup period:
warmup = 30 ; % (days)

% Load data:
load 01055500.txt
% Define initial and final day of simulation horizon:
date=[ 1948 10 1; 1950 09 30]; 
% Extract time series:
t_start=find(ismember(X01055500(:,1:3),date(1,:),'rows')>0);
t_end  =find(ismember(X01055500(:,1:3),date(2,:),'rows')>0);
time   =t_start:t_end;
prec = X01055500(time,4);
ept  = X01055500(time,5);
flow = X01055500(time,6);
temp = mean(X01055500(time,7:8),2);

% Define model output function:
myfun = 'hbv_snow_objfun';

%% Step 3 - define input ranges

% Snow routine parameters:
% 1. Ts     = threshold temperature [C]
% 2. CFMAX  = degree day factor [mm/C]
% 3. CFR    = refreezing factor [-]
% 4. CWH    = Water holding capacity of snow [-]       

% HBV parameters (soil moisture account):
% 5. BETA   = Exponential parameter in soil routine [-]
% 6. LP     = evapotranspiration limit [-]
% 7. FC     = field capacity [mm] 

% HBV parameters (flow routing):
% 8. PERC   = maximum flux from Upper to Lower Zone [mm/Dt]
% 9. K0     = Near surface flow coefficient (ratio) [1/Dt]  
% 10. K1    = Upper Zone outflow coefficient (ratio) [1/Dt]  
% 11. K2    = Lower Zone outflow coefficient (ratio) [1/Dt]  
% 12. UZL   = Near surface flow threshold [mm]
% 13. MAXBAS= Flow routing coefficient [Dt]

% Define parameter ranges (from Kollat et al., 2012):
xmin =[-3  0 0 0   0 0.3    1   0 0.05 0.01 0.05   0 1 ];
xmax =[ 3 20 1 0.8 7   1 2000 100 2    1    0.1  100 6 ];

% Define parameter distributions:
Mi = length(xmin) ;
DistrFun  = cell(1,Mi) ; 
for i=1:Mi-1;  DistrFun{i}='unif';  end; % all uniform...
DistrFun{end}='unid'; % ... but the last one, which is discrete uniform

% Define distribution parameters:
DistrPar = cell(Mi,1); for i=1:Mi-1; DistrPar{i} = [ xmin(i) xmax(i) ]; end
% all taken from the ranges, but the last one:
DistrPar{end} = xmax(end) ; % parameter MAXBAS is an integer


%% Step 4 - define groups

% Define grouping:
groups = { [1,2,3,4] , [5,6,7] , [8,9,10,11,12,13] } ;

% Choose names of groups (needed for plots):
G_labels = {'snow','soil','routing'} ;

Mg = length(G_labels) ; % number of uncertain groups

%% Step 5 - build catalogue

% Choose length of each catalogue:
Ng = [ 10000 , 10000 , 100000 ] ;

% Build the catalogue (via sampling):
samp_strat = 'rsu' ;
CAT = cell(1,Mg) ;
for i=1:Mg
    CAT{i} = AAT_sampling(samp_strat,length(groups{i}),{DistrFun{groups{i}}},{DistrPar{groups{i}}},Ng(i));
end

%% Step 6 - perform sampling 

% Choose sampling strategy:
SampStrategy = 'lhs' ;

% Choose base sample size:
N = 3000 ;

% Set distribution of the catalouge index variables to discrete uniform:
DistrFunG = cell(1,Mg) ; for i=1:Mg; DistrFunG{i}='unid';  end;
% Set distribution parameter to number of elements in the catalogue:
DistrParG = cell(1,Mg) ; for i=1:Mg; DistrParG{i} = Ng(i) ; end;

% Perform sampling (over indices of catalogue):
X = AAT_sampling(SampStrategy,Mg,DistrFunG,DistrParG,2*N);

% Create additional samples through resampling strategy:
[ XA, XB, XC ] = vbsa_resampling(X) ;

% Transform back to original inputs by reading the catalogue:
XAp = nan(N,Mi) ;
XBp = nan(N,Mi) ;
XCp = nan(N*Mg,Mi) ;
for i=1:Mg
    CAT_i = CAT{i}; 
    XAp(:,groups{i}) = CAT_i(XA(:,i),:) ;
    XBp(:,groups{i}) = CAT_i(XB(:,i),:) ;
    XCp(:,groups{i}) = CAT_i(XC(:,i),:) ;
end

%% Step 7 - run the model

YA = model_execution(myfun,XAp,prec,temp,ept,flow,warmup,Case) ; % size (N,1)
YB = model_execution(myfun,XBp,prec,temp,ept,flow,warmup,Case) ; % size (N,1)
YC = model_execution(myfun,XCp,prec,temp,ept,flow,warmup,Case) ; % size (N*M,1)

%% Step 8 - compute sensitivity indices

% Choose which model output to analyze:
j=2; % 1:AME, 2:NSE, 3:BIAS, 4:TRMSE, 5:SFDCE, 6:RMSE,
% see help of 'hbv_snow_objfun.m' for more information

% Compute main (first-order) and total effects:
[ Si, STi ] = vbsa_indices(YA(:,j),YB(:,j),YC(:,j));

% Plot results (plot main and total separately):
figure 
subplot(121); boxplot1(Si,G_labels,'main effects')
subplot(122); boxplot1(STi,G_labels,'total effects')

% Plot results (both in one plot):
figure
boxplot2([Si; STi],G_labels)
legend('main effects','total effects')% add legend



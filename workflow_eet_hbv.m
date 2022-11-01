% This script provides an application example of the 
% Elementary Effects Test to the HBV rainfall-runoff model.
% Useful to learn about how to use the EET when one of the parameter takes
% discrete values.
%
% METHOD
%
% see description in 'workflow_eet_hymod.m'
%
% MODEL AND STUDY AREA
%
% The model under study is the HBV rainfall-runoff model,
% the inputs subject to SA are the 13 model parameters,
% and the outputs for SA are a set of performance metrics.
% See help of function 'hbv_snow_objfun.m' for more details
% and references about the HBV model and the objective functions.
%
% The case study area is is the Nezinscot River at Turner center, Maine, 
% USA (USGS 01055500, see http://waterdata.usgs.gov/nwis/nwismap)

% This script prepared by Francesca Pianosi and Fanny Sarrazin
% University of Bristol, 2014
% mail to: francesca.pianosi@bristol.ac.uk


%% Step 1 (add paths)
my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%% Step 2 (setup the HBV model)

% Load data:
load 01055500.txt


Case  =  1; % Case=1: interflow is dominant / Case = 2: percolation is dominant
warmup=365; % Model warmup period (days)

date=[1948 10 1; 1953 09 30]; % 5-years simulation
t_start=find(ismember(X01055500(:,1:3),date(1,:),'rows')>0);
t_end  =find(ismember(X01055500(:,1:3),date(2,:),'rows')>0);
time   =t_start:t_end;
prec = X01055500(time,4);
ept  = X01055500(time,5);
flow = X01055500(time,6);
temp = mean(X01055500(time,7:8),2);

X_labels = {'TS','CFMAX','CFR','CWH','BETA','LP','FC','PERC','K0','K1','K2','UZL','MAXBAS'} ; % uncertain parameters
M = length(X_labels) ;% number of uncertain parameters

% Parameter ranges (from Kollat et al.(2012)) 
xmin =[-3 0 0 0 0 0.3 1 0 0.05 0.01 0.05 0 1];
xmax =[3 20 1 0.8 7 1 2000 100 2 1 0.1 100 6];
DistrPar = cell(M,1); 
for i=1:M-1; DistrPar{i} = [ xmin(i) xmax(i) ] ; end
DistrPar{end} = xmax(end) ; % parameter MAXBAS is an integer
% Parameter distributions
DistrFun  = cell(1,M) ; 
for i=1:M-1;  DistrFun{i}='unif';  end; % uniform
DistrFun{end}='unid'; % discrete uniform

% Define output:
myfun = 'hbv_snow_objfun';

%% Step 3 (sample inputs space)

r = 100 ; % Number of samples

% % Option 1: use the sampling method originally proposed by Morris (1991).
% %
% % This requires specifying the number of levels in the uniform grid (L).
% % In this specific case, because one of the parameters (MAXBAS) is 
% % discrete-valued, the value of 'L' must coincide with the number of
% % feasible values for the discrete parameter 
% % (6 in this case, ranging from 1 to 6).
% L = 6  ;
% design_type  = 'trajectory'; % (note used here but required later)
% X = Morris_sampling(r,xmin,xmax,L); % (r*(M+1),M)

% Option 2: Latin Hypercube sampling strategy.
%
% In this case, we do not need to make any specific choice for the tuning
% parameters of the sampling strategy, however it may happen that during
% the sampling we will get several 'warning' messages due to the fact that
% the 'OAT_sampling' function checks that each consecutive sampled input
% vectors differ at least by one component, and if they don't, randomly
% change one of the two.
SampStrategy = 'lhs' ; % Latin Hypercube
design_type='radial'; % other option is 'trajectory'
X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);

%% From now on just the same as in 'workflow_eet_hymod.m'...

% Step 4 (run the model) 
Y = model_execution(myfun,X,prec,temp,ept,flow,warmup,Case) ;

% Step 5 (Computation of the Elementary effects)

% Choose one among multiple outputs for subsequent analysis:
Yi = Y(:,2); 

% Compute indices:
[ mi, sigma ] = EET_indices(r,xmin,xmax,X,Yi,design_type); 
% Plot:
EET_plot(mi, sigma,X_labels )

% Use bootstrapping to derive confidence bounds:
Nboot=100;
[mi, sigma, EE, mi_sd, sigma_sd, mi_lb, sigma_lb, mi_ub, sigma_ub, mi_all, sigma_all] = EET_indices(r,xmin,xmax,X,Yi,design_type,Nboot);
% Plot results in the plane:
EET_plot(mi,sigma,X_labels,mi_lb,mi_ub,sigma_lb,sigma_ub)

% Convergence analysis:
rr = [ r/10:r/10:r ] ;
[m_r,s_r,m_lb_r,m_ub_r,s_lb_r,s_ub_r,m_sd_r,s_sd_r]=EET_convergence(EE,rr);
% Plot:
figure
plot_convergence(m_r,rr*(M+1),[],[],[],'no of model evaluations','mean of EEs',X_labels)

% Convergence analysis with bootstrapping:
Nboot = 100;
rr = [ r/5:r/5:r ] ;
[m_r,s_r,m_lb_r,m_ub_r,s_lb_r,s_ub_r,m_sd_r,s_sd_r]=EET_convergence(EE,rr,Nboot);
% Plot:
figure
plot_convergence(m_r,rr*(M+1),m_lb_r,m_ub_r,[],'no of model evaluations','mean of EEs',X_labels)

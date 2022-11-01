% This script provides an application example of Regional Sensitivity
% Analysis (RSA)
%
% METHODS
%
% This script provides an example of application of Regional Sensitivity
% Analysis (RSA). RSA is applied in two different ways:
%
% - according to the original formulation
% where the input samples are split into two datasets depending on whether
% the corresponding output satisfies a threshold condition 
% >> matlab functions: RSA_indices_thres, RSA_plot_thres, RSA_convergence_thres
% (see help of 'RSA_indices_thres' for more details and references)
%
% - splitting the input samples into 'ngroup' datasets corresponding to 
% equally spaced ranges of the output
% >> matlab functions: RSA_indices_groups, RSA_plot_groups, RSA_convergence_groups
% (see help of 'RSA_indices_groups' for more details and references)
%
% MODEL AND STUDY AREA
%
% The model under study is the rainfall-runoff model Hymod
% (see help of function hymod_sim.m for more details) 
% applied to the Leaf catchment in Mississipi, USA
% (see header of file LeafCatch.txt for more details).
% The inputs subject to SA are the 5 model parameters, and the scalar 
% output for SA is (one or multiple) performance metric.
%
% INDEX
%
% Steps:
% 1. Add paths to required directories
% 2. Load data and set-up the Hymod model
% 3. Sample inputs space
% 4. Run the model against input samples
% 5a. Perform RSA with thresholds 
% or
% 5b. Perform RSA with groups

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

%% Step 2 (setup the Hymod model)

% Load data:
load -ascii LeafCatch.txt
rain = LeafCatch(1:365,1)   ;
evap = LeafCatch(1:365,2)   ;
flow = LeafCatch(1:365,3)   ;

% Define inputs:
DistrFun  = 'unif'  ; % Parameter distribution
DistrPar  = { [ 0 400 ]; [ 0 2 ]; [ 0 1 ]; [ 0 0.1 ] ; [ 0.1 1 ] } ; % Parameter ranges (from literature)
x_labels = {'Sm','beta','alfa','Rs','Rf'} ;

% Define output:
myfun = 'hymod_MulObj' ;
warmup = 30 ; % lenght of model warmup period (days)
% (sensitivity indices will not be computed for the warmup period)
%% Step 3 (sample inputs space)
SampStrategy = 'lhs' ; % Latin Hypercube
N = 3000 ; % Number of samples
M = length(DistrPar) ; % Number of inputs
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,N);

%% Step 4 (run the model) 
Y = model_execution(myfun,X,rain,evap,flow) ;

%% Step 5a (Regional Sensitivity Analysis with threshold)

% Visualize input/output samples (this may help finding a reasonable value
% for the output threshold):
figure; scatter_plots(X,Y(:,1),[],'rmse',x_labels) ;
figure; scatter_plots(X,Y(:,2),[],'bias',x_labels) ;

% Set output threshold:
rmse_thres = 3   ; %  threshold for the first obj. fun.
bias_thres = 0.5 ; % behavioural threshold for the second obj. fun.

% RSA (find behavioural parameterizations):
threshold = [ rmse_thres bias_thres ] ;
[mvd,idxb] = RSA_indices_thres(X,Y,threshold) ;

% Highlight the behavioural parameterizations in the scatter plots:
figure; scatter_plots(X,Y(:,1),[],'rmse',x_labels,idxb) ;
figure; scatter_plots(X,Y(:,2),[],'bias',x_labels,idxb) ;

% Plot parameter CDFs:
RSA_plot_thres(X,idxb,[],x_labels,{'behav','non-behav'}); % add legend

% Check the ranges of behavioural parameterizations by
% Parallel coordinate plot:
parcoor(X,x_labels,[],idxb)

% Plot the sensitivity indices (maximum vertical distance between
% parameters CDFs):
figure; boxplot1(mvd,x_labels,'mvd')

% Assess robustness by bootstrapping:
Nboot = 100;
[mvd,idxb,mvd_lb,mvd_ub] = RSA_indices_thres(X,Y,threshold,[],Nboot);
% Plot results:
boxplot1(mvd,x_labels,'mvd',mvd_lb,mvd_ub)

% Repeat computations using an increasing number of samples so to assess
% convergence:
NN = [ N/5:N/5:N ] ;
mvd = RSA_convergence_thres(X,Y,NN,threshold) ;
% Plot the sensitivity measures (maximum vertical distance between
% parameters CDFs) as a function of the number of samples:
figure; plot_convergence(mvd,NN,[],[],[],'no of samples','mvd',x_labels)

% Repeat convergence analysis using bootstrapping to derive confidence
% bounds:
Nboot = 100;
[mvd_mean,mvd_lb,mvd_ub]  = RSA_convergence_thres(X,Y,NN,threshold,[],Nboot) ;
figure
plot_convergence(mvd_mean,NN,mvd_lb,mvd_ub,[],'no of samples','mvd',x_labels)

%% Step 5b (Regional Sensitivity Analysis with groups)

% RSA (find behavioural parameterizations):
[stat,idx,Yk] = RSA_indices_groups(X,Y(:,1)) ;

% Plot parameter CDFs:
RSA_plot_groups(X,idx,Yk) ;
% Customize labels and legend:
RSA_plot_groups(X,idx,Yk,[],x_labels,'rmse'); % add legend


% This script provides an application example of the Generalized Likelihood 
% Uncertainty Estimation (GLUE) method (see help GLUE.m for more details
% about this Uncertainty Analysis method and references).
%
% MODEL AND STUDY AREA
%
% The model under study is the rainfall-runoff model Hymod
% (see help of function hymod_sim.m for more details) 
% applied to the Leaf catchment in Mississipi, USA
% (see header of file LeafCatch.txt for more details).
% The inputs subject to SA are the 5 model parameters, and the scalar 
% output for SA is (one or multiple) performance metric.

% This script prepared by Francesca Pianosi and Fanny Sarrazin
% University of Bristol, 2014
% mail to: francesca.pianosi@bristol.ac.uk

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

% Load data:
load -ascii LeafCatch.txt
rain = LeafCatch(1:365,1)   ;
evap = LeafCatch(1:365,2)   ;
flow = LeafCatch(1:365,3)   ;

% Plot data:
figure
subplot(311)
plot(rain)
subplot(312)
plot(evap)
subplot(313)
plot(flow)
xlabel('time (days)')

%%%% EXPERIMENT SETUP

M  = 5 ; % number of uncertain parameters [ Sm beta alfa Rs Rf ]
DistrFun  = 'unif'  ; % Parameter distribution
% Parameter ranges:
xmin = [   0 0 0 0   0.1 ] ; % minimum values
xmax = [ 400 2 1 0.1 1   ] ; % maximum values

%%%% SAMPLING INPUT SPACE

% Sample the input space using the 'AAT_sampling' function
% (see help AAT_sampling.m for more details)

addpath('sampling')

SampStrategy = 'lhs' ; % sampling strategy (other options is 'rsu')
N = 3000 ; % sample size
% Create data structure for parameter ranges
% as required by AAT_sampling
for i=1:M; DistrPar{i} = [xmin(i) xmax(i) ] ; end % this is to resort to
% the format required by AAT_sampling function.
% Perform sampling:
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,N);
% Plot results (for instance samples of param.1 vs param. 2):
figure; plot(X(:,1),X(:,2),'ob')

%%%% MODEL EVALUATION

% Define output function:
myfun = 'hymod_sim' ;
% (in the first place we just want to simulate the model and look at the
% time series of simulated flows)
warmup = 30 ; % lenght of model warmup period (days)
% (sensitivity indices will not be computed for the warmup period)

% Run the model and compute selected model output at sampled parameter
% sets:
Qsim = model_execution(myfun,X,rain,evap) ;
% Check size of 'Qsim':
size(Qsim)

% Plot simulated time series:
figure
plot(Qsim','b')
hold on; plot(flow,'ok') % add flow observations as black circles
xlabel('time (days)'); ylabel('flow (m3/s)')

%%%% GLUE

% Find 'behavioural' parameterizations and associated prediction limits 
% by the 'GLUE' function
% (see help GLUE.m for more details)

addpath('GLUE')

% Compute Root Mean Squared Error for each time series:
Y = RMSE(Qsim(:,warmup:end),flow(warmup:end)') ;
% Check size of 'Y':
size(Y)

% Plot the distribution of the output samples:
figure
plot_cdf(Y,'RMSE (m3/s)')

% Apply GLUE
threshold = 2.5 ; % threshold value
% (tip: make sure the sign of 'Y' and 'threshold' is consistent
% with the inequality assumption (above/below threshold) 
% of the GLUE function)
[ idx, Llim, Ulim ] = GLUE(Y,threshold,Qsim) ;

% Plot 'behavioural' time series:
figure
plot(Qsim','b')
hold on
plot(Qsim(idx,:)','r')
plot(flow,'ok')
xlabel('time (days)'); ylabel('flow (m3/s)')

% Plot prediction limits by GLUE:
figure
plot([Llim Ulim],'m','LineWidth',2)
hold on
plot(flow,'ok')
xlabel('time (days)'); ylabel('flow (m3/s)')

%%%% VISUALIZATION TOOLS

% Use some visualization functions to gain further insights
% (see help scatter_plots.m 
% and help parcoor.m)

% Scatter plots:
figure
scatter_plots(X,Y)
% customize the plot:
X_Labels = {'Sm','beta','alfa','Rs','Rf'} ;
figure
scatter_plots(X,Y,[],'RMSE',X_Labels)
% Highlight 'behavioural parameterizations' in different colour: 
figure
scatter_plots(X,Y,[],'RMSE',X_Labels,idx)

% Parallel coordinate plots:
parcoor(X,X_Labels,[],idx);


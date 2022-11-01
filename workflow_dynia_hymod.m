% This script provides an application example of the 
% dynamic identifiability analysis (DYNIA) by Wagener et at. (2003).
%
% MODEL AND STUDY AREA
%
% The model under study is the rainfall-runoff model Hymod
% (see help of function hymod_sim.m for more details) 
% applied to the Leaf catchment in Mississipi, USA
% (see header of file LeafCatch.txt for more details).
% The inputs subject to SA are the 5 model parameters, 
% and the time-varying model output is a metric of model performance
% computed over a moving window around each time step of the simulation.
%
% REFERENCES
%
% Wagener, T., McIntyre, N., Lees, M., Wheater, H., Gupta, H., 2003. 
% Towards reduced uncertainty in conceptual rainfall-runoff modelling: 
% dynamic identifiability analysis. Hydrol. Process. 17, 455-476.

% This script prepared by Francesca Pianosi
% University of Bristol, 2015
% mail to: francesca.pianosi@bristol.ac.uk


%% Step 1: Set current directory to 'my_dir' and add path to sub-folders:
my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

cd(my_dir)
addpath(genpath(my_dir))

%% Step 2 (setup the Hymod model)

% Load data:
load -ascii LeafCatch.txt
T = 2*365 ; % length of simulation horizon
rain = LeafCatch(1:T,1)   ;
evap = LeafCatch(1:T,2)   ;
flow = LeafCatch(1:T,3)   ;

% Define inputs:
DistrFun  = 'unif'  ; % Parameter distribution
DistrPar  = { [ 0 400 ]; [ 0 2 ]; [ 0 1 ]; [ 0 0.1 ] ; [ 0.1 1 ] } ; % Parameter ranges (from literature)
x_labels = {'Sm','beta','alfa','Rs','Rf'} ;
M = length(DistrPar) ; % number of SA inputs (parameters)

%% Step 3 (Define the model output and GSA method)

% Choose the function for model simulation:
simfun = 'hymod_sim' ;
warmup = 30 ; % lenght of model warmup period (days)
% (sensitivity indices will not be computed for the warmup period)

% Choose the objective function among the 3 available options: 
%obj_fun = 'MAE' ; % mean( | obs-sim | ) 
obj_fun = 'RMSE'; % mean( (obs-sim).^2 ) 
%obj_fun ='BIAS'; % | mean(obs) - mean(sim) |

% Define semi-length of the window for averaging performances (days):
W = 10 ; 
% (in other words, the objective function will be computed using 2W+1 days)

% Define percentage of samples that will be retained 
% to estimate the a posteriori frequency distribution of the parameters:
perc = 10 ; % % of the input/output sample size

% Define number of bins to compute the a posteriori frequency distribution:
nbins = 10 ; 


%% Step 3 Sample inputs space and evaluate the model

% Choose number of samples:
N = 1000 ; 

% Choose sampling strategy:
SampStrategy = 'lhs';

% Perform sampling:
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,N); % matrix (N,M)

% Evaluate model output:
Q_sim = model_execution(simfun,X,rain,evap,warmup);  % matrix (N,T) holding the 
% time series of simulated flows
Y = tv_objfun(Q_sim,flow',obj_fun,W) ; % matrix (N,T) holding the time
% series of time-varying objective functions

% Plot temporal pattern of objective function:
figure
clrs = autumn ;
clrs = clrs(end:-1:1,:);
colormap(clrs)
imagesc(Y)
xlabel('time')
ylabel('model run')
colorbar
title(['Obj.Fun.:' obj_fun ])
% add flow:
hold on 
Nplot = size(Y,1); % scaling factor
plot(Nplot+0.5-flow/max(flow)*Nplot,'k','LineWidth',2)

%% Step 4 DYNIA


for i=1:M
    
    % Compute posterior probability:
    [ xi, fi ] = dynia_histograms(X(:,i),Y,perc,nbins) ;
    
    % Plot results:
    dynia_plot(xi,fi,x_labels{i},flow) ;
    title(['Sensitivity indices - Obj.Fun: ' obj_fun ])

end 




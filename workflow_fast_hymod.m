% This script provides an application example of the 
% Fourier Amplitude Sensitivity Test (FAST).
% FAST uses the Fourier decomposition of the model output
% to approximate the variance-based first-order sensitivity indices.
% See help of 'FAST_indices.m' for more details and references.
%
% In this workflow, FAST is applied to the rainfall-runoff Hymod model
% (see help of 'hymod_sim.m' for more details) 
% applied to the Leaf catchment in Mississipi, USA
% (see header of file LeafCatch.txt for more details).
% The inputs subject to SA are the 5 model parameters, and the scalar 
% output for SA is a performance metric.
% 
% FAST estimates are compared to those obtained by the 'conventional'
% resampling approach used in Variance-Based SA
% [see help of 'vbsa_indices.m'].

% This script prepared by Francesca Pianosi and Fanny Sarrazin
% University of Bristol, 2014
% mail to: francesca.pianosi@bristol.ac.uk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: set paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: setup the model and define input ranges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data:
addpath([ my_dir '/example/hymod'])
load -ascii LeafCatch.txt
rain = LeafCatch(1:365,1)   ;
evap = LeafCatch(1:365,2)   ;
flow = LeafCatch(1:365,3)   ;

% Number of uncertain parameters subject to SA:
M    = 5 ; 
% Parameter ranges (from literature):
xmin = [  0 0 0 0   0.1 ];
xmax = [400 2 1 0.1 1   ];
% Parameter distributions:
DistrFun  = 'unif'  ; 
DistrPar = cell(M,1); 
for i=1:M; DistrPar{i} = [ xmin(i) xmax(i) ] ; end
% Name of parameters (will be used to costumize plots):
X_labels = {'Sm','beta','alfa','Rs','Rf'} ; 
% Define output:
myfun = 'hymod_nse' ;
warmup = 30 ; % lenght of model warmup period (days)
% (sensitivity indices will not be computed for the warmup period)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: approximate first-order sensitivity indices by FAST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FAST sampling:
[X,s] = FAST_sampling(DistrFun,DistrPar,M);

% Run the model and compute model output at sampled parameter sets:
Y = model_execution(myfun,X,rain,evap,flow,warmup) ;

% Estimate indices: 
Si_fast = FAST_indices(Y,M) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Convergence analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The 'FAST_sampling' function used above automatically set the sample size
% to the minimum possible given the number M of inputs (see help of
% FAST_sampling to learn about this). In our case this is:
N_fast = length(Y) ;

% We can now assess whether FAST estimates would change if using a larger
% number of samples:
NNfast = [N_fast:1000:N_fast+5000] ;
Si_fast_conv = nan(length(NNfast),M) ;
Si_fast_conv(1,:) = Si_fast ;
for n=2:length(NNfast)
    [Xn,sn] = FAST_sampling(DistrFun,DistrPar,M,NNfast(n));
    Yn = model_execution(myfun,Xn,rain,evap,flow,warmup) ;
    Si_fast_conv(n,:) = FAST_indices(Yn,M) ;
end
% Plot results:
figure
plot_convergence(Si_fast_conv,NNfast,[],[],[],'model evals','1st-order sensitivity',X_labels)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Comparison with VBSA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we compare FAST estimates with those obtained by the 'conventional'
% resampling approach used in Variance-Based SA
% [see help of 'vbsa_indices.m'].

addpath([ my_dir '/vbsa'])

% Set the base sample size for VBSA in such a way that
% the total number of model evaluations be the same as FAST:
Nvbsa = ceil(max(NNfast)/(M+2)) ;
% VBSA sampling:
SampStrategy = 'lhs' ;
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,2*Nvbsa);
[ XA, XB, XC ] = vbsa_resampling(X) ;
% Run the model and compute model output at sampled parameter sets:
YA = model_execution(myfun,XA,rain,evap,flow,warmup) ;
YB = model_execution(myfun,XB,rain,evap,flow,warmup) ;
YC = model_execution(myfun,XC,rain,evap,flow,warmup) ;
% Use output samples to estimate the indices at different sub-sample sizes:
NNvbsa = floor(NNfast/(M+2)) ;
Si_vbsa_conv = vbsa_convergence([YA;YB;YC],M,NNvbsa);

% Compare
xm = min(NNfast(1),NNvbsa(1)*(M+2)); 
xM = max(NNfast(end),NNvbsa(end)*(M+2)) ; 
ym = -0.1 ; 
yM = 1.1 ;
figure
subplot(211); plot_convergence(Si_fast_conv,NNfast      ,[],[],[],'model evals','FAST',X_labels)
axis([xm,xM,ym,yM])
subplot(212); plot_convergence(Si_vbsa_conv,NNvbsa*(M+2),[],[],[],'model evals','VBSA',X_labels); 
axis([xm,xM,ym,yM])


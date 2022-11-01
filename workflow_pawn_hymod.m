% This script provides an application example
% of the PAWN sensitivity analysis approach (Pianosi and Wagener, 2015)
%
% MODEL AND STUDY AREA
%
% The model under study is the rainfall-runoff model Hymod
% (see help of function hymod_sim.m for more details) 
% applied to the Leaf catchment in Mississipi, USA
% (see header of file LeafCatch.txt for more details).
% The inputs subject to SA are the 5 model parameters, and the scalar 
% output for SA is a statistic of the simulated time series
% (e.g. the maximum flow over the simulation horizon)
% 
% REFERENCES
%
% Pianosi, F. and Wagener, T. (2015), A simple and efficient method 
% for global sensitivity analysis based on cumulative distribution 
% functions, Env. Mod. & Soft., 67, 1-11.

% This script prepared by Francesca Pianosi and Fanny Sarrazin
% University of Bristol, 2015
% mail to: francesca.pianosi@bristol.ac.uk

%% Step 1: Add paths

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%% Step 2: setup the Hymod model
load -ascii LeafCatch.txt
rain = LeafCatch(1:730,1)   ;
evap = LeafCatch(1:730,2)   ;
flow = LeafCatch(1:730,3)   ;
warmup = 30; % model warmup period days

% Define uncertain inputs (parameters):
M = 5 ; % number of inputs
labelparams={ 'SM','beta','alfa','Rs','Rf' } ; % input names
% parameter ranges:
xmin = [  0 0 0 0   0.1 ];
xmax = [400 2 1 0.1 1   ];
distrpar=cell(M,1); for i=1:M; distrpar{i}=[xmin(i) xmax(i)]; end

% Define model output:
fun_test = 'hymod_max';

%% Step 3: Apply PAWN

NU = 150 ; % number of samples to estimate unconditional CDF
NC = 100 ; % number of samples to estimate conditional CDFs
n  = 10 ; % number of conditioning points

% Create input/output samples to estimate the unconditional output CDF:
Xu = AAT_sampling('lhs',M,'unif',distrpar,NU); % matrix (NU,M)
Yu = model_execution(fun_test,Xu,rain,evap,warmup)  ; % vector (1,M)

% Create input/output samples to estimate the conditional output CDFs:
[ XX, xc ] = pawn_sampling('lhs',M,'unif',distrpar,n,NC);
YY = pawn_model_execution(fun_test,XX,rain,evap,warmup) ;

% Estimate unconditional and conditional CDFs:
[ YF, Fu, Fc  ] = pawn_cdfs(Yu,YY) ;

% Plot CDFs:
figure
for i=1:M
   subplot(1,M,i)
   pawn_plot_cdf(YF, Fu, Fc(i,:),[],'y (max flow)')
end

% Further analyze CDF of one input:
i = 3 ;
figure;
pawn_plot_cdf(YF, Fu, Fc(i,:),xc{i},'y (max flow)',labelparams{i}) % same
% function as before but exploiting more optional input arguments

% Compute KS statistics:
KS = pawn_ks(YF,Fu,Fc) ;

% Plot KS statistics:
figure
for i=1:M
   subplot(1,M,i)
   pawn_plot_kstest(KS(:,i),NC,NU,0.05,xc{i},labelparams{i})
end

% Compute PAWN index by taking a statistic of KSs (e.g. max):
Pi = max(KS);

% Plot:
figure 
boxplot1(Pi,labelparams)

% Use bootstrapping to assess robustness of PAWN indices:
stat = 'max' ; % statistic to be applied to KSs
Nboot = 100  ; % number of boostrap resamples
[ T_m, T_lb, T_ub ] = pawn_indices(Yu,YY,stat,[],Nboot);

% Plot:
figure; boxplot1(T_m,labelparams,[],T_lb,T_ub)

% Convergence analysis:
stat = 'max' ; % statistic to be applied to KSs
NCb = [ NC/10 NC/2 NC ] ;
NUb = [ NU/10 NU/2 NU ] ;

[ T_m_n, T_lb_n, T_ub_n ] = pawn_convergence( Yu, YY, stat, NUb, NCb,[],Nboot );
NN = NUb+n*NCb ;
figure; plot_convergence(T_m_n,NN,T_lb_n,T_ub_n,[],'no of evals',[],labelparams)

%% Step 4: Apply PAWN to sub-region of the output range

% Compute the PAWN index over a sub-range of the output distribution, for
% instance only output values above a given threshold
thres = 50 ;
[ T_m2, T_lb2, T_ub2 ]= pawn_indices( Yu, YY, stat,[], Nboot,[],'above',thres ) ;

% Plot:
figure; boxplot1(T_m2,labelparams,[],T_lb2,T_ub2)



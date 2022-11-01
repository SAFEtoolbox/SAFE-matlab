% This script applies Variance-Based Sensitivity Analysis to the
% Ishigami-Homma function.
% This function is commonly used to test approximation procedures
% of variance-based indices because its output variance, first-order and 
% total-order indices (or 'main' and 'total' effects) can be analytically
% computed. Therefore this script mainly aims at analyzing the accuracy
% and convergence of the function 'vbsa_indices'. 

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

% Setup the model and define input ranges

fun_test  = 'ishigami_homma_function' ;
M = 3 ;
distr_fun = 'unif' ;
distrpar =  [ -pi pi ];

% Compute the exact values of the output variance (V) and of 
% the first-order (Si_ex) and total-order (STi_ex) 
% variance-based sensitivity indices (this is possible in this
% very specific case because V, Si_ex and STi_ex can be computed
% analytically)
[ tmp, V, Si_ex, STi_ex ] = ishigami_homma_function(rand(1,M)) ;

% Sample parameter space:
SampStrategy = 'lhs' ;
N = 3000 ; 
X = AAT_sampling(SampStrategy,M,distr_fun,distrpar,2*N);
% Apply resampling strategy for the efficient approximation of the indices:
[ XA, XB, XC ] = vbsa_resampling(X) ;

% Run the model and compute selected model output at sampled parameter
% sets:
YA = model_execution(fun_test,XA) ; % size (N,1)
YB = model_execution(fun_test,XB) ; % size (N,1)
YC = model_execution(fun_test,XC) ; % size (N*M,1)

% Compute main (first-order) and total effects:
[ Si, STi ] = vbsa_indices(YA,YB,YC);

% Plot main and total effects and compare the values estimated 
% by the function 'vbsa_indices' with the exact values
figure 
subplot(121); boxplot2([ Si; Si_ex])
subplot(122); boxplot2([STi; STi_ex])
legend('estimated','exact')% add legend

% Analyze convergence of sensitivity indices:
NN = [N/5:N/5:N] ;
[ Si, STi ] = vbsa_convergence([YA;YB;YC],M,NN);
figure
subplot(121); plot_convergence(Si,NN*(M+2),[],[],Si_ex,'model evals','main effect')
subplot(122); plot_convergence(STi,NN*(M+2),[],[],STi_ex,'model evals','total effect')


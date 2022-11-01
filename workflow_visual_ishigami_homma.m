% This script provides an application example of how to use several
% visualization tools (scatter plots, coloured scatter plots, parallel
% coordinate plots, Andres' plot) to learn about sensitivity.
% The application example is the Hishigami-Homma function, which is a
% standard benchmark function in the Sensitivity Analysis literature.
% (see help of function 'ishigami_homma_function.m' for more details and
% references).

% This script prepared by Francesca Pianosi and Fanny Sarrazin
% University of Bristol, 2014
% mail to: francesca.pianosi@bristol.ac.uk

%% Step 1 (set directories)

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%% Step 2 (setup the model)
fun_test  = 'ishigami_homma_function' ;
M = 3 ;
distr_fun = 'unif' ;
distrpar =  [ -pi pi ];

%% Step 3 (sampling and model evaluation)

N = 3000;
X = AAT_sampling('lhs',M,'unif',distrpar,N);
Y = model_execution(fun_test,X)          ;

%% Step 4 (Scatter plots)

% Use scatter plots of inputs againts output to visually assess 
% direct effects:
scatter_plots(X,Y)

% Use coloured scatter plots of one input against another on to assess
% interactions:
i1 = 1 ;
i2 = 2 ;
figure; scatter_plots_col(X,Y,i1,i2) % plot x(i1) against x(i2)
% Change i2:
i2 = 3 ;
scatter_plots_col(X,Y,i1,i2)
% Put all possible combinations of i1,i2 into one figure:
scatter_plots_interaction(X,Y)
% Customize titles:
scatter_plots_interaction(X,Y,[],{'x(1)','x(2)','x(3)'})

%% Step 5 (Parallel Coordinate Plot)

% Use Parallel Coordinate Plot to find patterns of input combinations 
% mapping into specific output condition
idx = Y>30 ; % output condition to be highlighted
parcoor(X,{'x(1)','x(2)','x(3)'},[],idx)

%% Step 6 (Andres' visualization test)

Xref= [ 2 2 2 ] ;
idx = 2 ;
Andres_plots(X,Y,Xref,idx,fun_test) 



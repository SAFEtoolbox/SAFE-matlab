% This script provides an application example of time-varying 
% sensitivity analysis. 
% The user can choose one of three GSA methods (FAST, VBSA or EET) 
% and this script will guide through the application of that method
% to a time-varying model output.
% In this example, the time-varying model output is the model prediction
% (estimated flow) at each time step of the simulation.
%
% MODEL AND STUDY AREA
%
% The model under study is the rainfall-runoff model Hymod
% (see help of function hymod_sim.m for more details) 
% applied to the Leaf catchment in Mississipi, USA
% (see header of file LeafCatch.txt for more details).
% The inputs subject to SA are the 5 model parameters, 
% and the time-varying model output is the model prediction
% (estimated flow) at each time step of the simulation.
% For an example of how to use this type of analysis, 
% see for instance Reusser and Zehe (2011)
%
% REFERENCES
%
% Reusser, D. E., Zehe, E., 2011. Inferring model structural deficits 
% by analyzing temporal dynamics of model performance and parameter 
% sensitivity. Water Resources Research 47 (7).

% This script prepared by Francesca Pianosi
% University of Bristol, 2015
% mail to: francesca.pianosi@bristol.ac.uk

%% Step 1: Set current directory to 'my_dir' and add path to sub-folders:
my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%% Step 2 Setup the Hymod model and define input variability space

% Load data:
load -ascii LeafCatch.txt
rain = LeafCatch(1:730,1)   ;
evap = LeafCatch(1:730,2)   ;
flow = LeafCatch(1:730,3)   ;

% Define inputs:
DistrFun  = 'unif'  ; % Parameter distribution
% Parameter ranges (from literature):
xmin = [   0 0 0 0   0.1 ] ;
xmax = [ 400 2 1 0.1 1   ] ;
% Parameter names (needed for plots):
x_labels = {'Sm','beta','alfa','Rs','Rf'} ;

%% Step 3 Define the model output and GSA method

myfun = 'hymod_sim' ;
warmup = 30 ; % lenght of model warmup period (days)
% (sensitivity indices will not be computed for the warmup period)

GSA_met = 'FAST' ;
%GSA_met = 'EET' ;
%GSA_met = 'VBSA' ;

%% Step 3 Define choices for sampling

r = 20   ; % number of EEs (needed for EET only)

N = 3001 ; % sample size (needed for all other methods) 
% Notes:
% Case VBSA: N is the size of the base sample: the total number of model 
% evaluations will be N*(M+2) [more on this: see help vbsa_resampling]
% Case FAST: N must be odd [more on this: see help FAST_sampling]
% Case EET: N is not used, total number of model evaluations depend on 'r'
% and precisely it will be r*(M+1) [more on this: see help OAT_sampling]

SampStrategy = 'lhs' ; % Sampling strategy (needed for EET and VBSA)

%% Step 4 Sample inputs space and evaluate the model


M = length(xmin) ; % Number of inputs
DistrPar = cell(M,1) ; for i=1:M; DistrPar{i} = [ xmin(i) xmax(i) ] ; end
 

if strcmp(GSA_met,'FAST')
    
    X = FAST_sampling(DistrFun,DistrPar,M,N) ;
    Y = model_execution(myfun,X,rain,evap,warmup); % (N,T) 


elseif strcmp(GSA_met,'VBSA')
    
    X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,2*N);
    [ XA, XB, XC ] = vbsa_resampling(X) ;
    YA = model_execution(myfun,XA,rain,evap,warmup) ; % size (N,T)
    YB = model_execution(myfun,XB,rain,evap,warmup) ; % size (N,T)
    YC = model_execution(myfun,XC,rain,evap,warmup) ; % size (N*M,T)
    Y = [ YA; YB; YC ] ; % (N*(2+M),T) ... only needed for the next plot!

elseif strcmp(GSA_met,'EET')

    design_type = 'radial';
    X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);
    % More options are actually available for EET sampling,
    % see workflow_eet_hymod
    Y = model_execution(myfun,X,rain,evap,warmup); % (r*(M+1),T)
    
else
    error('No method called %s',GSA_met)
    
end

% Plot results:
figure
plot(Y','b')
xlabel('time')
ylabel('flow')
axis([1,size(Y,2),0,max(max(Y))])

%% Step 5 Compute time-varying Sensitivity Indices

T = size(Y,2) ;

if strcmp(GSA_met,'FAST')
    
    Si = nan(T,M) ;
    for t=warmup:T
        Si(t,:) = FAST_indices(Y(:,t),M) ;
    end
    S_plot = Si' ;
    
elseif strcmp(GSA_met,'VBSA')
    
    Si  = nan(T,M) ;
    STi = nan(T,M) ;
    for t=warmup:T
        [ Si(t,:), STi(t,:) ] = vbsa_indices(YA(:,t),YB(:,t),YC(:,t));
    end  
    % select sensitivity index to be plotted in the next figure:
    S_plot = Si' ; 
    %S_plot = STi' ;

elseif strcmp(GSA_met,'EET')
    
    mi    = nan(T,M) ;
    sigma = nan(T,M) ;
    for t=warmup:T
        [ mi(t,:), sigma(t,:) ] = EET_indices(r,xmin,xmax,X,Y(:,t),design_type) ;
    end  
    % select sensitivity index to be plotted in the next figure:
    S_plot = mi' ; 
    %S_plot = sigma' ;
    
end

% Plot results:

figure
clrs = autumn ;
clrs = clrs(end:-1:1,:);
colormap(clrs)
imagesc(S_plot)
xlabel('time')
ylabel('inputs')
set(gca,'YTick',1:M,'YTickLabel',x_labels)
colorbar
title(['Sensitivity indices - GSA meth: ' GSA_met ])
hold on
plot(M+0.5-flow/max(flow)*M,'k','LineWidth',2)




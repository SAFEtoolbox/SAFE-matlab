% This script provides an application example of the 
% Fourier Amplitude Sensitivity Test (FAST).
% FAST uses the Fourier decomposition of the model output
% to approximates the variance-based first-order sensitivity indices.
% See help of 'FAST_indices.m' for more details and references.
%
% In this workflow, FAST is applied to the Sobol g-function
% with varying number of inputs and different parameterizations for the
% function parameters not subject to SA
% [see help of 'sobol_g_function.m'].
% 
% FAST estimates are compared to those obtained by the 'conventional'
% resampling approach used in Variance-Based SA
% [see help of 'vbsa_indices.m'].
%
% This script can be used to reproduce Figure 2 in:
% Saltelli and Bolado (1998), An alternative way to compute Fourier
% amplitude sensitivity test (FAST), Computational Statistics & Data
% Analysis, 26, 445-460.

% This script prepared by Francesca Pianosi and Fanny Sarrazin
% University of Bristol, 2014
% mail to: francesca.pianosi@bristol.ac.uk

%% Step 1: set paths

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

% Define output function:
myfun  = 'sobol_g_function' ;
% Define input distribution and ranges:
DistrFun = 'unif' ;
DistrPar =  [ 0 1 ];
    
aa = [0 1 9 99] ; % options for the (fixed) parameters
MM = 5:11 ; % options for the number of inputs

ms = 16;
figure; hold on

for param = 1:4
        
    SM = nan(length(MM),3) ;
    Nfast = nan(1,length(MM)) ;
    Nvbsa = nan(1,length(MM)) ;
    for m=1:length(MM)
        
        M = MM(m);
        a = ones(1,M)*aa(param) ;
        % Analytic:
        [ tmp, V_ex, Si_ex ] = sobol_g_function(rand(1,M),a) ;                
        % FAST:
        [X,s] = FAST_sampling(DistrFun,DistrPar,M);
        % Run the model and compute selected model output at sampled parameter
        % sets:
        Y = model_execution(myfun,X,a) ; % size (Ns,1)
        [Si_fast,V_fast,A,B] = FAST_indices(Y,M) ;       
        % VBSA:
        SampStrategy = 'lhs' ;
        Nfast(m) = length(Y) ;        
        % Option 1: set the base sample size for VBSA in such a way that
        % the total number of model evaluations be the same as FAST:
        Nvbsa(m) = ceil(Nfast(m)/(M+2)) ;
        % Option 2: set the base sample size to 4096 independently by the
        % sample size used for FAST (as done in Saltelli and Bolado, 1998)
        % [THIS OPTION MAY BE TIME-CONSUMING!]
   %     Nvbsa(m) = 4096 ;
        X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,2*Nvbsa(m));
        [ XA, XB, XC ] = vbsa_resampling(X) ;
        YA = model_execution(myfun,XA,a) ; % size (N,1)
        YB = model_execution(myfun,XB,a) ; % size (N,1)
        YC = model_execution(myfun,XC,a) ; % size (N*M,1)
        Si_vbsa = vbsa_indices(YA,YB,YC);    
        SM(m,:) = [ mean(Si_ex) mean(Si_vbsa) mean(Si_fast) ] ;

    end
    
    subplot(2,2,param); hold on
    % next three lines just for the legend:
    plot(MM(1),SM(1,1),'sk','MarkerSize',ms,'MarkerFaceColor',[126 126 126]/256); hold on
    plot(MM(1),SM(1,2),'^k','MarkerSize',ms,'LineWidth',2)
    plot(MM(1),SM(1,3),'ok','MarkerSize',ms,'LineWidth',2)
    % next three for plotting:
    plot(MM,SM(:,1),'sk','MarkerSize',ms,'MarkerFaceColor',[126 126 126]/256)
    plot(MM,SM(:,2),'^k','MarkerSize',ms,'LineWidth',2)
    plot(MM,SM(:,3),'ok','MarkerSize',ms,'LineWidth',2)
    % customize picture:
    set(gca,'XTick',MM,'XLim',[MM(1)-1,MM(end)+1],'FontSize',20)
    set(gca,'YLim',[-0.1,1.1])
    box on
    xlabel('Number of inputs M','FontSize',20)
    ylabel('1st-order sensitivity','FontSize',20)
    legend('Analytic','VBSA','FAST')
    title(['a(i)=' num2str(aa(param)) ])
    
end

% Plot number of model evaluations against number of inputs:
figure; hold on
plot(MM,Nvbsa.*(MM+2),'^k','MarkerSize',ms,'LineWidth',2)
plot(MM,Nfast,'ok','MarkerSize',ms,'LineWidth',2)
set(gca,'XTick',MM,'FontSize',20)
xlabel('Number of inputs M','FontSize',20)
ylabel('Number of model evaluations N','FontSize',20)
legend('VBSA','FAST')
box on
grid on











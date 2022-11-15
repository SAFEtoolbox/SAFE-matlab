% This workflow scripts demonstrates the approximation strategy 
% of PAWN indices from generic input-output dataset
% by application to the Ishigami-homma function, 
% as described in Section 3.1 of:
%
% Pianosi, F. and Wagener, T. (2018), Distribution-based sensitivity 
% analysis from a generic input-output sample, Env. Mod. & Soft.
%
% The code works together with the SAFE Toolbox. More information and
% dowload at: www.bris.ac.uk/cabot/resources/safe-toolbox/
%
% This code developed by Francesca Pianosi, University of Bristol
% email: francesca.pianosi@bristol.ac.uk

%%%%%%%%%%%%%%%%%%%%% SET PATH %%%%%%%%%%%%%%%%%%%%%

% Set the path to the SAFE Toolbox directory
% (where you must have also copied the additional functions/data 
% for applying PAWN to generic datasets)

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.1';

addpath(genpath(my_dir))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute variance-based sensitivity indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 3 ;
fun_test  = 'ishigami_homma_function' ;
DistrFun = 'unif' ;
DistrPar =  [ -pi pi ];
[ y, V, Si_ex, STi_ex ] = ishigami_homma_function(rand(1,M));
Si_ex
STi_ex
% what do we infer from these in terms of input ranking?
% x1 most influential factor; 
% x2 has (little) influence but no interactions
% x3 only influential through interactions with x1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute distribution-based (PAWN) sensitivity indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tailored sampling (Pianosi and Wagener, 2015):
Nu = 150; Nc = 50; n =3;
Xu = AAT_sampling('lhs',M,DistrFun,DistrPar,Nu); % matrix (NU,M)
Yu = model_evaluation(fun_test,Xu)  ; % vector (1,M)
[ XX, xc ] = pawn_sampling('lhs',M,DistrFun,DistrPar,n,Nc); % XX= (M,n)
YY = pawn_model_evaluation(fun_test,XX)  ; % vector (1,M)

% Generic sampling + splitting (Pianosi and Wagener, 2018):
N = 150; n = 3;
X = AAT_sampling('lhs',M,DistrFun,DistrPar,N);
Y = model_evaluation(fun_test,X)  ;
[KS,KS_dummy,YF,Fu,Fc,YY2,xc2,NC,XX2,YU,idx_bootstrap] = pawn_ks_givendata(X,Y,n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +++ FIGURE 1: compare sampling approaches ++++
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i =1; % input factor under investigation
fs = 22; ms = 10 ; % font size and marker size

% ^^^^ tailored sampling ^^^

figure; 
% (a): unconditional input samples (tailored sampling)
subplot(1,2,1); hold on
plot3(Xu(:,1),Xu(:,2),Xu(:,3),'ok','MarkerFaceColor','r','MarkerSize',ms); box on
set(gca,'FontSize',fs,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{})
xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); 
axis([-pi,pi,-pi,pi,-pi,pi])

% (b): conditional input samples (tailored sampling)
col=gray(n+1); % Color for colormap, (n+1) so that the last line is not white
map = colormap(col(1:n,:));% Remove white color
xci = xc{i}; [tmp,idx] = sort(xci);
subplot(1,2,2)
for k=1:n
    plot3(XX{i,idx(k)}(:,1),XX{i,idx(k)}(:,2),XX{i,idx(k)}(:,3),'ok','MarkerFaceColor',map(k,:),'MarkerSize',ms); hold on; box on
end
set(gca,'FontSize',fs,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{})
xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); 
axis([-pi,pi,-pi,pi,-pi,pi])

% (c): scatter plot:
figure; hold on; 
plot(Xu(:,i),Yu,'ok','MarkerFaceColor','r','MarkerSize',ms); hold on
for k=1:n
    plot(XX{i,idx(k)}(:,i),YY{i,idx(k)},'ok','MarkerFaceColor',map(k,:),'MarkerSize',10); hold on; box on
end
set(gca,'FontSize',fs,'XTick',sort(xc{i}),'XTickLabel',{},'YTickLabel',{},'XLim',[-pi,pi],'XGrid','on')
xlabel('x_1'); ylabel('y')

% ^^^^ generic sampling ^^^

figure; hold on; 
% (d): unconditional input samples (generic sampling)
subplot(1,2,1)
plot3(X(:,1),X(:,2),X(:,3),'ok','MarkerSize',10); hold on; box on
plot3(X(idx_bootstrap,1),X(idx_bootstrap,2),X(idx_bootstrap,3),'ok','MarkerFaceColor','r','MarkerSize',10)
set(gca,'FontSize',fs,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{})
xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); 
axis([-pi,pi,-pi,pi,-pi,pi])

% (e): conditional input samples (generic sampling)
subplot(1,2,2)
xtick2 = [];
for k=1:n
    plot3(XX2{i,k}(:,1),XX2{i,k}(:,2),XX2{i,k}(:,3),'ok','MarkerFaceColor',map(k,:),'MarkerSize',10); hold on; box on
    xtick2 = [ xtick2, min(XX2{i,k}(:,i)), max(XX2{i,k}(:,i)) ] ;
end
set(gca,'FontSize',fs,'XTick',xtick2,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{},'XGrid','on')
xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); 
axis([-pi,pi,-pi,pi,-pi,pi])

% (f): scatter plot:
figure; hold on; 
xtick2 = [];
for k=1:n
    plot(XX2{i,k}(:,i),YY2{i,k},'ok','MarkerFaceColor',map(k,:),'MarkerSize',10); hold on; box on
    xtick2 = [ xtick2, min(XX2{i,k}(:,i)), max(XX2{i,k}(:,i)) ] ;
end
plot(X(idx_bootstrap,i),Y(idx_bootstrap),'ok','MarkerFaceColor','r','MarkerSize',ms); hold on
set(gca,'FontSize',fs,'XTick',xtick2,'XTickLabel',{},'YTickLabel',{},'XLim',[-pi,pi],'XGrid','on')
xlabel('x_1'); ylabel('y')

% ^^^^ CDFs and KSs ^^^

figure; hold on;
% (g): plot CDFs
subplot(1,2,1)
for k=1:n
    plot(YF,Fc{i,k},'Color',map(k,:),'LineWidth',3.5); hold on
end
plot(YF,Fu,'r','LineWidth',3.5)
set(gca,'FontSize',fs); xlabel('y'); ylabel('cdf'); 
% (h): plot KS
subplot(1,2,2)
for k=1:n
    plot(xc2{i}(k),KS(k,i),'ok','MarkerFaceColor',map(k,:),'MarkerSize',ms); hold on
end
set(gca,'FontSize',fs,'XTick',round(xc2{i},1),'YLim',[0,1],'XLim',[-pi,pi],'XGrid','on')
xlabel('x_1'); ylabel('KS'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +++ FIGURE 3: IH results for varying (N,n) +++
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NN = [ 100 500 2000 ] ; % sample size to be tested
nn = [ 1:14 ] ; % number of conditioning intervals to be tested
alfa  = 0.05 ; % confidence level for intervals around sensitivity indices

% For each combination of N and n,
% compute the PAWN indices (using stat=median)
% and associated confidence intervals 
% as well as KS of the 'dummy parameter'
% (this may take some computing time):
PAWN_median    = cell(length(NN),1);
PAWN_median_lb = cell(length(NN),1);
PAWN_median_ub = cell(length(NN),1);
KS_dummy_mean  = nan(length(NN),length(nn)) ;
for jj= 1:length(NN)    
    PAWN_median{jj} = nan(length(nn),M) ;
    for k = 1:length(nn)
        X = AAT_sampling('lhs',M,DistrFun,DistrPar,NN(jj));
        Y = model_evaluation(fun_test,X)  ;
        nboot = 50;
        [KS_median,KS_mean,KS_max,KS_dummy] = pawn_indices_givendata(X,Y,nn(k),nboot) ;
        median_lb = sort(KS_median) ; median_lb = median_lb(max(1,round(nboot*alfa/2)),:);
        median_ub = sort(KS_median) ; median_ub = median_ub(round(nboot*(1-alfa/2)),:);
        PAWN_median_lb{jj}(k,:) = median_lb ;
        PAWN_median_ub{jj}(k,:) = median_ub ;        
        PAWN_median{jj}(k,:) = mean(KS_median); % (n,M)
        KS_dummy_mean(jj,k) = mean(KS_dummy);
    end
end

% Plot results:
hfig = figure; hold on
ms = 10; fs = 24; % marker and font size
legend_NN = cell(length(NN),1); for jj=1:length(NN); legend_NN{jj} = [ 'N=' num2str(NN(jj)) ] ; end % data for the legend
pert = [-0.2,0,0.2] ; col = {'b','m','g'} ; % positioning and color line
old_res = [0.48,0.14,0.3]; n_old_res = 15 ; % values of PAWN indices from Pianosi and Wagener (2015)
for jj= 1:length(NN)
    for k = 1:length(nn)
        for i=1:M         
            
            %subaxis(1,3,i, 'SpacingVert', 0.03,'SpacingHoriz',0.01,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.15,'MarginLeft',0.1); hold on
            subplot(1,3,i); hold on
            
            % add legend to the middle plot:
            if (i==2)&&(jj==1)&&(k==1); 
                for tt=1:length(NN); plot(100,100,'ok','MarkerSize',ms,'MarkerFaceColor',col{tt}); hold on; end; 
                legend(legend_NN );
            end
            % plot PAWN indices:
            plot(nn+pert(jj),PAWN_median{jj}(:,i),'ok','MarkerSize',ms,'MarkerFaceColor',col{jj})
            
            % plot PAWN index for dummy parameter:
            plot(nn+pert(jj),KS_dummy_mean(jj,:),':','Color',col{jj},'LineWidth',2)
            
            % plot confidence interval as a vertical line connecting the
            % upper and lower bound of the PAWN indices:
            for tt = 1:length(nn) 
                plot([nn(tt) nn(tt)]+pert(jj),[PAWN_median_lb{jj}(tt,i) PAWN_median_ub{jj}(tt,i)],'Color',col{jj},'LineWidth',2)
            end
            
            % plot old resultsn (Pianosi and Wagener, 2015) for comparison:
            plot(n_old_res,old_res(i),'ok','MarkerSize',ms,'MarkerFaceColor','k')
            
            % set axis limits and labels        
            axis([min(nn)-1,max(max(nn),n_old_res)+1,0,1])
            set(gca,'XTick',nn,'FontSize',16)
            xlabel('n','FontSize',fs)
            if i == 1; ylabel('median KS','FontSize',fs); else; set(gca,'YTickLabel',{}); end            
            if jj==1; if k==1; text(nn(2),0.8,['input x_' num2str(i) ''],'FontSize',fs); end; end
            
        end
    end
end
set(hfig, 'Position', [0 0 1200 400])


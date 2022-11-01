function [X, AB] = OAT_sampling(r,M,distr_fun,distr_par,samp_strat,des_type)
%
% Build a matrix X of input samples to be used for the Elementary Effects
% Test, using a One-At-the-Time sampling strategy as described in
% Campolongo et al. (2011).
%
% Usage:
%
% X = OAT_sampling(r,M,distr_fun,distr_par,samp_strat,des_type)             
%
%
% Input:
%          r = number of elementary effects  - positive integer number 
%          M = number of inputs              - positive integer number 
%  distr_fun = probability distribution function of each input 
%              - string (eg: 'unif') if all inputs have the same pdf
%              - cell array of M strings (eg:{'unif','norm'}) otherwise
%              See help of AAT_sampling to check supported PDF types
%  distr_par = parameters of the probability distribution function
%              - row vector if all input pdfs have the same parameters
%              - cell array of M vectors otherwise
% samp_strat = sampling strategy             - string                  
%              Options: 'rsu': random uniform                     
%                       'lhs': latin hypercube                    
%   des_type = design type                   - string
%              Options: 'trajectory','radial'
%
% Output:
%          X = matrix of sampling datapoints where EE must be computed.   
%              This is a matrix with r*(M+1) rows and M columns.         
%              Each row is a point in the input space. Rows are sorted in
%              'r' blocks, each including 'M+1' rows. Within each block, 
%              points (rows) differ in one component at the time.        
%              Thus, each block can be used to compute one Elementary    
%              Effect (EE_i) for each model input (i=1,...,M).           
%
% Examples:
%
% % Example 1: 2 inputs, both from Unif[0,3]
% r = 10 ;
% M = 2 ;
% distr_fun = 'unif' ;
% distr_par = [ 0 3 ];
% samp_strat = 'lhs';
% des_type = 'trajectory';
% X = OAT_sampling(r,M,distr_fun,distr_par,samp_strat,des_type);
% % Plot results:
% clrs = jet(r) ;
% figure; hold on; j = 0 ;
% for k=1:r
%     idx =  j+1:j+(M+1) ; j=j+M+1;
%     plot(X(idx,1),X(idx,2),'o:k','MarkerFaceColor',clrs(k,:))
% end
% xlabel('x_1'); ylabel('x_2')
% %
% % Example 2: 2 inputs, one from Unif[0,3], one from Unif[1,5]
% distr_fun = 'unif' ;
% distr_par = {[ 0 3 ],[1 5]};
% X = OAT_sampling(r,M,distr_fun,distr_par,samp_strat,des_type);
% % (use above code to plot results)
% %
% % Example 3: 2 inputs, one from Unif[0,3], one from discrete, uniform in [1,5]
% distr_fun = {'unif','unid'} ;
% distr_par = {[ 0 3 ],5};
% X = OAT_sampling(r,M,distr_fun,distr_par,samp_strat,des_type);
%
% References:
%
% Campolongo F., Saltelli, A. and J. Cariboni (2011), From screening to
% quantitative sensitivity analysis. A unified approach, Computer Physics 
% Communications, 182(4), 978-988.

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isscalar(r); error('''r'' must be a scalar'); end
if r<=0; error('''r'' must be positive' ); end
if abs(r-round(r)); error('''r'' must be integer'); end
if ~isscalar(M); error('''M'' must be a scalar'); end
if M<=0; error('''M'' must be positive' ); end
if abs(M-round(M)); error('''M'' must be integer'); end
% 'distr_fun'
% 'distr_par'
% 'samp_strat'
% will be checked later by the AAT_sampling function
if ~ischar(des_type); error('''des_type'' must be a string'); end

%%%%%%%%%%%%%%%%%%
% Perform sampling
%%%%%%%%%%%%%%%%%%

AB = AAT_sampling(samp_strat,M,distr_fun,distr_par,2*r) ;

X  = nan(r*(M+1),M); % sampling points
k = 1 ;
for i=1:r
    % Sample datapoints:
    ab  = AB((i-1)*2+1:(i-1)*2+2,:);
    a = ab(1,:) ;
    b = ab(2,:) ;
    for j=1:M
        if a(j)==b(j)
            if strcmp(distr_fun{j},'unid') % resample this component
                    while(a(j)==b(j))
                        tmp = AAT_sampling(samp_strat,M,distr_fun,distr_par,2);
                        b(j) = tmp(1,j) ;
                    end
                    fprintf('\n WARNING: b(i=%d,j=%d) was randomly changed not to overlap with a(i=%d,j=%d) \n',[i,j,i,j]);
            else % just print a warning message
                    fprintf('\n WARNING: b(i=%d,j=%d) and a(i=%d,j=%d) are the same! \n',[i,j,a(j),b(j)]);
            end
        end
    end
    
    X(k,:) = a ;
    k=k+1 ;
    if strcmp(des_type,'radial')
        for j=1:M
            x_   = a       ;
            x_(j)= b(j)    ;
            X(k,:) = x_ ;
            if (abs( X(k,j) - X(k-1,j) )==0); fprintf('\n WARNING \n: X(%d,%d) and X(%d,%d) are equal\n',[k,j,k-1,j]); end
            k=k+1 ;
        end
    elseif strcmp(des_type,'trajectory')
        x = a ;
        for j=1:M
            x_   = x     ;
            x_(j)= b(j)  ;
            X(k,:) = x_  ;
            if (abs( X(k,j) - X(k-1,j) )==0); fprintf('\n WARNING \n: X(%d,%d) and X(%d,%d) are equal\n',[k,j,k-1,j]); end
            x = x_  ; % move one step forward in the trajectory
            k = k+1 ;
        end
    else
        error('''des_type'' must be one among {''radial'',''trajectory''}')
    end
end


function [Xext, Xnew] = OAT_sampling_extend(X,r_ext,distr_fun,distr_par,des_type,varargin)
%
% This function create an expanded sample 'Xext' starting from a sample 
% 'X' to be used for the Elementary Effects Test. One-At-the-Time sampling 
% strategy as described in Campolongo et al.(2011) is used. The matrix of
% baseline and auxiliary points is built using latin hypercube and the 
% maximin criterion.
%
% Usage:
%
% Xext = OAT_sampling_extend(X,r_ext,M,distr_fun,distr_par,des_type)             
% Xext = OAT_sampling_extend(X,r_ext,M,distr_fun,distr_par,des_type,n_rep)  
%
% Input:
%          X = initial sample build as described        - matrix(r*(M+1),M)
%              in OAT_sampling                                                                                  
%      r_ext = new number of elementary effects   - positive integer number
%                                                           
%  distr_fun = probability distribution function of each input 
%              - string (eg: 'unif') if all inputs have the same pdf
%              - cell array of M strings (eg:{'unif','norm'}) otherwise
%              See help of AAT_sampling to check supported PDF types
%  distr_par = parameters of the probability distribution function
%              - row vector if all input pdfs have the same parameters
%              - cell array of M vectors otherwise               
%   des_type = design type                         - string
%              Options: 'trajectory','radial'
%       nrep = number of replicate to select the maximin hypercube - scalar
%             (default value: 10)
%
% Output:
%      Xext = expanded sample                      - matrix(r_ext*(M+1),M)
%      Xnew = newly added samples                  - matrix(r_new*(M+1),M)            
%             (where    Xext = [ X; Xnew ]
%                       r_ext = r + r_new )

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

if ~isnumeric(X) ; error('input ''X'' must be a matrix of size (N,M)'); end
[N,M]=size(X) ;
if mod(N,M+1)~=0;error('''X'' must have r*(M+1) lines and M columns');end
r=N/(M+1);
if ~isscalar(r_ext); error('''r_ext'' must be a scalar'); end
if (r_ext - round(r_ext))~=0 ; error('''r_ext'' must be an integer'); end
if r_ext<=0 ; error('''r_ext'' must be a positive integer'); end
if r_ext < r ; error('''r_ext'' must be larger than ''r''');end
% 'distr_fun'
% 'distr_par'
% will be checked later by the AAT_sampling function
if ~ischar(des_type); error('''des_type'' must be a string'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
nrep=10;

% Recover and update optional arguments:
if nargin > 5
    if ~isempty(varargin{1})
        nrep = varargin{1};
        if ~isscalar(nrep); error('''nrep'' must be a scalar'); end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the new matrix of baseline and auxiliary points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AB_old=nan(2*r,M); %old matrix of baseline and auxiliary points in sample X
if strcmp(des_type,'radial')
    for i=1:r
        AB_old(2*(i-1)+1,:)=X((i-1)*(M+1)+1,:);
        for j=1:M
            AB_old(2*i,j)=X((i-1)*(M+1)+j+1,j);
        end
    end
elseif strcmp(des_type,'trajectory')
    for i=1:r
        AB_old(2*(i-1)+1,:)=X((i-1)*(M+1)+1,:);
        AB_old(2*i,:)=X(i*(M+1),:);
    end
end

ABext=AAT_sampling_extend(AB_old,distr_fun,distr_par,2*r_ext,nrep);
%new matrix of baseline and auxiliary points for extended sample Xext

ABnew=ABext(2*r+1:end,:);
r_new=r_ext-r;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the intermediate points to the extension of the sample 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xnew  = nan(r_new*(M+1),M); % sampling points
k = 1 ;
for i=1:r_new
    % Sample datapoints:
    ab  = ABnew((i-1)*2+1:(i-1)*2+2,:);
    a = ab(1,:) ;
    b = ab(2,:) ;
    for j=1:M
        if a(j)==b(j)
            if strcmp(distr_fun{j},'unid') % resample this component
                    while(a(j)==b(j))
                        tmp = AAT_sampling('lhs',M,distr_fun,distr_par,2);
                        b(j) = tmp(1,j) ;
                    end
                    fprintf('\n WARNING: b(i=%d,j=%d) was randomly changed not to overlap with a(i=%d,j=%d) \n',[i,j,i,j]);
            else % just print a warning message
                    fprintf('\n WARNING: b(i=%d,j=%d) and a(i=%d,j=%d) are the same! \n',[i,j,a(j),b(j)]);
            end
        end
    end
    
    Xnew(k,:) = a ;
    k=k+1 ;
    if strcmp(des_type,'radial')
        for j=1:M
            x_   = a       ;
            x_(j)= b(j)    ;
            Xnew(k,:) = x_ ;
            if (abs( Xnew(k,j) - Xnew(k-1,j) )==0); fprintf('\n WARNING \n: X(%d,%d) and X(%d,%d) are equal\n',[k,j,k-1,j]); end
            k=k+1 ;
        end
    elseif strcmp(des_type,'trajectory')
        x = a ;
        for j=1:M
            x_   = x     ;
            x_(j)= b(j)  ;
            Xnew(k,:) = x_  ;
            if (abs( Xnew(k,j) - Xnew(k-1,j) )==0); fprintf('\n WARNING \n: X(%d,%d) and X(%d,%d) are equal\n',[k,j,k-1,j]); end
            x = x_  ; % move one step forward in the trajectory
            k = k+1 ;
        end
    else
        error('''des_type'' must be one among {''radial'',''trajectory''}')
    end
end

Xext=[X;Xnew];
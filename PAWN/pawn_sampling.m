function [ XX, xc ] = pawn_sampling(samp_strat,M,distr_fun,distr_par,n,NC)
%
% This function performs random sampling to implement the PAWN method
% (Pianosi and Wagener, 2015).
%
% Usage:
% [ XX, xc ] = pawn_sampling(samp_strat,M,distr_fun,distr_par,n,NC)
%
% Input:
% samp_strat = sampling strategy                                   - string                  
%              Options: 'rsu': random uniform                     
%                       'lhs': latin hypercube                    
%          M = number of inputs                                    - scalar
%  distr_fun = probability distribution function of each input 
%              - string (eg: 'unif') if all inputs have the same pdf
%              - cell array of M strings (eg:{'unif','norm'}) otherwise
%  distr_par = parameters of the probability distribution function
%              - row vector if all input pdfs have the same parameters
%              - cell array of M vectors otherwise
%          n = number of conditioning points                       - scalar
%         NC = number of samples to be generated                   - scalar
%              at each conditioning point
% 
% Output:
% XX = cell array of size (M,n). The element X{i,k} is the subsamples
%      of input datapoints to compute the sensitivity of factor i in 
%      subsample k, and it is a matrix of size (NC,M)                
% xc = cell array of size (M,1). The element xc{i} is the list of 
%      subsample centers used to compute the sensitivity of factor i
%      and it is a matrix of size (n,1).
%
% REFERENCES
%
% Pianosi, F. and Wagener, T. (2015), A simple and efficient method 
% for global sensitivity analysis based on cumulative distribution 
% functions, Env. Mod. & Soft., 67, 1-11.

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

if ~iscell(distr_fun); tmp=cell(1,M);for i=1:M; tmp{i}=distr_fun; end; distr_fun = tmp; end
if ~iscell(distr_par); tmp=cell(1,M);for i=1:M; tmp{i}=distr_par; end; distr_par = tmp; end

if ~isscalar(M); error('''M'' must be a scalar'); end
if M<=0; error('''M'' must be positive' ); end
if abs(M-round(M)); error('''M'' must be integer'); end

if ~isscalar(n); error('''n'' must be a scalar'); end
if n<=0; error('''n'' must be positive' ); end
if abs(n-round(n)); error('''n'' must be integer'); end

if ~isscalar(NC); error('''NC'' must be a scalar'); end
if NC<=0; error('''NC'' must be positive' ); end
if abs(NC-round(NC)); error('''NC'' must be integer'); end

XX = cell(M,n);
xc = cell(M,1);

for i=1:M

    ifixed    = false(1,M);
    ifixed(i) = true ;
    [ XXi, xci ] = conditional_sampling(samp_strat,M,distr_fun,distr_par,n,NC,ifixed) ;
    xc{i} = xci ;
    for k=1:n; XX{i,k}=XXi{k}; end

end
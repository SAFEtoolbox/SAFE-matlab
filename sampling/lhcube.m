function [X,d] = lhcube(N,M,varargin)
% 
% Generate a latin hypercube of N datapoints in the M-dimensional hypercube
% [0,1]x[0,1]x...x[0,1]. If required, generation can be repeated for a
% prescribed number of times and the maximin latin hypercube is returned.
% 
% Usage:
% 
% X = lhcube(N,M)
% X = lhcube(N,M,nrep)
% [X,d] = lhcube(N,M,nrep)
% 
% Input:
%    N = number of samples                               - positive scalar
%    M = number of inputs                                - positive scalar
% nrep = number of repetition (default: 5)               - positive scalar
%        
% Output:
%    X = sample points                                   - matrix (N,M)
%    d = minimum distance between two points (rows) in X - scalar
% 
% Example:
% 
% N = 10 ;
% M =  2 ;
% X = lhcube(N,M);
% figure
% plot(X(:,1),X(:,2),'.')
% set(gca,'XTick',[0:1/N:1],'XTickLabel',{},'YTick',[0:1/N:1],'YTickLabel',{})
% grid on
% 
% References:
% 
% this is the matlab/octave version of the code in:
%         J.C. Rougier, calibrate package 
%         http://www.maths.bris.ac.uk/~mazjcr/#software
% 
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

if ~isscalar(N); error('''N'' must be a scalar'); end
if N<=1; error('''N'' must be larger than 1' ); end
if abs(N-round(N)); error('''N'' must be integer'); end

if ~isscalar(M); error('''M'' must be a scalar'); end
if M<=0; error('''M'' must be positive' ); end
if abs(M-round(M)); error('''M'' must be integer'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrep=5 ;
if nargin>2
    if ~isempty(varargin{1})
        nrep = varargin{1} ;
        if ~isscalar(nrep); error('''nrep'' must be a scalar'); end
        if nrep<=0; error('''nrep'' must be positive' ); end
        if abs(nrep-round(nrep)); error('''nrep'' must be integer'); end
    end
end

d = 0 ;

% Generate nrep hypercube sample X and keep the one that maximises the 
% minimum inter-point Euclidean distance between any two sampled points:
for k=1:nrep
    
    % Generate a latin-hypercube:
    ran = rand(N,M) ;
    Xk  = zeros(N,M);
    for i=1:M
        idx=randperm(N);
        Xk(:,i)=(idx'-ran(:,i))/N;
    end
    
    % Compute the minimum distance between points in X:
    dk = min(pdist(Xk)) ; % Requires Statistical Toolbox
    % Alternative option:
    %dk=(ipdm(Xk,'metric',2)); dk=min(dk(dk>0));
    
    % If the current latin hypercube has minimum distance higher than
    % the best so far, it will be retained as the best.
    if (dk > d)
        X = Xk ;
        d = dk ;
    end
    
end







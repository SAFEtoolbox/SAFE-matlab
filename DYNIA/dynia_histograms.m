function [ xi, fi ] = dynia_histograms(X,Y,perc,varargin)
%
% Estimate probability distribution of the subset of samples in X
% corresponding to the top 'perc' values in Y (probability
% distribution is approximated by histogram).
% The deviation of this 'posterior' distribution
% from the 'prior' distribution of the entire sample X,
% provides an indication of the strength of the relationship
% between X and Y (or, the sensitivity of Y to X) (Wagener et al, 2003).
%
% Usage:
% [ xi, fi ] = dynia_histograms(X,Y,perc)
% [ xi, fi ] = dynia_histograms(X,Y,perc,nbins)
%
% Input:
%    X = set of input samples                                - vector (N,1)
%    Y = set of output samples                               - matrix (N,T)
% perc = percentage of samples that will                     - scalar
%        be retained to estimate the a posteriori distribution (*)
% nbins = number of bins used to approximate the a posteriori      - scalar
%         distribution by histogram (default: 10)
%
% Output:
% xi = bins centers                                       - vector (nbin,1)
% fi = frequency associated to each bin at each time step - matrix (nbin,T)
%
% (*) The function retains the 'perc' samples with lower-valued output.
%
% References:
% 
% Wagener, T., McIntyre, N., Lees, M., Wheater, H., Gupta, H., 2003. 
% Towards reduced uncertainty in conceptual rainfall-runoff modelling: 
% dynamic identifiability analysis. Hydrol. Process. 17, 455?476.

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

if ~isnumeric(X) ; error('input ''X'' must be a vector of size (N,1)'); end
if ~isnumeric(Y) ; error('input ''Y'' must be a matrix of size (N,T)'); end
[N,M]=size(X) ;
if M>1 ; error('input ''X'' must be a vector of size (N,1)'); end
[n,T]=size(Y) ;
if N~=n; error('input ''X'' and ''Y'' must have the same number of rows'); end
if T<=1; error('input ''Y'' must have at least 2 columns'); end

if ~isnumeric(perc) ; error('input ''perc'' must be a scalar number'); end
if ~isscalar(perc); error('''perc'' must be scalar'); end
if perc <=0 ; error('''perc'' must be higher than 0'); end
if perc >=100 ; error('''perc'' must be lower than 100'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
nbins  = 10 ;

% Recover and update optional arguments:
if nargin > 3
    if ~isempty(varargin{1})
        nbins = varargin{1} ;
        if ~isnumeric(nbins) ; error('input ''perc'' must be a scalar number'); end
        if ~isscalar(nbins); error('''nbins'' must be scalar'); end
        if nbins<1; error('''nbins'' must be larger than 1' ); end
        if abs(nbins-round(nbins)); error('''nbins'' must be an integer'); end
    end
end

%%%%%%%
% DYNIA
%%%%%%%

% Find top-perc samples:

[Y_sorted,idx] = sort(Y) ;
n = floor(N*perc/100) ;

% Sort samples according to associated performances:
X_sorted = X(idx) ; % (N,T)
% Take the first 'n' parameter sets at each time t=1,...,T
X_best = X_sorted(1:n,:)  ; % (n,T)

% Because we want function 'hist' to define the bins in the same way
% at all time steps, i.e. covering the a priori variability range of
% the input, we add two fake rows filled in with the minimum and 
% maximum possible values of that input:
xmin = min(X) ;
xmax = max(X) ;
X_best_mod = [ ones(1,T)*xmin ; X_best ; ones(1,T)*xmax ] ;

% Estimate and frequency distribution of the input samples that
% generated the 'perc' top performances:
[ni,xi] = hist(X_best_mod,nbins) ;

% Then, remove the two 'fake' counts from the first and last bin:
ni(1,:)   = ni(1,:)-1 ;
ni(end,:) = ni(end,:)-1 ;

% Finally, we compute the frequency:
fi = ni / n ;


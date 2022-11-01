function [ idx, Llim, Ulim ] = GLUE(GLF,threshold,Y_sim,varargin)
%
% This function implements the Generalized Likelihood Uncertainty 
% Estimation (GLUE) method, as first proposed by Beven and Binley (1992) 
% and Beven and Freer (2001)).
% It retains as 'behavioural' the simulations associated to values
% of the Generalized Likelihood Function(s) (GLFs) below a given
% threshold ('below' because we assume that GLFs measure the distance
% between observations and predictions) i.e.:
%
%                  GLF(i,j)<threshold(j)     for j=1,...,P
% 
% Then, the function computes the prediction limits for each time step 
% using the GLF values of 'behavioural' simulations as a measure of
% 'probablity'
%
% Usage:
% [ idx, Llim, Ulim ] = GLUE(GLF,threshold,Y_sim)
% [ idx, Llim, Ulim ] = GLUE(GLF,threshold,Y_sim,alfa)
%
%
% Input:
%     Y_sim = time series of modelled variable               - matrix (N,T)
%             (N is the number of sampled parameterizations,
%              T is the length of each simulated time series)
%         GLF = values of generalized likelihood function(s) - matrix (N,P)
%             (P>1 for multiple GLFs are considered simoultaneously)
% threshold = threshold for the GFLs                         - vector (1,P)
%             (if not specified: threshold=median(Y) )   
%      alfa = significance level for the predictions               - scalar 
%             limits (default: 0.05)
%
% Output:
%
%    idxb = indices of samples statisfying the condition     - vector (N,1)
%           (logical vector: 1 if behavioural, 0 otherwise)
%    Llim = time series of lower prediction limit            - vector (T,1)
%    Ulim = time series of upper prediction limit            - vector (T,1)
%
% REFERENCES
% 
% Beven, K. and Binley, A. (1992), The future of distributed models: 
% Model calibration and uncertainty prediction. Hydrol. Process.,6,279-298
%
% Beven, K.GLF and Freer, J.E. (2001), Equifinality, data assimilation, and 
% uncertainty estimation in mechanistic modelling of complex environmental 
% systems using the GLUE methodology. Journal of Hydrology, 249(1-4),11-29

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

if ~isnumeric(GLF) ; error('GLF must be numeric'); end
if ~isnumeric(threshold) ; error('''threshold'' must be numeric'); end
if ~isnumeric(Y_sim) ; error('Y_sim must be numeric'); end
[N,T] = size(Y_sim) ;
[N2,P] = size(GLF) ;
[f,P2] = size(threshold) ;
if N2~=N
    error('GLF and Y_sim must have the same number of rows')
end
if P2~=P
    error('GLF and threshold must have the same number of columns')
end
if f>1
    error('threshold must be a row vector')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
alfa=0.05;

% Recover and update optional arguments:

if nargin > 3
    if ~isempty(varargin{1})
        alfa=varargin{1};
        if ~isscalar(alfa); error('''alfa'' must be scalar'); end
        if any([alfa<0,alfa>1]); error('''alfa'' must be in [0,1]' ); end
    end
end

%%%%%%%%%%%%%
% GLUE
%%%%%%%%%%%%

idx = sum(GLF<+repmat(threshold,N,1),2)==P ;% indices of behavioural samples
% idx = GLF >= threshold ; % old (scalar case)

if isempty(idx)
    
    warning('No samples satisfying the condition GLF>=%g',threshold);    
    Llim = []; Ulim= [];
    
else
    
    % Define "likelihood" measure:
    Jb = GLF ; % take the objective function 'GLF' as likelihood function
    Jb(~idx) = 0 ; % set to 0 the likelihood of non-behavioural samples
    Jb = Jb/(sum(Jb)); % rescale
    
    % Find "CDFs" (and prediction limits) at all time steps:
    Jbt  = Jb(idx) ; % (Nb,1)
    Llim = nan(T,1);
    Ulim = nan(T,1);
    for t=1:T
        [y_sorted,idx_sort] = sort(Y_sim(idx,t));
        CDF_t = cumsum(Jbt(idx_sort));
        % Find lower limit:
        Llim_ = y_sorted(CDF_t<alfa);
        if ~isempty(Llim_);
            Llim(t) = Llim_(end);
        else
            Llim(t) = y_sorted(1);
        end
        % Find upper limit:
        Ulim_ = y_sorted(CDF_t>1-alfa);
        if ~isempty(Ulim_);
            Ulim(t) = Ulim_(1);
        else
            Ulim(t) = y_sorted(end);
        end
    end
    
end


function X = Morris_sampling(r,xmin,xmax,L)
%
% Build a matrix X of input samples to be used for the Elementary Effects
% Test, using the One-At-the-Time sampling strategy originally proposed by
% Morris (1991). It implicitely assumes that all the inputs be uncorrelated
% and drawn from a continuous, uniform distribution function. 
% 
% Usage:
% X = Morris_sampling(r,xmin,xmax,L)                                       
% 
% Input:
%    r = number of elementary effects          - positive integer number 
% xmin = lower bounds of input ranges          - vector (1,M)            
% xmax = upper bounds of input ranges          - vector (1,M)            
%    L = number of levels in the sampling grid - positive, even number   
% 
% Output:
%    X = matrix of sampling datapoints where EE must be computed.   
%        This is a matrix with r*(M+1) rows and M columns.         
%        Each row is a point in the input space. Rows are sorted in
%        'r' blocks, each including 'M+1' rows. Within each block, 
%        points (rows) differ in one component at the time.        
%        Thus, each block can be used to compute one Elementary    
%        Effect (EE_i) for each model input (i=1,...,M).           
% 
% References:
% 
% Morris, M.D. (1991), Factorial sampling plans for preliminary          
% computational experiments, Technometrics, 33(2).                       
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

if ~isscalar(r); error('''r'' must be a scalar'); end
if r<=0; error('''r'' must be positive' ); end
if abs(r-round(r)); error('''r'' must be integer'); end
[N,M] = size(xmin)   ;
[n,m] = size(xmax) ;
if N~=1 ;error('''xmin'' must be a row vector'); end
if n~=1 ;error('''xmax'' must be a row vector'); end
if M~=m ;error('''xmin'' and ''xmax'' must be the same size'); end
Dr = xmax - xmin ;
if any(Dr<=0)
    error('all components of ''xmax'' must be higher than the corresponding ones in ''xmin''')
end

if ~isscalar(L); error('''L'' must be scalar'); end
if L<=0; error('''L'' must be positive' ); end
if abs(L-round(L)); error('''L'' must be integer'); end
if abs(L-ceil(L/2)*2);
    L = ceil(L/2)*2  ; %
    fprintf('WARNING: ''L'' must be even!\n Using L=%g instead of user-defined value\n',L);
end

%%%%%%%%%%%%%%%%%%
% Perform sampling 
%%%%%%%%%%%%%%%%%%

X  = nan(r*(M+1),M); % sampling points
k = 1 ;
for i=1:r
    % Sample datapoints:
    Bstar = Morris_orientation_matrix(M,L); % (M+1,M)
    % Resort to original ranges:
    Bstar = repmat(xmin,M+1,1) + Bstar .* repmat(Dr,M+1,1);
    for j=1:M+1
        X(k,:) = Bstar(j,:) ; k=k+1 ;
    end
end

function Bstar = Morris_orientation_matrix(k,p)
%
% Bstar = Morris_orientation_matrix(k,p)
%
%     k = number of inputs    - integer positive scalar
%     p = number of levels    - integer positive even scalar
% Bstar = matrix of (k+1) datapoints in the k-dimensional 
%         input space (to be used for computing one Elementary
%         Effect for each of the 'k' inputs)
%
% Example in two-dimensional space (k=2):
%
% p = 4 ;
% figure
% Bstar = Morris_orientation_matrix(2,p); plot(Bstar(:,1),Bstar(:,2),'.r-')
% % if you want to generate more datapoints:
% hold on
% Bstar = Morris_orientation_matrix(2,p); plot(Bstar(:,1),Bstar(:,2),'xb-')
% Bstar = Morris_orientation_matrix(2,p); plot(Bstar(:,1),Bstar(:,2),'oc-')
% Bstar = Morris_orientation_matrix(2,p); plot(Bstar(:,1),Bstar(:,2),'sm-')
% set(gca,'XTick',[0:1/(p-1):1],'YTick',[0:1/(p-1):1])
% axis([-0.1,1.1,-0.1,1.1]); box on; grid on

if ~isscalar(k); error('''k'' must be a scalar'); end
if k<=0; error('''k'' must be positive' ); end
if abs(k-round(k)); error('''k'' must be integer'); end
if ~isscalar(p); error('''p'' must be a scalar'); end
if p<=0; error('''p'' must be positive' ); end
if abs(p-round(p)); error('''p'' must be integer'); end
if abs(p-ceil(p/2)*2);
    p = ceil(p/2)*2  ; %
    fprintf('WARNING: input ''p'' must be even!\n Using p=%g instead of user-defined value\n',p);
end

m=k+1;
Delta=p/(2*(p-1));
% sampling matrix:
B = tril(ones(m,k),-1); % (m,k)

% Create diagonal matrix with {-1,1} elements
tmp=randi([0 1],k,1) ; % random numbers from {0,1}
tmp(tmp==0)=-1; % random numbers from {-1,1}
D=diag(tmp);

% Create base value vector
In = nan(1,p/2);
for i=0:(p/2-1)
    In(i+1)=i/(p-1);
end
tmp=randi([1 length(In)],1,k);
x=In(tmp);

% Create random permutation matrix
P = eye(k);
idx = randperm(k);
P = P(idx,:);

Jmk=ones(m,k);
Jm1=ones(m,1);

% Create a random orientation of B:
Bstar = ( Jm1*x + Delta/2*((2*B-Jmk)*D+Jmk) )*P;

function [KS_median,KS_mean,KS_max,KS_dummy,YY,xc,NC,XX] = pawn_indices_givendata(X,Y,n,nboot)
%
% Compute PAWN sensitivity indices as a statistic of the KS
% between the unconditional and conditional output distributions,
% according to the approximation strategy by Pianosi and Wagener (2018).
%
% Usage:
% [KS_median,KS_mean,KS_max,KS_dummy,YY,xc,NC,XX] = pawn_indices_givendata(X,Y,n,nboot)
%
% Input:
%         X = set of inputs samples                        - matrix (N,M)
%         Y = set of output samples                        - matrix (N,1)
%         n = number of conditional intervals              - scalar
%     nboot = number of bootstrap resamples to derive      - scalar
%             confidence intervals
%
% Output:
% KS_median = median KS across n conditioning points       - size (nboot,M)    
%             (one value for each input and each bootstrap resample)                                                             
% KS_mean   = mean KS across n conditioning points         - size (nboot,M)
% KS_max    = max KS across n conditioning points          - size (nboot,M)
% KS_dummy  = KS of dummy parameter                        - size (nboot,1)
%             (one value for each bootstrap resample)
%  YY = output samples for the conditional distributions   - cell array (M,n) 
%       The element Y{i,k} is the k-th subsample 
%       of output evaluations obtained by fixing the i-th input 
%       (or group of inputs) to its k-th conditioning interval,
%       and it is a matrix of size (NC(i,k),1).
% xc = mean value of conditioning intervals                - cell array (M,1)
%      The element xc{i} is the list of 
%      subsample centers (mean points of the conditioning intervals) 
%      used to ccondition input i, and it is a matrix of size (n,1).
% NC = number of conditional output samples                - matrix (M,n)
%      for each input and each conditional interval 
% XX = input samples to derive the conditional samples     - cell array (M,n)
%      The element X{i,k} is the k-th subsample
%      of input values where the i-th input is fixed to
%      the k-th conditioning interval (while other inputs vary freely), 
%      and it is a matrix of size (NC(i,k),M)
%
% REFERENCES
%
% Pianosi, F. and Wagener, T. (2018), Distribution-based sensitivity 
% analysis from a generic input-output sample, Env. Mod. & Soft.

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

% Check inputs
if ~isnumeric(X) ; error('input ''X'' must be numeric'); end
if ~isnumeric(Y) ; error('input ''Y'' must be numeric'); end
[N,M]=size(X) ;
[N2,m]=size(Y) ;
if N~=N2; error('input ''X'' and ''Y'' must have the same number of rows'); end
if m~=1; error('input ''Y'' must be a column vector'); end
if ~isscalar(n) ; error('input ''n'' must be scalar'); end
if abs(n-round(n)); error('''n'' must be integer'); end
if n<1  ; error('input ''n'' must be a positive integer'); end
if ~isscalar(nboot) ; error('input ''nboot'' must be scalar'); end
if abs(nboot-round(nboot)); error('''nboot'' must be integer'); end
if nboot<1  ; error('input ''nboot'' must be a positive integer'); end

%
[YY,xc,NC,XX] = pawn_split_sample(X,Y,n);
bootsize = round(mean(mean(NC)));
KS_median = nan(nboot,M) ;
KS_mean = nan(nboot,M)   ;
KS_max = nan(nboot,M)    ;
KS_dummy = nan(nboot,1)  ;
  
for k=1:nboot
    
    % Bootrstrap from unconditional sample:
    idx_bootstrap = randperm(N,bootsize) ;
    YU = Y(idx_bootstrap)        ; % (bootsize,1)    
    % Compute empirical CDFs of unconditional and conditional samples:
    [ YF, Fu, Fc  ] = pawn_cdfs(YU,YY) ;
    % Compute KS among CDFs:
    KS  = pawn_ks(YF,Fu,Fc)      ; % (n,M)
    
    % Take a statistic of KS across x_i values:
    KS_median(k,:) = median(KS)  ; % (1,M)
    KS_mean(k,:)   = mean(KS)    ; % (1,M)
    KS_max(k,:)    = max(KS)     ; % (1,M)    
    
    % Compute KS statistic for dummy parameter:   
    % Bootrstrap again from unconditional sample:
    YU2 = Y(randperm(N,bootsize))  ; % (bootsize,1)
    % Compute empirical CDFs of the two unconditional samples:
    [ YF, Fu, Fc  ] = pawn_cdfs(YU,{YU2}) ;
    % Compute KS among CDFs:
    KS_dummy(k) = pawn_ks(YF,Fu,Fc)       ; 

    %    figure(100); hold on; plot(YF,Fu,'.-b',YF,Fc{1},'.-m')

end

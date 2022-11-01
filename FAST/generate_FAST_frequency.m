function omega = generate_FAST_frequency(M)
%
% Generates a sequence of M frequency values (omega) for sampling 
% according to the FAST method
% (See also FAST_sampling_unif.m for details and references about FAST 
% sampling strategy)
%
% omega = generate_FAST_frequency(M)
%
%     M = number of inputs (integer scalar between 4 and 50)
% omega = frequency set free of interferences through (at least) 4th order
%         (vector (1,M))
%
% For M>4, frequencies are computed based on the recursive algorithm by:
%
% Cukier et al. (1975) Study of the sensitivity of coupled reaction systems
% to uncertainties in rate coefficients. III. Analysis of the 
% approximations, J. Chem. Phys. 63, 1140
%
% which is free of interferences through the 4th order.
% For M<=4, we use values from the literature that guarantee higher order
% interferences free (see comments in the code for specific references)

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

if ~isscalar(M); error('''M'' must be a scalar'); end
if (M - round(M))~=0 ; error('''M'' must be an integer'); end
if M<2 ; error('''M'' must be >4'); end
if M>50; error('''M'' must be <25'); end

if M==2 % Use values from Sec. 3.1 in:
% Xu, C. and G. Gertner (2007), Extending a global sensitivity analysis 
% technique to models with correlated parameters, Computational Statistics
% and Data Analysis, 51, 5579-5590.
    omega = [ 5 23 ] ;
    % (free of interference through 10th order)
elseif M==4 % Use values from Table III in Cukier et al. (1975)
    % (free of interferences through 6th order)
    omega = [13 31 37 41 ];
else % Use recursive algorithm in the same paper
    Omega = [ 0 0 1 5 11 1 17 23 19 25 41 31 23 87 67 73 85 143 149 99 119 ...
        237 267 283 151 385 157 215 449 163 337 253 375 441 673 773 875 873 ...
        587 849 623 637 891 943 1171 1225 1335 1725 1663 2019 ] ;
    d     = [ 4 8 6 10 20 22 32 40 38 26 56 62 46 76 96 60 86 126 134 112 ...
        92 128 154 196 34 416 106 208 328 198 382 88 348 186 140 170 284 ...
        568 302 438 410 248 448 388 596 216 100 488 166 ] ;
    % above values taken from Table VI  
    omega = nan(1,M) ;
    omega(1) = Omega(M) ;
    for i=2:M
        omega(i)=omega(i-1)+d(M+1-i);
        % equation (5.1)
    end
end

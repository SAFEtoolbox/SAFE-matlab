function [Q_sim,STATES,FLUXES] = hymod_sim(param,prec,evap)
%
% This function simulates the Hymod rainfall-runoff model
% (Boyle, 2001; Wagener et al., 2001).
%
% Usage:
%
% [Q_sim,STATES,FLUXES] = hymod_sim(param,prec,evap)
%
% param = 5 elements vector of model parameters:             - vector (1,5)
%         1. Sm   = maximum soil moisture [mm]
%         2. beta = exponent in the soil moisture 
%                    routine [-]
%         3. alfa = partition coefficient [-]
%         4. Rs   = slow reservoir coefficient [1/Dt]
%         5. Rf   = fast reservoir coefficient [1/Dt] 
%  prec = time series of rainfall                            - vector (T,1)
%  evap = time series of potential evaporation               - vector (T,1)
%
%  Q_sim = time series of simulated flow                     - vector (T,1)
% STATES = time series of simulated states                   - matrix (T,5)
% FLUXES = time series of simulated fluxes                   - matrix (T,4)
%          (see comments in the code for the definition of state and flux
%          variables)
% 
% Recommended parameter values:  Smax should fall in [0,400] (mm)
%                                beta should fall in [0,2] (-)
%                                alfa should fall in [0,1] (-)
%                                  Rs should fall in [0,k] (-)
%                                  Rf should fall in [k,1] (-) with 0<k<1
% References:
%
% Boyle, D. (2001). Multicriteria calibration of hydrological models. 
% PhD thesis, Dep. of Hydrol. and Water Resour., Univ. of Ariz., Tucson.
%
% Wagener, T., Boyle, D., Lees, M., Wheater, H., Gupta, H., and Sorooshian, 
% S. (2001). A framework for development and application of hydrological 
% models. Hydrol. Earth Syst. Sci., 5, 13-26.

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/ 

% --------------------------
% Recover model parameters:
% --------------------------

Sm   =max(eps,param(1)) ; % Maximum Soil Moisture (cannot be zero! 
                          % See lines 35 and 37)
beta =param(2) ; % Exponential parameter in soil routine [-]
alfa =param(3) ; % Partitioning factor [-]
Rs   =param(4) ; % Slow reservoir outflow coefficient (ratio) [1/Dt]  
Rf   =param(5) ; % Fast reservoir outflow coefficient (ratio) [1/Dt] 


N = length(prec);% number of time samples

% -----------------------
% Initialize variables:
% ----------------------

Pe = zeros(N,1);  % Recharge from the soil [mm/Dt]
Ea = zeros(N,1);  % Actual Evapotranspiration [mm/Dt]
sm = zeros(N+1,1);  % Soil Moisture [mm]
sL = zeros(N+1,1);  % Slow reservoir moisture [mm]
sF1 = zeros(N+1,1); % Fast reservoir 1 moisture [mm]
sF2 = zeros(N+1,1); % Fast reservoir 2 moisture [mm]
sF3 = zeros(N+1,1); % Fast reservoir 3 moisture [mm]
QsL = zeros(N,1); % Slow flow [mm/Dt]
QsF = zeros(N,1); % Fast flow [mm/Dt]

for t=1:N
    
% --------------------------
%   Soil Moisture Dynamics:
% --------------------------
    
    F  = 1 - ( 1 - sm(t)/Sm )^beta ;
    Pe(t)  = F*prec(t) ; % Compute the value of the outflow 
    %(we assumed that this process is faster than evaporation)
    sm_temp  = max(min(sm(t) + prec(t) - Pe(t),Sm),0);% Compute the water 
    % balance with the value of the outflow 
    Pe(t)=Pe(t)+ max(sm(t) + prec(t) - Pe(t)-Sm,0)+min(sm(t) + prec(t) - Pe(t),0);
    % adjust Pe by an amount equal to the possible negative sm amount or 
    % to the possible sm amount above Sm.
    
    W = min(abs( sm(t)/Sm ),1) ; %Correction factor for evaporation
    Ea(t)= W*evap(t) ; % Compute the evaporation
    sm(t+1) = max(min(sm_temp-Ea(t),Sm),0); % Compute the water balance 
    Ea(t)= Ea(t)+ max(sm_temp-Ea(t)-Sm,0)+min(sm_temp-Ea(t),0);%adjust Ea 
    % by an amount equal to the possible negative sm amount or to the 
    % possible sm amount above Sm 
    
% -------------------------
%    Groundwater Dynamics:
% -------------------------

    % slow flow
    QsL(t)   = Rs * sL(t) ;
    sL(t+1)  = sL(t) +  (1-alfa)*Pe(t) - QsL(t) ;
    % fast flow
    sF1(t+1) = sF1(t) +  alfa*Pe(t) - Rf * sF1(t) ;
    sF2(t+1) = sF2(t) +  Rf*sF1(t) - Rf * sF2(t) ;
    sF3(t+1) = sF3(t) +  Rf*sF2(t) - Rf * sF3(t) ;
    QsF(t)   = Rf * sF3(t) ;

end

Q_sim = QsL+QsF;

STATES=[sm,sL,sF1,sF2,sF3];
FLUXES=[Pe,Ea,QsL,QsF];
    
end


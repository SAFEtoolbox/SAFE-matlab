function [P,STATES,FLUXES] = snow_routine(temp,prec,param)
%
% This function simulates a simple, conceptual snow accumulation/melting
% model based on a degree day approach.
%
% Usage:
%
% [P,STATES,FLUXES] = snow_routine(temp,prec,param)
%
% Input:
%   temp = time series of temperature                        - vector (T,1)
%   prec = time series of precipitation                      - vector (T,1)
%  param = model parameters                                  - vector (1,4)
%            1. Ts    = threshold temperature [C]
%            2. CFMAX = degree day factor [mm/C]
%            3. CFR   = refreezing factor [-]
%            4. CWH   = Water holding capacity of snow [-]       
%
% Output:
%      P = time series of simulated flow exiting from the    - vector (T,1)
%          snowpack (as a result of melt-refreezing) [mm/Dt]        
% STATES = time series of simulated storages (all in mm)     - matrix (T,2)
%          1st column: water content of snowpack (snow component)
%          2nd column: water content of snowpack (liquid component)
% FLUXES = time series of simulated fluxes (all in mm/Dt)    - matrix (T,2)
%          1st column: refreezing
%          2nd column: snowmelt

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

Ts    = param(1); % threshold temperature [C]
CFMAX = param(2); % degree day factor [mm/C]
CFR   = param(3); % refreezing factor [-]
CWH   = param(4); % Water holding capacity of snow [-]


N = length(prec);% number of time samples

% ----------------
% SNOWPACK ROUTINE
% ----------------

P  = zeros(N,1); % snowmelt leaving the snowpack/recharge to the soil [mm/Dt]
  
rain = prec; rain(temp<Ts) = 0  ; % [mm/Dt]
snow = prec; snow(temp>=Ts) = 0 ; % [mm/Dt]
Ta   = temp-Ts; Ta(temp<Ts) = 0 ; % Active Temperature for snowmelt
Tn   = Ts-temp; Tn(temp>=Ts) = 0 ; % Active Temperature for refreezing
m  = zeros(N,1); % snowmelt [mm/Dt]
rfz= zeros(N,1); % refreezing [mm/Dt]
v  = zeros(N+1,1); % snowpack depth [mm]: solid component
vl = zeros(N+1,1); % snowpack depth [mm]: liquid component   

for t=1:N
        
    m(t)  = min(CFMAX*Ta(t),v(t)) ;                        
    rfz(t)= min(CFR*CFMAX*Tn(t) ,vl(t)) ;                   
    %   snowpack dynamics: solid component
    v(t+1)  = v(t)- m(t)+snow(t)+rfz(t); 
    %   snowpack dynamics: liquid component
    vl(t+1) = vl(t)+m(t)+rain(t)-rfz(t); 
    if vl(t+1)>CWH*v(t+1) % if the liquid component exceed the snow pack
        % holding capacity
        P(t)  = vl(t+1)-CWH*v(t+1) ;
        vl(t+1) = CWH*v(t+1) ;
    else
        P(t)=0;
    end
end
STATES=[v,vl];
FLUXES=[rfz,m];
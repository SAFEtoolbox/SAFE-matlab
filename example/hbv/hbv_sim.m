function [Q_sim,STATES,FLUXES] = hbv_sim(P,ept,param,Case)
%
% This function simulates the HBV rainfall-runoff model (Seibert, 1997).
%
% [Q_sim,STATES,FLUXES] = hbv_sim(P,ept,param,Case)
%
%     P = time series of precipitation                      - vector (T,1)
%   ept = time series of potential evapotranspiration       - vector (T,1)
% param = vector of model parameters                        - vector (1,9)
%            1. BETA  = Exponential parameter in soil routine [-]
%            2. LP    = evapotranspiration limit [-]
%            3. FC    = field capacity [mm] 
%            4. PERC  = maximum flux from Upper to Lower Zone [mm/Dt]
%            5. K0    = Near surface flow coefficient (ratio) [1/Dt]  
%            6. K1    = Upper Zone outflow coefficient (ratio) [1/Dt]  
%            7. K2    = Lower Zone outflow coefficient (ratio) [1/Dt]  
%            8. UZL   = Near surface flow threshold [mm]
%            9. MAXBAS= Flow routing coefficient [Dt]
% Case = flag for preferred path in the Upper Zone dynamics - scalar
%        flag=1 -> Preferred path is runoff 
%        flag=2 -> Preferred path is percolation
%
%  Q_sim = time series of simulated flow                    - vector (T,1)
% STATES = time series of simulated storages (all in mm)    - matrix (T,3)
%          1: water content of soil (soil moisture)
%          2. water content of upper reservoir of flow routing routine 
%          3. water content of lower reservoir of flow routing routine
% FLUXES = time series of simulated fluxes (all in mm/Dt)   - matrix (T,5)
%          1: actual evapotranspiration
%          2: recharge (water flux from soil moisture accounting module
%             to flow routing module)
%          3: percolation (water flux from upper to lower reservoir of the 
%             flow routing module)
%          4: runoff from upper reservoir
%          5: runoff from lower reservoir
%
% References: 
%
% Seibert, J.(1997), Estimation of Parameter Uncertainty in the HBV Model,
% Nordic Hydrology, 28(4/5), 247-262.
%
% Comments:
% * The Capillary flux (from upper tank to soil moisture accounting module)
% is not considered
% * The recharge from the soil to the upper zone is considered to be a
% faster process than evapotranspiration.
% * The preferential path from the upper zone can be modified 
%           - Case 1: interflow is dominant
%           - Case 2: percolation is dominant

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

BETA  = param(1); % Exponential parameter in soil routine [-]
LP    = param(2); % evapotranspiration limit [-]
FC    = max(eps,param(3)); % field capacity [mm] cannot be zero

PERC  = param(4); % maximum flux from Upper to Lower Zone [mm/Dt]
K0    = param(5); % Near surface flow coefficient (ratio) [1/Dt]  
K1    = param(6); % Upper Zone outflow coefficient (ratio) [1/Dt]  
K2    = param(7); % Lower Zone outflow coefficient (ratio) [1/Dt]  
UZL   = param(8); % Near surface flow threshold [mm]

MAXBAS= max(1,round(param(9)));% Flow routing coefficient [Dt]


N = length(ept);% number of time samples


% ---------------------
% SOIL MOISTURE ROUTINE
% ---------------------

EA = zeros(N,1) ; % Actual Evapotranspiration [mm/Dt]
SM = zeros(N+1,1) ; % Soil Moisture [mm]
R  = zeros(N,1) ; % Recharge (water flow from Soil to Upper Zone) [mm/Dt]
UZ = zeros(N+1,1) ; % Upper Zone moisture [mm]
LZ = zeros(N+1,1) ; % Lower Zone moisture [mm]
RL = zeros(N,1) ; % Recharge to the lower zone [mm]
Q0 = zeros(N,1) ; % Outflow from Upper Zone [mm/Dt]
Q1 = zeros(N,1) ; % Outflow from Lower Zone [mm/Dt]

for t=1:N
    
% --------------------------
%    Soil Moisture Dynamics:
% --------------------------

    R(t)= P(t)*(SM(t)/FC)^BETA;  % Compute the value of the recharge to the 
    %upper zone (we assumed that this process is faster than evaporation)
    SM_dummy = max(min(SM(t)+P(t)-R(t),FC),0); % Compute the water balance 
    % with the value of the recharge  
    R(t)=R(t)+ max(SM(t)+P(t)- R(t)-FC,0)+min(SM(t)+P(t)-R(t),0);%adjust R 
    % by an amount equal to the possible negative SM amount or to the 
    % possible SM amount above FC
    
    EA(t)=ept(t)*min(SM_dummy/(FC*LP),1);% Compute the evaporation
    SM(t+1) = max(min(SM_dummy-EA(t),FC),0); % Compute the water balance 
    
    EA(t)=EA(t)+ max(SM_dummy-EA(t)-FC,0)+min(SM_dummy-EA(t),0);%adjust EA
    % by an amount equal to the possible negative SM amount or to the 
    % possible SM amount above FC
     
% --------------------
% Upper Zone dynamics:
% --------------------

    if Case==1
    % Case 1: Preferred path = runoff from the upper zone 
  	    Q0(t) = max(min(K1*UZ(t)+K0*max(UZ(t)-UZL,0),UZ(t)),0) ;     
        RL(t) = max(min(UZ(t)-Q0(t),PERC),0);

    elseif Case==2
    % Case 2: Preferred path = percolation
        RL(t)= max(min(PERC,UZ(t)),0);
        Q0(t) = max(min(K1*UZ(t)+K0*max(UZ(t)-UZL,0),UZ(t)-RL(t)),0);
    else
        error('Case must equal to 1 or 2 ')
    end

  UZ(t+1) = UZ(t)+R(t)-Q0(t)-RL(t);

% --------------------
% Lower Zone dynamics: 
% --------------------

    Q1(t) = max(min(K2*LZ(t),LZ(t)),0);
    LZ(t+1) = LZ(t)+RL(t)-Q1(t);
   
end
 
Q = Q0 + Q1 ; % total outflow (mm/Dt)

% --------------------
% FLOW ROUTING ROUTINE
% --------------------

%c = trimf(1:MAXBAS,[0 (MAXBAS+1)/2 MAXBAS+1]) ; %(Seibert,1997)
c = mytrimf(1:MAXBAS,[0 (MAXBAS+1)/2 MAXBAS+1]); %(Seibert,1997)


c = c/sum(c); % vector of normalized coefficients - (1,MAXBAS)
Q_sim = zeros(N,1);
for t=MAXBAS:N
    Q_sim(t) = c * Q(t-MAXBAS+1:t) ; % (Seibert,1997)
end


STATES=[SM,UZ,LZ];
FLUXES=[EA,R,RL,Q0,Q1];

function f = mytrimf(x,param)
% implements triangular-shaped membership function
% (available in Matlab Fuzzy Logic Toolbox as 'trimf')

f = zeros(size(x)) ; 
idx = (x>param(1)) & (x<=param(2)) ;
f(idx) = (x(idx)-param(1))/(param(2)-param(1)) ;
idx = (x>param(2)) & (x<=param(3)) ;
f(idx) = (param(3)-x(idx))/(param(3)-param(2)) ;

    


function [f,Q_sim,STATES,FLUXES]  = hbv_snow_objfun(param,prec,temp,ept,flow,warmup,Case)
%
% This function simulates the snow accumulation/melting process
% (via the function ''snow_routine'') and the rainfall-runoff process
% (via the HBV model by Seibert (1997)) and returns 6 objective functions
% (see Kollat et al, 2012).
%
% [f,Q_sim,STATES,FLUXES]  = hbv_snow_objfun(param,prec,temp,ept,flow,warmup,Case)
% 
% Input:
% param = vector of model parameters                        - vector (1,13)
%         Snow routine parameters:
%            1. Ts    = threshold temperature [C]
%            2. CFMAX = degree day factor [mm/C]
%            3. CFR   = refreezing factor [-]
%            4. CWH   = Water holding capacity of snow [-]       
%         HBV parameters:
%            5. BETA  = Exponential parameter in soil routine [-]
%            6. LP    = evapotranspiration limit [-]
%            7. FC    = field capacity [mm] 
%            8. PERC  = maximum flux from Upper to Lower Zone [mm/Dt]
%            9. K0    = Near surface flow coefficient (ratio) [1/Dt]  
%           10. K1    = Upper Zone outflow coefficient (ratio) [1/Dt]  
%           11. K2    = Lower Zone outflow coefficient (ratio) [1/Dt]  
%           12. UZL   = Near surface flow threshold [mm]
%           13. MAXBAS= Flow routing coefficient [Dt]
%   prec = time series of precipitation                      - vector (T,1)
%   temp = time series of temperature                        - vector (T,1)
%    ept = time series of evapotranspiration                 - vector (T,1)
%   flow = time series of observed flow                      - vector (T,1)
% warmup = number of time steps for model warm-up            - scalar
%   Case = flag [see hbv_sim.m for more details]             - scalar
%
% Output:
%      f = vector of objective functions                     - vector (1,6)
%          1:AME, 2:NSE, 3:BIAS, 4:TRMSE, 5:SFDCE, 6:RMSE
%  Q_sim = time series of simulated flow                     - vector (T,1)
% STATES = time series of simulated storages (all in mm)     - matrix (T,5)
%          1: water content of snowpack (snow component)
%          2: water content of snowpack (liquid component)
%          3: water content of soil (soil moisture)
%          4. water content of upper reservoir of flow routing routine 
%          5. water content of lower reservoir of flow routing routine
% FLUXES = time series of simulated fluxes (all in mm/Dt)    - matrix (T,8)
%          1: refreezing
%          2: snowmelt
%          3: actual evapotranspiration
%          4: recharge (water flux from soil moisture accounting module
%             to flow routing module)
%          5: percolation (water flux from upper to lower reservoir of the 
%             flow routing module)
%          6: runoff from upper reservoir
%          7: runoff from lower reservoir
%
% References: 
%
% Seibert,J.(1997)."Estimation of Parameter Uncertainty in the HBV Model".
% Nordic Hydrology.28(4/5).247-262.
% 
% Kollat,J.B.,Reed,P.M.,Wagener,T.(2012)."When are multiobjective 
% calibration trade-offs in hydrologic models meaningful?". Water resources
% research,VOL.48,W03520,doi:10.1029/2011WR011534.

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/ 

% Comments:
% * Model components: snow routine (optional)- soil moisture, upper zone
% and lower routine - flow routing routine

N = length(prec);
STATES=nan(N+1,5);
FLUXES=nan(N,7);

[P,STATES(:,1:2),FLUXES(:,1:2)] = snow_routine(temp,prec,param(1:4)) ;
[Q_sim,STATES(:,3:5),FLUXES(:,3:7)] = hbv_sim(P,ept,param(5:13),Case);

f=zeros(1,6);

Qs=Q_sim(warmup+1:end);
Qo=flow(warmup+1:end);

N=length(Qs);

N_67=floor(N*0.67);
N_33=floor(N*0.33);
Qs_sort=sort(Qs);
Qo_sort=sort(Qo);


lambda=0.3;
Zs=((1+Qs).^lambda-1)/lambda;
Zo=((1+Qo).^lambda-1)/lambda;

f(1)= mean(abs(Qs-Qo)); %AME
f(2)= 1 - sum((Qs - Qo).^2)/sum((mean(Qo) - Qo).^2); %NSE
f(3) = abs(mean(Qs - Qo)) ; %BIAS
f(4)=sqrt(mean((Zs - Zo).^2)); % TRMSE
f(5) = abs((Qs_sort(N_33)-Qs_sort(N_67))/(Qo_sort(N_33)-Qo_sort(N_67))-1)*100;  %SFDCE
f(6) = sqrt(mean((Qs - Qo).^2));%RMSE

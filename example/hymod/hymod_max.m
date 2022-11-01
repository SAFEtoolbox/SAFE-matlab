function [y,Q_sim,STATES,FLUXES] = hymod_max(x,rain,evap,warmup)
%
% This function runs the rainfall-runoff Hymod model
% and returns the maximum flow in the simulated time series
%
% [y,Q_sim,STATES,FLUXES] = hymod_max(param,rain,evap,warmup)
% 
% Input:
% param = vector of model parameters (Smax,beta,alfa,Rs,Rf)  - vector (1,5)
%  rain = time series of rainfall                            - vector (T,1)
%  evap = time series of potential evaporation               - vector (T,1)
% warmup = number of time steps for model warm-up            - scalar
%
% Output:
%      y = maximum flow over the simulation horizon          - scalar
%  Q_sim = time series of simulated flow                     - vector (T,1)
% STATES = time series of simulated storages (all in mm)     - matrix (T,5)
% FLUXES = time series of simulated fluxes (all in mm/Dt)    - matrix (T,4)
%
% See also hymod_sim about the model parameters, simulated variables,
% and references.

% y = hymod_max( x )
%
% x = model parameters
% y = maximum flow over the simulation horizon

M = 5 ; % number of model parameters
x = x(:);
if ~isnumeric(x); error('input argument ''param'' must be numeric'); end
if length(x)~=M; error('input argument ''param'' must have %d components',M); end

[Q_sim,STATES,FLUXES] = hymod_sim(x,rain,evap) ;

y = max(Q_sim(warmup+1:end)) ; 




function [y,Q_sim,STATES,FLUXES] = hymod_nse(x,rain,evap,flow,warmup)
%
% This function runs the rainfall-runoff Hymod model
% and returns the associated Nash-Sutcliffe Efficiency
%
% [y,Q_sim,STATES,FLUXES] = hymod_nse(param,rain,evap,flow,warmup)
% 
% Input:
% param = vector of model parameters (Smax,beta,alfa,Rs,Rf)  - vector (1,5)
%  rain = time series of rainfall                            - vector (T,1)
%  evap = time series of potential evaporation               - vector (T,1)
%  flow = time series of observed flow                       - vector (T,1)
% warmup = number of time steps for model warm-up            - scalar
%
% Output:
%      y = Nash-Sutcliffe Efficiency                         - scalar
%  Q_sim = time series of simulated flow                     - vector (T,1)
% STATES = time series of simulated storages (all in mm)     - matrix (T,5)
% FLUXES = time series of simulated fluxes (all in mm/Dt)    - matrix (T,8)
%
% See also hymod_sim about the model parameters, simulated variables,
% and references.


M = 5 ; % number of model parameters
x = x(:);
if ~isnumeric(x); error('input argument ''param'' must be numeric'); end
if length(x)~=M; error('input argument ''param'' must have %d components',M); end

[Q_sim,STATES,FLUXES] = hymod_sim(x,rain,evap) ;

Qs = Q_sim(warmup+1:end);
Qo = flow(warmup+1:end);

y = NSE(Qs', Qo');





function [zt, action] = integrate_bd(z0, T, direction, potential_and_gradient, beta, D, dt)
% Integrate Brownian dynamics trajectory from given initial configuration.
% 
% zt = integrate_bd(z0, nsteps, direction, potential, gradient, beta, D, dt)
%
% ARGUMENTS
%
% z0 (Kx1) - initial configuration
% T - total trajectory length (including initial given point)
% direction - +1 if forward in time, -1 if backward
% potential_and_gradient (function that takes Kx1 and returns [scalar, Kx1]) - potential and gradient
% beta - inverse temperature
% D - diffusion constant
% dt - timestep
%
% RETURNS
%
% zt (K x nsteps) - trajectory starting at z0.

% Determine dimensionality of configuration space.
[K,L] = size(z0);
if (L ~= 1)
  error('z0 must be Kx1 vector');
end

% Generate random variables.
R = randn(K,T);   %%% random variables, unit variance

% Initialize storage for trajectory.
zt = zeros(K,T);
zt(:,1) = z0;

% Initialize action.
[potential, gradient] = potential_and_gradient(z0);
action = beta * potential;

% Integrate using Brownian dynamics.
for t = 2:T
  % compute force
  [potential, gradient] = potential_and_gradient(zt(:,t-1));  
  force = - gradient;
  % compute displacement
  dx = sqrt(2*D*dt)*R(:,t) + beta*dt*D*force;
  % update position
  zt(:,t) = zt(:,t-1) + dx;
  % update action
  action = action + sum(R(:,t).^2/2);
end

% Return the trajectory.
return;



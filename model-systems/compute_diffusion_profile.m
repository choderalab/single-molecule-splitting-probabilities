function [xbins, Dx] = compute_diffusion_profile(x_t, x_A, x_B, dt)
% Estimate diffusion profile.
%
% Dx = compute_diffusion_profile(x_t, x_A, x_B)
%
% ARGUMENTS
%   x_t - trajectory
%   x_A - leftmost boundary
%   x_B - rightmost boundary

import timeseries

nbins = 50; % number of bins to compute diffusion profile along
xbins = linspace(x_A, x_B, nbins); % positions at which to estimate diffusion profile (FIXME)
dx = (x_B - x_A) / nbins; % width of bin to compute diffusion profile

tau_min = 50;
tau_max = 100;
tskip = 1;

Dx = zeros(1,nbins); % diffusion constant

for i = 1:nbins
  i
  % Compute time-correlation function.
  %C_t = timeseries.diffusion_pande(x_t, xbins(i), dx, tau_max, tskip)';
  C_t = timeseries.diffusion(x_t, xbins(i), dx, tau_max)';

  % Compute diffusion coefficient.
  tvec = (1:tau_max) * dt;
  indices = tau_min:tau_max; % indices for fit
%  P = polyfit(tvec(indices), log(C_t(indices)), 1);
%  Dx(i) = -dt / P(1); % diffusion constant


  % DEBUG PLOTTING
%  clf;
%  plot(tvec, log(C_t), 'k.');
%  hold on;
%  plot(tvec(indices), log(C_t(indices)), 'r.');
%  plot(tvec, polyval(P, tvec), 'r-');

  P = polyfit(tvec(indices), C_t(indices), 1);
  Dx(i) = 0.5 * P(1) / dt; % diffusion constant
%  plot(tvec, C_t, 'k.');
%  hold on;
%  plot(tvec(indices), C_t(indices), 'r.');
%  plot(tvec, polyval(P, tvec), 'r-');
%  pause

end


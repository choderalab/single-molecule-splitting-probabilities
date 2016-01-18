% Estimate diffusion profile.
% NOTE: This must have already been read in via 'analyze.m'.

import timeseries

dx_diffusion = 5.0; % instead of bin width, use larger width for estimating diffusion
tau_min = 1000;
tau_max = 4000;
tskip = 100;

D_i = zeros(1,nbins); % diffusion constant

for i = 1:nbins
  i
  % Compute time-correlation function.
  C_t = timeseries.diffusion_pande(x_t, xbins(i), dx_diffusion, tau_max, tskip)';
  % Compute diffusion coefficient.
  tvec = (1:tau_max) * dt;
  indices = tau_min:tau_max; % indices for fit
  P = polyfit(tvec(indices), log(C_t(indices)), 1);
  D_i(i) = -dt / P(1); % diffusion constant
end

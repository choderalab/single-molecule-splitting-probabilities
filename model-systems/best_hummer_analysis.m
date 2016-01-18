% Perform Bayesian analysis of Best and Hummer on trajectory data to estimate PMF, rates, and diffusion constant.

% Limits of data to use for analysis.
% NOTE: We can bring these in a bit
%xmin = min(x_t);
%xmax = max(x_t);
%nbins = 10; % number of bins to discretize reaction coordinate into

%xmin = 525;
%xmax = 555;
%nbins = 60;

% Lag time for estimating rates.
%tau = 200; % lag time, in number of samples

coordinate = 'x';
if coordinate == 'x'
  % x projection
  xmin = -2; xmax = +5.5; 
  x_A = -0.7; x_B = 5; % kT = 1
  x_t = xt;
  tau = 10;
  nbins = 50;
  dt = 0.01;
elseif coordinate == 'y'
  % y projection
  xmin = -2; xmax = +5.5; 
  x_A = -0.1; x_B = 5.5; % kT = 1
  xmin = x_A; xmax = x_B;
  x_t = yt;
  tau = 10;
  nbins = 50;
  dt = 0.01;
elseif coordinate == 'q'
  % q projection
  xmin = -3; xmax = +3; 
  x_A = -2.5; x_B = +2.5; % kT = 1
  xmin = x_A; xmax = x_B;
  x_t = qt;
  tau = 10;
  nbins = 50;
  dt = 0.01;
end

% Timeseries length.
T = length(x_t);

% Discretize trajectory into bins.
bin_edges = linspace(xmin, xmax, nbins+1); % bin edges
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2.0; % bin centers
%SMALL = (xmax-xmin)/100.0; % a small value
%bin_edges(end) = max(x_t) + SMALL; % move rightmost edge to contain all data
%bin_edges(1) = min(x_t) - SMALL; % move leftmost edge to contain all data
bin_widths = bin_edges(2:end) - bin_edges(1:end-1); % bin widths
[N_i, s_t] = histc(x_t, bin_edges);
N_i = N_i(1:end-1); % last entry will be zero

% Compute MLE of stationary probabilities.
p_i_mle = N_i / sum(N_i); % equilibrium populations

% Compute transition count matrix.
javaaddpath .
import bayesian.*
dx = (xmax - xmin)/nbins;
Nij = bayesian.compute_restricted_transition_counts(x_t, xmin, dx, nbins, tau);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate committor probabilities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Estimating committor probabilities...');
javaaddpath .
import timeseries.*
pAx_i = timeseries.observed_splitting(x_t, bin_edges, xmin, xmax);
pBx_i = 1 - pAx_i;
disp('Done.');

% Compute better rate matrix guess from transition matrix.
Tij = Nij;
for i = 1:nbins
  Tij(i,:) = Nij(i,:) / sum(Nij(i,:));
end
Kij = rate_matrix_guess(Tij, dt*tau);

% Estimate effective number of statistically independent samples using bin lifetimes.
Kij_mle = Kij;
lifetimes_i = -1.0 ./ diag(Kij) / dt; % estimated lifetimes
for i = 1:nbins
  Nij(i,:) = Nij(i,:) / lifetimes_i(i);
end
Neff = sum(sum(Nij));
disp(sprintf('%.1f effective samples', Neff));

% Sample from rate matrix.
nequil = 10; % number of equilibration samples to discard
nsamples = 20 + nequil; % number of samples to produce
Kij_samples = zeros(nsamples,nbins,nbins);
p_i_samples = zeros(nsamples,nbins);
F_i_samples = zeros(nsamples,nbins);
D_i_samples = zeros(nsamples,nbins);
D0_samples = zeros(nsamples);
s_i_samples = zeros(nsamples,nbins);
Fs_i_samples = zeros(nsamples,nbins);
pAx_pmf_i_samples = zeros(nsamples,nbins);
pBx_pmf_i_samples = zeros(nsamples,nbins);
for sample = 1:nsamples
  disp(sprintf('sample %d / %d', sample, nsamples));
  
  % Sample rate matrix.
  Kij = hummer_rate_matrix_update(Kij, Nij, tau*dt);
  %Kij = hummer_rate_matrix_update_njp(Kij, Nij, tau*dt);

  % Store sample.
  Kij_samples(sample,:,:) = Kij(:,:);

  % Compute stationary probability of this sample.  
  Tij = expm(Kij * dt * tau);
  options = struct('disp', 0);
  [V,D] = eigs(Tij',1,'LR',options);
  p_i = (V / sum(V))';

  % Compute PMF estimate.
  F_i = - log(p_i / dx);
  
  % Compute diffusion constant at intermediate bin indices.
  D_i_midpoint = zeros(1,nbins-1);
  for i = 1:(nbins-1)
    D_i_midpoint(i) = dx^2 * Kij(i,i+1) * sqrt(p_i(i) / p_i(i+1)); % hummer suggestion
  end  
  % Interpolate or extrapolate diffusion constant at bin indices.
  D_i = 0*F_i;
  D_i(1) = D_i_midpoint(1);
  D_i(end) = D_i_midpoint(end);
  D_i(2:end-1) = 0.5*(D_i_midpoint(1:end-1) + D_i_midpoint(2:end));

  % Store samples.
  p_i_samples(sample,:) = p_i;
  F_i_samples(sample,:) = F_i;
  D_i_samples(sample,:) = D_i;

  % Compute remapped coordinates.
  D0 = (dx*sum(D_i.^(-1/2))).^(-2); % effective diffusion constant
  s_i = dx*cumsum((D0./D_i).^(1/2)); % rescaled s-coordinate %%% FIXME
  Fs_i = F_i + 0.5 * log(D0./D_i); % corrected free energy in s-coordinate
  Fs_i = Fs_i - log(sum(exp(-Fs_i)));

  % Store samples.
  D0_samples(sample) = D0;
  s_i_samples(sample,:) = s_i;
  Fs_i_samples(sample,:) = Fs_i;

  % Compute Pfold from PMF.
  pAx_pmf_i = zeros(1,nbins);
  pBx_pmf_i = zeros(1,nbins);
  for i = 1:nbins
    pBx_pmf_i(i) = sum(exp(+Fs_i(1:i))) / sum(exp(+Fs_i(1:end)));
    pAx_pmf_i(i) = 1 - pBx_pmf_i(i);
  end

  % Store samples.
  pAx_pmf_i_samples(sample,:) = pAx_pmf_i;
  pBx_pmf_i_samples(sample,:) = pAx_pmf_i;

end

% Discard initial samples
p_i_samples = p_i_samples(nequil+1:end,:);
F_i_samples = F_i_samples(nequil+1:end,:);
D_i_samples = D_i_samples(nequil+1:end,:);
D0_samples = D0_samples(nequil+1:end);
s_i_samples = s_i_samples(nequil+1:end,:);
Fs_i_samples = Fs_i_samples(nequil+1:end,:);
pAx_pmf_i_samples = pAx_pmf_i_samples(nequil+1:end,:);
pBx_pmf_i_samples = pBx_pmf_i_samples(nequil+1:end,:);

% DEBUG: Plot
clf;
hold on;
plot(s_i, pAx_i, 'r-', 'linewidth', 2)
plot(xpoints, pAx_pmf_i_samples, 'k-');
legend('PMF', 'observed');
xlabel('s');
ylabel('p_A(x)');
axis([0 1 0 1]);

filename = sprintf('best-hummer-%s.eps', coordinate);
exportfig(gcf, filename, 'width', 10, 'height', 2.5, 'color', 'cmyk', 'resolution', 600);
system(sprintf('epstopdf %s', filename));

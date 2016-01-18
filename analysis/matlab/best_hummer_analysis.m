% Perform Bayesian analysis of Best and Hummer on trajectory data to estimate PMF, rates, and diffusion constant.

% Limits of data to use for analysis.
% NOTE: We can bring these in a bit
xmin = min(x_t);
xmax = max(x_t);
nbins = 25; % number of bins to discretize reaction coordinate into

xmin = 525;
xmax = 555;
nbins = 60;

% Timeseries length.
T = length(x_t);

% Lag time for estimating rates.
tau = 200; % lag time, in number of samples

% Estimate statistical inefficiency.
disp('Computing statistical inefficiency...');
g = statistical_inefficiency_mex(x_t, x_t);
Neff = length(x_t) / g;
disp(sprintf('g = %.1f, Neff = %.1f', g, Neff));

% Discretize trajectory into bins.
bin_edges = linspace(xmin, xmax, nbins+1); % bin edges
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2.0; % bin centers
SMALL = (xmax-xmin)/100.0; % a small value
bin_edges(end) = max(x_t) + SMALL; % move rightmost edge to contain all data
bin_edges(1) = min(x_t) - SMALL; % move leftmost edge to contain all data
bin_widths = bin_edges(2:end) - bin_edges(1:end-1); % bin widths
[N_i, s_t] = histc(x_t, bin_edges);
N_i = N_i(1:end-1); % last entry will be zero

% Explore statistical inefficiency of transitions.
gi = zeros(1,nbins);
gij = zeros(nbins,nbins);
for i = 1:nbins
  A_t = (s_t(1:T)==i) + 0.0; 
  gi(i) = statistical_inefficiency_mex(A_t, A_t);
  disp(sprintf('gi(%3d) = %.1f', i, gi(i)));
%  
%  for j = 1:nbins
%    A_t = ((s_t(1:T-tau)==i) & (s_t(1+tau:T)==j)) + 0.0;
%    gij(i,j) = statistical_inefficiency_mex(A_t, A_t);
%    disp(sprintf('gij(%3d,%3d) = %.1f', i, j, gij(i,j)));
%  end
end


% Compute MLE of stationary probabilities.
p_i_mle = N_i / sum(N_i); % equilibrium populations

% Compute transition count matrix.
javaaddpath .
import bayesian.*
Nij = bayesian.compute_transition_counts(s_t, tau);
% Correct for statistical inefficiency.
%Nij = Nij / g;
for i = 1:nbins
  Nij(i,:) = Nij(i,:) / gi(i);
end
Neff = sum(sum(Nij));
disp(sprintf('%.1f effective samples', Neff));

% Compute rate matrix guess.
Kij = (real(logm(Nij))/tau);
Kij = zeros(nbins,nbins) + diag(abs(diag(Kij,1)),1) + diag(abs(diag(Kij,-1)),-1);
Kij = Kij - diag(sum(Kij,2),0);

% Sample from rate matrix.
nsamples = 100; % number of samples to produce
p_i_samples = zeros(nsamples,nbins);
F_i_samples = zeros(nsamples,nbins);
D_i_samples = zeros(nsamples,nbins-1);
for sample = 1:nsamples
  disp(sprintf('sample %d / %d', sample, nsamples));
  
  % Sample rate matrix.
  Kij = hummer_rate_matrix_update(Kij, Nij, tau);

  % Compute stationary probability of this sample.  
  Tij = expm(Kij * tau);
  [V,D] = eigs(Tij',1,'LR');
  p_i = (V / sum(V))';

  % Compute PMF estimate.
  F_i = - log(p_i ./ bin_widths);
  
  % Compute diffusion constant at intermediate bin indices.
  D_i = zeros(1,nbins-1);
  for i = 1:(nbins-1)
    D_i(i) = bin_widths(i)^2 * Kij(i,i+1) * sqrt(p_i(i) / p_i(i+1));
  end  
  D_i

  % Store samples.
  p_i_samples(sample,:) = p_i;
  F_i_samples(sample,:) = F_i;
  D_i_samples(sample,:) = D_i;
end



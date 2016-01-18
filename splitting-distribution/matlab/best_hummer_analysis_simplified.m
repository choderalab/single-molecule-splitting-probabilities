function best_hummer_analysis_simplified(coordinate, x_t, xmin, xmax)
% Perform Bayesian analysis of Best and Hummer on trajectory data to estimate PMF and diffusion constant,
% as well as remapping from x to s where diffusion constant is uniform.

% Parameters.

tau = 1; % Lag time for estimating rates.
nbins = 50;
dt = 0.01;

% Timeseries length.
T = length(x_t);

% Discretize trajectory into bins.
bin_edges = linspace(xmin, xmax, nbins+1); % bin edges
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2.0; % bin centers
bin_widths = bin_edges(2:end) - bin_edges(1:end-1); % bin widths
[N_i, s_t] = histc(x_t, bin_edges);
N_i(nbins) = N_i(nbins) + N_i(end); N_i = N_i(1:nbins); 

% Compute MLE of stationary probabilities.
p_i_mle = N_i / sum(N_i); % equilibrium populations

% Compute statistical inefficiency.
javaaddpath .
import timeseries.*
g = timeseries.statistical_inefficiency(x_t, x_t);

% Compute transition count matrix.
javaaddpath .
import bayesian.*
dx = (xmax - xmin)/nbins;
Nij = bayesian.compute_restricted_transition_counts(x_t, xmin, dx, nbins, tau);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate committor probabilities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Estimating committor probabilities from timeseries data...');
javaaddpath .
import timeseries.*
pAx_i = timeseries.observed_splitting(x_t, bin_edges, xmin, xmax);
pBx_i = 1 - pAx_i;
disp('Done.');

% Compute estimate of rate matrix from transition matrix.
Tij = Nij;
for i = 1:nbins
  Tij(i,:) = Nij(i,:) / sum(Nij(i,:));
end
Kij = rate_matrix_guess(Tij, dt*tau);

% Estimate effective number of statistically independent samples using bin lifetimes.
%Kij_mle = Kij;
%lifetimes_i = -1.0 ./ diag(Kij) / dt; % estimated lifetimes (in sampling periods)
%lifetimes_i
%for i = 1:nbins
%  Nij(i,:) = Nij(i,:) / lifetimes_i(i);
%end

disp('Computing statistical inefficiencies...');
g_i = ones(nbins,1);
for i = 1:nbins
  A_t = (x_t > bin_edges(i)) & (x_t < bin_edges(i+1)) + 0.0;
  g_i(i) = timeseries.statistical_inefficiency(A_t, A_t);  
  Nij(i,:) = Nij(i,:) / g_i(i);
  disp(sprintf('N_i[%5d] = %8.1f', i, sum(Nij(i,:))));
end
g_i
Neff = sum(sum(Nij));
disp(sprintf('%.1f effective samples', Neff));

% Sample from rate matrix.
nequil = 20; % number of equilibration samples to discard
nsamples = 20 + nequil; % number of samples to produce
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
    pBx_pmf_i(i) = sum(exp(+F_i(1:i)) .* D_i(1:i).^(-1)) / sum(exp(+F_i(1:end)) .* D_i(1:end).^(-1));
    pAx_pmf_i(i) = 1 - pBx_pmf_i(i);
  end

  % Compute Pfold from remapped PMF.
  pAx_pmf_remapped_i = zeros(1,nbins);
  pBx_pmf_remapped_i = zeros(1,nbins);
  for i = 1:nbins
    pBx_pmf_remapped_i(i) = sum(exp(+Fs_i(1:i))) / sum(exp(+Fs_i(1:end)));
    pAx_pmf_remapped_i(i) = 1 - pBx_pmf_i(i);
  end

%  % Compare
%  disp('comparison:');
%  pAx_pmf_i
%  pAx_pmf_remapped_i
%  pAx_pmf_i - pAx_pmf_remapped_i

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot PMF along observed coordinate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;

hold on;
plot(bin_centers, F_i_samples, 'r-');
%legend('PMF', 'observed');
xlabel(sprintf('%s / nm', coordinate));
xlabel(coordinate);
ylabel('F / kT');
%set(gca, 'XTick', [0 1]);
%set(gca, 'YTick', [0 0.2 0.4 0.6 0.8 1]);
set(gca, 'XTick', [525 535 545 555]);
set(gca, 'XTick', [xmin xmax]);
box on;
%axis([0 1 0 1]);
axis square
axis([xmin xmax 0 max(max(F_i_samples))]);

% Set font properties for all axes.
fontsize = 7;
fontname = 'Arial';
fontweight = 'bold';
set(findall(0, '-property', 'FontSize'), 'FontSize', fontsize);
set(findall(0, '-property', 'FontName'), 'FontName', fontname);
set(findall(0, '-property', 'FontWeight'), 'FontWeight', fontweight);

addpath exportfig
filename = sprintf('pmf-%s.eps', coordinate);
exportfig(gcf, filename, 'width', 1.75, 'height', 1.75, 'color', 'cmyk', 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

% bw
filename = sprintf('pmf-%s-bw.eps', coordinate);
exportfig(gcf, filename, 'width', 1.75, 'height', 1.75, 'color', 'bw', 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot diffusion constant along observed coordinate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
clf;

plot(bin_centers, D_i_samples, 'r-');
%semilogy(bin_centers, D_i_samples, 'r-'); % log scale
%legend('PMF', 'observed');
xlabel(sprintf('%s / nm', coordinate));
xlabel(coordinate);
ylabel('D');
%set(gca, 'XTick', [0 1]);
%set(gca, 'YTick', [0 0.2 0.4 0.6 0.8 1]);
set(gca, 'XTick', [525 535 545 555]);
set(gca, 'XTick', [xmin xmax]);
%set(gca, 'YTick', [10^0 10^1 10^2 10^3 10^4]); % log scale
box on;
%axis([0 1 0 1]);
axis square
axis([xmin xmax 0 max(max(D_i_samples))]);

% Set font properties for all axes.
fontsize = 7;
fontname = 'Arial';
fontweight = 'bold';
set(findall(0, '-property', 'FontSize'), 'FontSize', fontsize);
set(findall(0, '-property', 'FontName'), 'FontName', fontname);
set(findall(0, '-property', 'FontWeight'), 'FontWeight', fontweight);

addpath exportfig
filename = sprintf('diffusion-%s.eps', coordinate);
exportfig(gcf, filename, 'width', 1.75, 'height', 1.75, 'color', 'cmyk', 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

% bw
filename = sprintf('diffusion-%s-bw.eps', coordinate);
exportfig(gcf, filename, 'width', 1.75, 'height', 1.75, 'color', 'bw', 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

figure(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot splitting probabilities in remapped s-coordinate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;

hold on;
plot(bin_centers, pAx_pmf_i_samples, 'r-');
plot(bin_centers, pAx_i, 'k--', 'linewidth', 2)
%legend('PMF', 'observed');
xlabel(sprintf('%s / nm', coordinate));
xlabel(coordinate);
ylabel('p_A');
set(gca, 'XTick', [525 535 545 555]);
set(gca, 'XTick', [xmin xmax]);
set(gca, 'YTick', [0 0.2 0.4 0.6 0.8 1]);
box on;
axis([xmin xmax 0 1]);
axis square

% Set font properties for all axes.
fontsize = 7;
fontname = 'Arial';
fontweight = 'bold';
set(findall(0, '-property', 'FontSize'), 'FontSize', fontsize);
set(findall(0, '-property', 'FontName'), 'FontName', fontname);
set(findall(0, '-property', 'FontWeight'), 'FontWeight', fontweight);

addpath exportfig
filename = sprintf('best-hummer-%s.eps', coordinate);
exportfig(gcf, filename, 'width', 1.75, 'height', 1.75, 'color', 'cmyk', 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

filename = sprintf('best-hummer-%s-bw.eps', coordinate);
exportfig(gcf, filename, 'width', 1.75, 'height', 1.75, 'color', 'bw', 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

return

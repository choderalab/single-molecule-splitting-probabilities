% Perform splitting analysis on Woodside single-molecule force-spectroscopy data.
% NOTE: This must have already been read in via 'read_data.m'.

% PARAMETERS

x_A = 525.0; % absorbing boundary (nm)
x_B = 555.0; % absorbing boundary (nm)
nbins = 100; % number of histogram bins

%x_A = 533.68;
%x_B = 549.561;
%nbins = 100;

% Estimate statistical inefficiency for timeseries.
disp('Computing statistical inefficiency...');
javaaddpath .;
import timeseries.*;
g = timeseries.statistical_inefficiency(x_t, x_t);
T = length(x_t);
Neff = T/g;
disp(sprintf('g = %.1f; Neff = %.1f', g, Neff));

% Find min and max
xmin = min(x_t);
xmax = max(x_t);

% Discretize trajectory into bins.
bin_edges = linspace(xmin, xmax, nbins+1); % bin edges
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2.0; % bin centers
bin_widths = bin_edges(2:end) - bin_edges(1:end-1); % bin widths
[N_i, s_t] = histc(x_t, bin_edges);
N_i(nbins) = N_i(nbins) + N_i(end); N_i = N_i(1:nbins); 

% Estimate population histograms.
xbins = bin_centers;
Tx = N_i'; % samples per bin
Nx = Tx / g; % number of uncorrelated samples/bin
N = sum(Nx); % total number of uncorrelated samples
dx = xbins(2) - xbins(1); % compute bin width
Px = Nx / N; % probability per bin estimate
lightgray = 0.8 * [1 1 1]; % color gray for shading error bars
darkgray = 0.6 * [1 1 1]; % color gray for shading error bars
% compute statistical uncertainties
dPx = sqrt(Px .* (1.0 - Px)) ./ sqrt(N);
% Incorporate Jacobian.
Px = Px * dx;
dPx = dPx * dx;

% Estimate potential of mean force.
fx = - log(Px);  fx = fx - min(fx);
dfx = dPx ./ Px; 

% Compute analytical committor from PMF
disp('Computing analytical committor...');
pAx_pmf = zeros(size(Px),'double');
pBx_pmf = zeros(size(Px),'double');
% find bins in allowed region
allowed_bins = find((xbins > x_A) & (xbins < x_B));
% find bins in absorbed regions
A_bins = find(xbins <= x_A);
B_bins = find(xbins >= x_B);
% mark bins in absorbed regions
pAx_pmf(A_bins) = 1.0;
pAx_pmf(B_bins) = 0.0;
pBx_pmf(A_bins) = 0.0;
pBx_pmf(B_bins) = 1.0;
% find bins just inside absorbed region
x_A_bin = allowed_bins(1) - 1;
x_B_bin = allowed_bins(end) + 1;
% compute splitting probability by integration
for i = allowed_bins
  pBx_pmf(i) = sum(exp(+fx(x_A_bin:i))) / sum(exp(+fx(x_A_bin:x_B_bin)));
end
% compute pB
pAx_pmf = 1.0 - pBx_pmf;

% Estimate committor probabilities.
disp('Estimating committor probabilities...');
javaaddpath .
import timeseries.*
retval = timeseries.observed_splitting_with_error(x_t, bin_edges, x_A, x_B)';
pAx = retval(:,1)'
dpAx = retval(:,2)'
pBx = 1 - pAx;
dpBx = dpAx;

save hairpin_analysis.mat
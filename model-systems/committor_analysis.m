function [xbins, pAx_pmf] = committor_analysis(x_t, x_A, x_B)

nbins = 40;

% Estimate statistical inefficiency.
disp('Computing statistical inefficiency...');
g = statistical_inefficiency_mex(x_t, x_t);
disp(sprintf('g = %.1f', g));

% Find min and max
xmin = min(x_t);
xmax = max(x_t);

% Estimate population histograms.
[Tx, xbins] = hist(x_t,nbins); % samples per bin
Nx = Tx / g; % number of uncorrelated samples/bin
N = sum(Nx); % total number of uncorrelated samples
dx = xbins(2) - xbins(1); % compute bin width
Px = Nx / N; % probability per bin estimate
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
% find bins just inside absorbed region
x_A_bin = allowed_bins(1) - 1;
x_B_bin = allowed_bins(end) + 1;
% compute splitting probability by integration
for i = allowed_bins
  pBx_pmf(i) = sum(exp(+fx(x_A_bin:i))) / sum(exp(+fx(x_A_bin:x_B_bin)));
end
% compute pB
pAx_pmf = 1.0 - pBx_pmf;

return

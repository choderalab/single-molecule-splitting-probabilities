% Perform splitting analysis on Woodside single-molecule force-spectroscopy data.
% NOTE: This must have already been read in via 'read_data.m'.

% PARAMETERS

x_A = 525.0; % absorbing boundary (nm)
x_B = 555.0; % absorbing boundary (nm)
nbins = 50; % number of histogram bins

%x_A = 533.68;
%x_B = 549.561;
%nbins = 100;

% Estimate statistical inefficiency for timeseries.
disp('Computing statistical inefficiency...');
g = statistical_inefficiency_mex(x_t, x_t);
T = length(x_t);
Neff = T/g;
disp(sprintf('g = %.1f; Neff = %.1f', g, Neff));

% Find min and max
xmin = min(x_t);
xmax = max(x_t);

% Estimate population histograms.
[Tx, xbins] = hist(x_t,nbins); % samples per bin
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
A_t = zeros(size(x_t), 'int8') - 1;
B_t = zeros(size(x_t), 'int8') - 1;
% mark all spans in the absorbed regions
disp('Marking absorbed spans...');
indices = find(x_t <= x_A);
A_t(indices) = 1;
B_t(indices) = 0;
indices = find(x_t >= x_B);
B_t(indices) = 1;
A_t(indices) = 0;
% iterate to mark
disp('Iterating to mark...')
unassigned = find(A_t == -1);
while (length(unassigned) > 0)
  % find A-committed states
  indices = find(A_t(unassigned + 1) == 1);
  indices = unassigned(indices);
  A_t(indices) = 1;
  B_t(indices) = 0;
  % find B-committed states
  indices = find(B_t(unassigned + 1) == 1);
  indices = unassigned(indices);
  B_t(indices) = 1;
  A_t(indices) = 0;
  % find remaining unassigned states
  unassigned = find(A_t == -1);
end
% compute committor probabilities
disp('Computing committor probabilities...');
pAx = zeros(size(Px),'double');
dpAx = zeros(size(Px),'double');
pBx = zeros(size(Px),'double');
dpBx = zeros(size(Px),'double');
for i = 1:nbins
  indices = find((x_t >= xbins(i) - dx/2) & (x_t < xbins(i) + dx/2)); % indices of samples in this bin
  Ti = length(indices); % number of correlated samples in this bin
  Ni = Ti / g; % number of uncorrelated samples in this bin
  disp(sprintf('Tx(%d) = %d, Ti = %d', i, Tx(i), Ti));
  pAx(i) = mean(A_t(indices)); % committor for A
  pBx(i) = mean(B_t(indices)); % committor for B
  dpAx(i) = sqrt(pAx(i) * (1.0 - pAx(i))) / sqrt(Ni);
  dpBx(i) = sqrt(pAx(i) * (1.0 - pAx(i))) / sqrt(Ni);
end

save analysis.mat

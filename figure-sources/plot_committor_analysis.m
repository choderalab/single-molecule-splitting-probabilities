function plot_committor_analysis(x_t, x_A, x_B, position, coordinate_label, time_label, dt, xx, Fx)
% Plot committor analysis of the trajectory.
%
% plot_committor_analysis(x_t, x_A, x_B, position, coordinate_label, time_label, dt, xx, Fx)
%
% ARGUMENTS
%   x_t - trajectory 
%   x_A - left absorbing boundary
%   x_B - right absorbing boundary
%   position - ignored
%   coordinate_label - text label for the coordinate x
%   time_label - text label for time coordinate t
%   dt - time interval between samples for coordinate t
%   xx - x coordinates at which PMF Fx is evaluted
%   Fx - PMF cooresponding to xx
%
% RETURNS
%   xbins - bin centers
%   pAx - committor for state A
%   dpAx - uncertainy for pAx
%   pAx_pmf - committor from PMF

T = length(x_t);
nbins = 100;

times = [0:(T-1)] * dt; % sampling times for trajectory

%clf;
figure_handle = gcf;

% PROPERTIES

fontsize = 7;
fontname = 'Arial';
fontweight = 'bold';

%figurewidth = 7.5; % inches
%figureheight = 2.0; % inches
%set(figure_handle, 'PaperUnits', 'inches', 'PaperSize', [figurewidth figureheight], 'PaperPosition', [0 0 figurewidth figureheight]);
%set(figure_handle, 'Units', 'inches', 'Position', [0 0 figurewidth figureheight]);

axes_handles = [];
label_handles = [];

% Estimate statistical inefficiency.
disp('Computing statistical inefficiency...');
javaaddpath .;
import timeseries;
g = timeseries.statistical_inefficiency(x_t, x_t);
disp(sprintf('g = %.1f', g));
disp(sprintf('Neff = %.1f', T / g));

% Find min and max
xmin = min(x_t);
xmax = max(x_t);

% Determine bin edges and centers.
bin_edges = linspace(xmin, xmax, nbins+1);
dx = (xmax - xmin) / nbins;
xbins = bin_edges(1:nbins) + dx/2;

% Estimate population histograms.
Tx = histc(x_t, bin_edges); % samples per bin
Tx(nbins) = Tx(nbins) + Tx(nbins+1); Tx = Tx(1:nbins); % move from nbins+1 overflow bin to nbins
Nx = Tx / g; % number of uncorrelated samples/bin
N = sum(Nx); % total number of uncorrelated samples
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
pBx_pmf = ones(size(Px),'double');
% find bins in allowed region
allowed_bins = find((xbins > x_A) & (xbins < x_B));
% find bins in absorbed regions
A_bins = find(xbins <= x_A);
B_bins = find(xbins >= x_B);
% mark bins in absorbed regions
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

% Plot histogram.
ax = subplot('Position', [0.1 0.15 0.15 0.75]);
axes_handles = [axes_handles ax]; % store axis handle
%[Tx, xbins] = hist(x_t,nbins); % samples per bin
Nx = Tx / g; % number of uncorrelated samples/bin
N = sum(Nx); % total number of uncorrelated samples
Px = Nx / N; % probability per bin estimate
lightgray = 0.8 * [1 1 1]; % color gray for shading error bars
darkgray = 0.6 * [1 1 1]; % color gray for shading error bars
% compute statistical uncertainties
dPx = sqrt(Px .* (1.0 - Px)) ./ sqrt(N);
% plot absorbing boundaries

Px = Px * dx;
dPx = dPx * dx;
hold on;
fill([(Px-2*dPx) fliplr(Px+2*dPx)], [xbins fliplr(xbins)], lightgray, 'LineWidth', 1.0e-6); % 2 sigma
fill([(Px-dPx) fliplr(Px+dPx)], [xbins fliplr(xbins)], darkgray, 'LineWidth', 1.0e-6); % 1 sigma
plot(Px, xbins, 'k-', 'LineWidth', 2.0); % MLE
axis([0 max(Px+2*dPx) xmin xmax]);
set(ax, 'XTick', [0]);
ylabel(coordinate_label);
h = xlabel('p(x)'); label_handles = [label_handles h];
set(ax, 'XDir', 'reverse'); % flip direction of X plot

% Plot trajectory
ax = subplot('Position', [0.25 0.15 0.4 0.75]);
axes_handles = [axes_handles ax]; % store axis handle
plot(times, x_t, 'k.', 'MarkerSize', 0.5);
%plot(times, x_t, 'k-');
h = xlabel(time_label); label_handles = [label_handles h];
axes_limits = axis;
axis([0 T*dt xmin xmax]);
set(ax, 'YTick', []);
set(ax, 'XTick', [0]);

% Plot free energy
ax = subplot('Position', [0.65 0.15 0.15 0.75]);
axes_handles = [axes_handles ax]; % store axis handle
fx = - log(Px);  fx = fx - min(fx);
dfx = dPx ./ Px; 
hold on;
fill([(fx-2*dfx) fliplr(fx+2*dfx)], [xbins fliplr(xbins)], lightgray, 'LineWidth', 1.0e-6); % 2 sigma
fill([(fx-dfx) fliplr(fx+dfx)], [xbins fliplr(xbins)], darkgray, 'LineWidth', 1.0e-6); % 1 sigma
plot(fx, xbins, 'k-', 'LineWidth', 2.0); % MLE
axis([0 6 xmin xmax]);
h = xlabel(sprintf('F(%s) / kT', coordinate_label)); label_handles = [label_handles h];
set(ax, 'YTick', []);
set(ax, 'XTick', [0 2 4 6]);
% plot analytical PMF if available
if length(xx) > 0
  hold on;
  plot(Fx,xx,'b-');
end

% Plot splitting probabilities.
plot_height = (x_B - x_A) / (xmax - xmin) * 0.75;
ax = subplot('Position', [0.8 0.15+(0.75-plot_height)/2 0.15 plot_height]);
axes_handles = [axes_handles ax]; % store axis handle
hold on;
fill([(pAx-2*dpAx) fliplr(pAx+2*dpAx)], [xbins fliplr(xbins)], lightgray, 'LineWidth', 1.0e-6); % 2 sigma % FIX ME
fill([(pAx-dpAx) fliplr(pAx+dpAx)], [xbins fliplr(xbins)], darkgray, 'LineWidth', 1.0e-6); % 1 sigma
plot(pAx, xbins, 'k-', 'LineWidth', 2.0); % MLE
axis([0 1 x_A x_B]);
h = xlabel('p_A'); label_handles = [label_handles h];
set(ax, 'YTick', []);
set(ax, 'XDir', 'reverse'); % flip direction of X plot

% Set font properties for all axes.
for ax = axes_handles
  set(ax, 'FontSize', fontsize);
  set(ax, 'FontName', fontname);
  set(ax, 'FontWeight', fontweight);
end

for ax = label_handles
  set(ax, 'FontSize', fontsize);
  set(ax, 'FontName', fontname);
  set(ax, 'FontWeight', fontweight);
end

% Plot pmf-based splitting
plot(pAx_pmf, xbins, 'r--', 'LineWidth', 2.0); % PMF-based

% Draw x_A and x_B lines.
ax = axes('Position', [0.1 0.15 0.70 0.75], 'Visible', 'off');
linewidth = 1;
hold on;
plot([0 1], x_A * [1 1], 'b:', 'LineWidth', linewidth);
plot([0 1], x_B * [1 1], 'b:', 'LineWidth', linewidth);
axis(ax, [0 1 xmin, xmax])

% Set font properties for all axes.
set(findall(0, '-property', 'FontSize'), 'FontSize', fontsize);
set(findall(0, '-property', 'FontName'), 'FontName', fontname);
set(findall(0, '-property', 'FontWeight'), 'FontWeight', fontweight);

% Draw boxes on all axes.
set(findall(0, '-property', 'box'), 'box', 'on');

return

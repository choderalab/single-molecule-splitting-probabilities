% Plot the trajectory

times = [0:(T-1)] * dt; % sampling times for trajectory

clf;
figure_handle = gcf;

% PROPERTIES

fontsize = 7;
fontname = 'Arial';
fontweight = 'bold';


figurewidth = 7.5; % inches
figureheight = 2.0; % inches
set(figure_handle, 'PaperUnits', 'inches', 'PaperSize', [figurewidth figureheight], 'PaperPosition', [0 0 figurewidth figureheight]);
set(figure_handle, 'Units', 'inches', 'Position', [0 0 figurewidth figureheight]);

%nbins = 100;

axes_handles = [];
label_handles = [];

xmin = min(x_t);
xmax = max(x_t);

% plot histogram
ax = subplot('Position', [0.1 0.15 0.15 0.75]);
axes_handles = [axes_handles ax]; % store axis handle
[Tx, xbins] = hist(x_t,nbins); % samples per bin
Nx = Tx / g; % number of uncorrelated samples/bin
N = sum(Nx); % total number of uncorrelated samples
dx = xbins(2) - xbins(1); % compute bin width
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
ylabel('extension (nm)');
h = xlabel('p(x)'); label_handles = [label_handles h];
set(ax, 'XDir', 'reverse'); % flip direction of X plot

% plot trajectory
ax = subplot('Position', [0.25 0.15 0.4 0.75]);
axes_handles = [axes_handles ax]; % store axis handle
plot(times, x_t, 'k.', 'MarkerSize', 0.5);
plot(times, x_t, 'k-');
h = xlabel('time (s)'); label_handles = [label_handles h];
axes_limits = axis;
axis([0 T*dt xmin xmax]);
set(ax, 'YTick', []);

% plot free energy
ax = subplot('Position', [0.65 0.15 0.15 0.75]);
axes_handles = [axes_handles ax]; % store axis handle
fx = - log(Px);  fx = fx - min(fx);
dfx = dPx ./ Px; 
hold on;
fill([(fx-2*dfx) fliplr(fx+2*dfx)], [xbins fliplr(xbins)], lightgray, 'LineWidth', 1.0e-6); % 2 sigma
fill([(fx-dfx) fliplr(fx+dfx)], [xbins fliplr(xbins)], darkgray, 'LineWidth', 1.0e-6); % 1 sigma
plot(fx, xbins, 'k-', 'LineWidth', 2.0); % MLE
axis([0 10 xmin xmax]);
h = xlabel('f(x) (kT)'); label_handles = [label_handles h];
set(ax, 'YTick', []);
set(ax, 'XTick', [0 2 4 6 8]);

% plot splitting probabilities
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

% plot analytical splitting probability
line(pAx_pmf, xbins, 'Color', 'r', 'LineStyle', '-', 'LineWidth', linewidth); % MLE

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

% write PDF
figurename = 'figure';
print('-depsc', sprintf('%s.eps', figurename)); % create EPS figure
unix(sprintf('epstopdf %s.eps', figurename)); % convert to PDF
unix(sprintf('ps2pdf13 -dPDFSETTINGS=/prepress %s.pdf %s-embedded.pdf', figurename, figurename)); % embed fonts in PDF




% Plot the trajectory

times = [0:(T-1)] * dt; % sampling times for trajectory

clf;
figure_handle = gcf;

% PROPERTIES

fontsize = 7;
fontname = 'Arial';
fontweight = 'bold';

x_A = 525.0; % absorbing boundary (nm)
x_B = 555.0; % absorbing boundary (nm)

figurewidth = 7.5; % inches
figureheight = 2.0; % inches
set(figure_handle, 'PaperUnits', 'inches', 'PaperSize', [figurewidth figureheight], 'PaperPosition', [0 0 figurewidth figureheight]);
set(figure_handle, 'Units', 'inches', 'Position', [0 0 figurewidth figureheight]);

nbins = 100;

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

% plot splitting probabilities
ax = subplot('Position', [0.8 0.15 0.15 0.75]);
axes_handles = [axes_handles ax]; % store axis handle
A_t = zeros(size(x_t), 'int8');
B_t = zeros(size(x_t), 'int8');
% mark all spans in the absorbed regions
A_t(find(x_t <= x_A)) = 1;
B_t(find(x_t >= x_B)) = 1;
% locate all spans in the transition region
disp('Locating all absorbed trajectories...')
start_index = find((x_t > x_A) & (x_t < x_B), 1); % find first sample in transition region
while ((start_index >= 1) & (start_index < T))
  start_index

  % find first sample in absorbed region
  first_absorbed_index = find((x_t(start_index:end) <= x_A) | (x_t(start_index:end) >= x_B), 1) + start_index - 1; 
  first_absorbed_index

  % record whether these samples are committed to A or B
  if (x_t(first_absorbed_index) <= x_A)
    A_t(start_index:first_absorbed_index) = 1;
  else
    B_t(start_index:first_absorbed_index) = 1;
  end

  % find next sample in transition region
  start_index = find((x_t(first_absorbed_index:end) > x_A) & (x_t(first_absorbed_index:end) < x_B), 1) + first_absorbed_index - 1; 
end
% compute committor probabilities
disp('Computing committor probabilities...');
pAx = zeros(size(Px),'double');
dpAx = zeros(size(Px),'double');
pBx = zeros(size(Px),'double');
dpBx = zeros(size(Px),'double');
for i = 1:nbins
  indices = find((x_t >= xbin(i) - dx/2) & (x_t < xbin(i) + dx/2)); % indices of samples in this bin
  Ti = length(indices); % number of correlated samples in this bin
  Ni = Ti / g; % number of uncorrelated samples in this bin
  disp(sprintf('Tx(%d) = %d, Ti = %d', i, Tx(i), Ti));
  pAx(i) = mean(A_t(indices)); % committor for A
  pBx(i) = mean(B_t(indices)); % committor for B
  dpAx(i) = sqrt(pAx(i) * (1.0 - pAx(i))) / sqrt(Ni);
  dpBx(i) = sqrt(pAx(i) * (1.0 - pAx(i))) / sqrt(Ni);
end
hold on;
fill([(pAx-2*dpAx) fliplr(pAx+2*dpAx)], [xbins fliplr(xbins)], lightgray, 'LineWidth', 1.0e-6); % 2 sigma
fill([(pAx-dpAx) fliplr(pAx+dpAx)], [xbins fliplr(xbins)], darkgray, 'LineWidth', 1.0e-6); % 1 sigma
plot(pAx, xbins, 'k-', 'LineWidth', 2.0); % MLE
axis([0 10 xmin xmax]);
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

% write PDF
figurename = 'figure';
print('-depsc', sprintf('%s.eps', figurename)); % create EPS figure
unix(sprintf('epstopdf %s.eps', figurename)); % convert to PDF
unix(sprintf('ps2pdf13 -dPDFSETTINGS=/prepress %s.pdf %s-embedded.pdf', figurename, figurename)); % embed fonts in PDF




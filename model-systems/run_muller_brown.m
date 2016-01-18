% Generate plots for simple two-dimensional model.

beta = 0.06; % inverse temperature

xmin = -1.5;
xmax = +1;
ymin = -0.5;
ymax = +2;

% Evaluate on grid, a row at a time.
npoints = 100;
xpoints = linspace(xmin, xmax, npoints);
ypoints = linspace(ymin, ymax, npoints);
Vxy = zeros(npoints, npoints);
for xindex = 1:npoints
  % Form row of same x point, all y points.
  q = [repmat(xpoints(xindex), 1, npoints); ypoints];
  [potential, gradient] = muller_brown(q);
  alpha = 0/beta;
  Vxy(:,xindex) = beta * potential - alpha * ypoints;
end

% Make contour plot.
clf;
contour_levels = linspace(-9, +9, 50);
contour(xpoints, ypoints, Vxy, contour_levels);
axis([xmin xmax ymin ymax]);
axis square;
colorbar;

% Simulate a brownian dynamics trajectory.
q0 = [0.5; 0];
T = 10000;
D = 1.0;
dt = 0.0001;
qt = integrate_bd(q0, T, +1, @(q) muller_brown(q), beta, D, dt);
hold on;
plot(qt(1,:), qt(2,:), 'k-');


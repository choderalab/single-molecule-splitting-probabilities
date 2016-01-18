% Generate trajectory and plots for simple two-dimensional Rhee-Pande model.

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kT = 5; % reduced temperature
beta = 1/kT; % inverse temperature

xmin = -3;
xmax = +7;
ymin = -2;
ymax = +7;
qmin = -3.5;
qmax = +3.5;

fontsize = 7;
fontname = 'Arial';
fontweight = 'bold';
markersize = 5;
columnwidth = 3.375;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate potential on grid, a row at a time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npoints = 50;
xpoints = linspace(xmin, xmax, npoints);
ypoints = linspace(ymin, ymax, npoints);
Vxy = zeros(npoints, npoints);
for xindex = 1:npoints
  % Form row of same x point, all y points.
  q = [repmat(xpoints(xindex), 1, npoints); ypoints];
  [potential, gradient] = rhee_pande(q);
  Vxy(:,xindex) = potential;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate a brownian dynamics trajectory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate.
q0 = [4; 1];
T = 1000000;
D = 1.0;
dt = 0.01;
qt = integrate_bd(q0, T, +1, @(q) rhee_pande(q), beta, D, dt);

% Extract coordinate trajectories.
xt = qt(1,:);
yt = qt(2,:);
qt = (xt-yt) / sqrt(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store trajectory data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save trajectory.mat;



% Generate plots for simple two-dimensional model.

clear;

kT = 5; % reduced temperature
%kT = 1; % reduced temperature
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

% Evaluate on grid, a row at a time.
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

% Simulate a brownian dynamics trajectory.
q0 = [4; 1];
T = 1000000;
D = 1.0;
dt = 0.01;
qt = integrate_bd(q0, T, +1, @(q) rhee_pande(q), beta, D, dt);

% Extract coordinate trajectories.
xt = qt(1,:);
yt = qt(2,:);
qt = (xt-yt) / sqrt(2);

% Make contour plot.
figure(1);
clf;

subplot(2,2,1);
contour_levels = [0 5 10 15 20 25];
[C,H] = contour(xpoints, ypoints, Vxy, contour_levels);
%hclabel = clabel(C,H,'labelspacing',1000); 
axis([xmin xmax ymin ymax]);
axis square;
%colorbar;
hold on;
xlabel('x');
ylabel('y');

% Label q coordinate.
x0 = 4;
y0 = 1;
hold on; 
plot([x0 x0-3], [y0 y0+3], 'r-', 'linewidth', 2);
text(x0+0.5, y0, 'q');



% pmf along y
subplot(2,2,2);
Px = zeros(npoints,1);
for xindex = 1:npoints
  Vx = Vxy(:,xindex);
  Px(xindex) = sum(exp(-beta*Vx));
end
Px = Px / sum(Px);
Fx = - log(Px);
Fx = Fx - min(Fx);
plot(xpoints, Fx, 'k-');
axis([xmin xmax 0 3.5]);
axis square;
xlabel('y');
ylabel('F(y) / kT');
set(gca,'YTick',[0 1 2 3]);

% pmf along x
subplot(2,2,3);
Py = zeros(npoints,1);
for yindex = 1:npoints
  Vy = Vxy(yindex,:);
  Py(yindex) = sum(exp(-beta*Vy));
end
Py = Py / sum(Py);
Fy = - log(Py);
Fy = Fy - min(Fy);
plot(ypoints, Fy, 'k-');
axis([ymin ymax 0 3.5]);
axis square;
xlabel('x');
ylabel('F(x) / kT');
set(gca,'YTick',[0 1 2 3]);

% pmf along q
subplot(2,2,4);
dx = xpoints(2) - xpoints(1);
qpoints = linspace(qmin, qmax, npoints);
Pq = zeros(npoints,1);
for qindex = 1:npoints
  q = qpoints(qindex);
  %[Vq, gradient] = rhee_pande(sqrt(2)*[xpoints - q; xpoints + q]);
  [Vq, gradient] = rhee_pande([xpoints; xpoints - sqrt(2)*q]);
  Pq(qindex) = sum(exp(-beta*Vq));
end
Pq = Pq / sum(Pq);
Fq = - log(Pq);
Fq = Fq - min(Fq);
plot(qpoints, Fq, 'k-');
axis([qmin qmax 0 3.5]);
axis square;
xlabel('q');
ylabel('F(q) / kT');
set(gca,'YTick',[0 1 2 3]);

% Set font properties for all axes.
set(findall(0, '-property', 'FontSize'), 'FontSize', fontsize);
set(findall(0, '-property', 'FontName'), 'FontName', fontname);
set(findall(0, '-property', 'FontWeight'), 'FontWeight', fontweight);
set(findall(0, '-property', 'box'), 'box', 'on');

% Resize contour labels.
%set(hclabel, 'FontSize', 4.5);

%print -depsc rhee_pande_pmf.eps
exportfig(gcf, 'rhee_pande_pmf.eps', 'width', columnwidth, 'height', columnwidth, 'color', 'cmyk', 'fontmode', 'scaled');
system('epstopdf rhee_pande_pmf.eps');

% DEBUG
%return

%% Make contour plot.
figure(2);
clf;

Tplot = 100000;
plotskip = 10;

%subplot('position', [0.05 0.05 0.30 0.9]);
%contour_levels = [0 5 10 15 20 25 30];
%contour(xpoints, ypoints, Vxy, contour_levels);
%axis([xmin xmax ymin ymax]);
%axis square;
%%colorbar;
%hold on;
%plot(xt, yt, 'k-');
%xlabel('x');
%ylabel('y');

subplot('position', [0.1 0.65 0.85 0.25])
plot(xt(plotskip:plotskip:Tplot), 'k.', 'markersize', markersize);
set(gca,'XTick',[]);
ylabel('x');
axis([0 Tplot/plotskip xmin xmax]);

subplot('position', [0.1 0.35 0.85 0.25])
plot(yt(plotskip:plotskip:Tplot), 'k.', 'markersize', markersize);
set(gca,'XTick',[]);
ylabel('y');
axis([0 Tplot/plotskip ymin ymax]);

subplot('position', [0.1 0.05 0.85 0.25])
plot(qt(plotskip:plotskip:Tplot), 'k.', 'markersize', markersize);
set(gca,'XTick',[]);
ylabel('q');
xlabel('time');
axis([0 Tplot/plotskip qmin qmax]);

% Set font properties for all axes.
set(findall(0, '-property', 'FontSize'), 'FontSize', fontsize);
set(findall(0, '-property', 'FontName'), 'FontName', fontname);
set(findall(0, '-property', 'FontWeight'), 'FontWeight', fontweight);
set(findall(0, '-property', 'box'), 'box', 'on');

%print -depsc rhee_pande_trajectory.eps
exportfig(gcf, 'rhee_pande_trajectory.eps', 'width', columnwidth, 'height', columnwidth, 'color', 'cmyk', 'fontmode', 'scaled');
system('epstopdf rhee_pande_trajectory.eps');

save trajectory.mat;

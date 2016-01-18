% Generate plots for simple two-dimensional Rhee-Pande model.

clear;

% Read trajectory data and parameters.
load trajectory.mat; 

% Add figure assistance packages.
addpath exportfig
addpath epscombine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot 2D potential and PMFs along x, y, and q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
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

% Plot pmf along y
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

% Plot pmf along x
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

% Plot pmf along q
subplot(2,2,4);
dx = xpoints(2) - xpoints(1);
qpoints = linspace(qmin, qmax, npoints);
Pq = zeros(npoints,1);
for qindex = 1:npoints
  q = qpoints(qindex);
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

% Write plot.
filename = 'rhee-pande-pmf.eps';
exportfig(gcf, filename, 'width', columnwidth, 'height', columnwidth, 'color', 'cmyk', 'fontmode', 'scaled');
system(sprintf('epstopdf %s', filename));

% Write plot.
filename = 'rhee-pande-pmf-bw.eps';
exportfig(gcf, filename, 'width', columnwidth, 'height', columnwidth, 'color', 'bw', 'fontmode', 'scaled');
system(sprintf('epstopdf %s', filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot trajectory projected along x, y, and q.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
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

filename = 'rhee-pande-trajectory.eps';
exportfig(gcf, filename, 'width', columnwidth, 'height', columnwidth, 'color', 'cmyk', 'fontmode', 'scaled');
system(sprintf('epstopdf %s', filename));

filename = 'rhee-pande-trajectory-bw.eps';
exportfig(gcf, filename, 'width', columnwidth, 'height', columnwidth, 'color', 'bw', 'fontmode', 'scaled');
system(sprintf('epstopdf %s', filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate committor plots in x.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
clf;

position = [1 1 1 1];
coordinate_label = 'x';
time_label = 'time';
%x_A = -2; x_B = +5.5; % kT = 5
x_A = -0.7; x_B = 5; % kT = 1
x_t = xt;
dt = 0.1;
filename = 'rhee-pande-full_committor-x.eps';
plot_committor_analysis(x_t, x_A, x_B, position, coordinate_label, time_label, dt, xpoints, Fx);
exportfig(gcf, 'zbuffer', 'color', 'cmyk', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'zbuffer', 'resolution', 300);
exportfig(gcf, 'painter', 'color', 'cmyk', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'painters');
epscombine('painter_t.eps','zbuffer.eps',filename);
system(sprintf('epstopdf %s', filename));

% bw
filename = 'rhee-pande-full_committor-x-bw.eps';
exportfig(gcf, 'zbuffer', 'color', 'bw', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'zbuffer', 'resolution', 300);
exportfig(gcf, 'painter', 'color', 'bw', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'painters');
epscombine('painter_t.eps','zbuffer.eps',filename);
system(sprintf('epstopdf %s', filename));

clf;
plot_committor_only(x_t, x_A, x_B, position, coordinate_label, time_label, dt, xpoints, Fx);
filename = 'rhee-pande-committor-x.eps';
exportfig(gcf, filename, 'color', 'cmyk', 'width', 1.5, 'height', 1.5, 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

% bw
filename = 'rhee-pande-committor-x-bw.eps';
exportfig(gcf, filename, 'color', 'bw', 'width', 1.5, 'height', 1.5, 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate committor plots in y.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
clf;

position = [1 1 1 1];
coordinate_label = 'y';
time_label = 'time';
%x_A = -1.3; x_B = +6.5; % kT = 5
x_A = -0.1; x_B = 5.5; % kT = 1
x_t = yt;
dt = 0.1;
filename = 'rhee-pande-full_committor-y.eps';
plot_committor_analysis(x_t, x_A, x_B, position, coordinate_label, time_label, dt, ypoints, Fy);
exportfig(gcf, 'zbuffer', 'color', 'cmyk', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'zbuffer', 'resolution', 300);
exportfig(gcf, 'painter', 'color', 'cmyk', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'painters');
epscombine('painter_t.eps','zbuffer.eps',filename);
system(sprintf('epstopdf %s', filename));

filename = 'rhee-pande-full_committor-y-bw.eps';
exportfig(gcf, 'zbuffer', 'color', 'bw', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'zbuffer', 'resolution', 300);
exportfig(gcf, 'painter', 'color', 'bw', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'painters');
epscombine('painter_t.eps','zbuffer.eps',filename);
system(sprintf('epstopdf %s', filename));

clf;
plot_committor_only(x_t, x_A, x_B, position, coordinate_label, time_label, dt, xpoints, Fx);
filename = 'rhee-pande-committor-y.eps';
exportfig(gcf, filename, 'color', 'cmyk', 'width', 1.5, 'height', 1.5, 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

% bw
filename = 'rhee-pande-committor-y-bw.eps';
exportfig(gcf, filename, 'color', 'bw', 'width', 1.5, 'height', 1.5, 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate committor plots in q.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
clf;

position = [1 1 1 1];
coordinate_label = 'q';
time_label = 'time';
%x_A = -3.3; x_B = +3.3; % kT = 5
x_A = -2.5; x_B = +2.5; % kT = 1
x_t = qt;
dt = 0.1;
filename = 'rhee-pande-full_committor-q.eps';
plot_committor_analysis(x_t, x_A, x_B, position, coordinate_label, time_label, dt, qpoints, Fq);
exportfig(gcf, 'zbuffer', 'color', 'cmyk', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'zbuffer', 'resolution', 300);
exportfig(gcf, 'painter', 'color', 'cmyk', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'painters');
epscombine('painter_t.eps','zbuffer.eps',filename);
system(sprintf('epstopdf %s', filename));

% bw
filename = 'rhee-pande-full_committor-q-bw.eps';
exportfig(gcf, 'zbuffer', 'color', 'bw', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'zbuffer', 'resolution', 300);
exportfig(gcf, 'painter', 'color', 'bw', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'painters');
epscombine('painter_t.eps','zbuffer.eps',filename);
system(sprintf('epstopdf %s', filename));

clf;
plot_committor_only(x_t, x_A, x_B, position, coordinate_label, time_label, dt, xpoints, Fx);
filename = 'rhee-pande-committor-q.eps';
exportfig(gcf, filename, 'color', 'cmyk', 'width', 1.5, 'height', 1.5, 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

filename = 'rhee-pande-committor-q-bw.eps';
exportfig(gcf, filename, 'color', 'bw', 'width', 1.5, 'height', 1.5, 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate committor plots for diffusion-remapped coordinate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x projection
xmin = -2; xmax = +5.5; 
x_A = -0.7; x_B = 5; % kT = 1
best_hummer_analysis_simplified('x', xt, x_A, x_B)

% y projection
xmin = -2; xmax = +5.5; 
x_A = -0.1; x_B = 5.5; % kT = 1
best_hummer_analysis_simplified('y', yt, x_A, x_B)

% q projection
xmin = -3; xmax = +3; 
x_A = -2.5; x_B = +2.5; % kT = 1
best_hummer_analysis_simplified('q', qt, x_A, x_B)


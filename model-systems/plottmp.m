% Generate committor analysis plots for model system.

load trajectory

addpath exportfig
addpath epscombine

%% Plot x

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

clf;
plot_committor_only(x_t, x_A, x_B, position, coordinate_label, time_label, dt, xpoints, Fx);
filename = 'rhee-pande-committor-x.eps';
exportfig(gcf, filename, 'color', 'cmyk', 'width', 1.5, 'height', 1.5, 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

%% Plot y

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

clf;
plot_committor_only(x_t, x_A, x_B, position, coordinate_label, time_label, dt, xpoints, Fx);
filename = 'rhee-pande-committor-y.eps';
exportfig(gcf, filename, 'color', 'cmyk', 'width', 1.5, 'height', 1.5, 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

%% Plot q

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

clf;
plot_committor_only(x_t, x_A, x_B, position, coordinate_label, time_label, dt, xpoints, Fx);
filename = 'rhee-pande-committor-q.eps';
exportfig(gcf, filename, 'color', 'cmyk', 'width', 1.5, 'height', 1.5, 'renderer', 'painters');
system(sprintf('epstopdf %s', filename));

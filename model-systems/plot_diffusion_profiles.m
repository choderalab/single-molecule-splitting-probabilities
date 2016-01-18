% Calculate and plot diffusion profiles.

%figure(6)
clf;

position = [1 1 1 1];
coordinate_label = 'x';
time_label = 'time';
x_A = -2; x_B = +5.5; % kT = 5
%x_A = -0.7; x_B = 5; % kT = 1
x_t = xt;
dt = 0.1;
%filename = 'rhee_pande_x.eps';
[xx, Dx] = compute_diffusion_profile(x_t, x_A, x_B, dt);
subplot(3,1,1); % DEBUG
plot(xx, Dx, '.');
%exportfig(gcf, filename, 'width', 10, 'height', 2.5, 'color', 'cmyk');
%system(sprintf('epstopdf %s', filename));

%figure(7)
%clf;
%
position = [1 1 1 1];
coordinate_label = 'y';
time_label = 'time';
x_A = -1.3; x_B = +6.5; % kT = 5
x_A = -0.1; x_B = 5.5; % kT = 1
x_t = yt;
dt = 0.1;
%filename = 'rhee_pande_y.eps';
%plot_committor_analysis(x_t, x_A, x_B, position, coordinate_label, time_label, dt, ypoints, Fy);
%exportfig(gcf, filename, 'width', 10, 'height', 2.5, 'color', 'cmyk');
%system(sprintf('epstopdf %s', filename));
subplot(3,1,2); % DEBUG
[xx, Dx] = compute_diffusion_profile(x_t, x_A, x_B, dt);
plot(xx, Dx, '.');


%
%figure(8)
%clf;
%
position = [1 1 1 1];
coordinate_label = 'q';
time_label = 'time';
x_A = -3.3; x_B = +3.3; % kT = 5
x_A = -2.5; x_B = +2.5; % kT = 1
x_t = qt;
dt = 0.1;
filename = 'rhee_pande_q.eps';
%plot_committor_analysis(x_t, x_A, x_B, position, coordinate_label, time_label, dt, qpoints, Fq);
%exportfig(gcf, filename, 'width', 10, 'height', 2.5, 'color', 'cmyk');
%system(sprintf('epstopdf %s', filename));
subplot(3,1,3); % DEBUG
[xx, Dx] = compute_diffusion_profile(x_t, x_A, x_B, dt);
plot(xx, Dx, '.');

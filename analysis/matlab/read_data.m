% Read Woodside DNA optical-trap dataset.

clear;

% PARAMETERS

data_filename = '../../original-data/optical-trap/Hpin37_trajectory.txt';
dt = 40e-6; % sampling interval (s)

% Read trajectory into x_t.
disp('Reading dataset...')
x_t = textread(data_filename); % x_t[t] is the extension coordinate of sample t (nm)

% Determine trajectory length.
T = length(x_t); % number of time-correlated samples
disp(sprintf('%d samples read (%.3f s sampled at rate %e s/sample', T, T*dt, dt));


% Read Woodside DNA optical-trap dataset.

clear;

% PARAMETERS

data_filename = '../original-data/optical-trap/Hpin37_trajectory.txt';
dt = 40e-6; % sampling interval (s)

% Read trajectory into x_t.
disp('Reading dataset...')
x_t = textread(data_filename); % x_t[t] is the extension coordinate of sample t (nm)

% Determine trajectory length.
T = length(x_t); % number of time-correlated samples
disp(sprintf('%d samples read (%.3f s sampled at rate %e s/sample', T, T*dt, dt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate committor plots for diffusion-remapped coordinate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x_A = 525.0; % absorbing boundary (nm)
x_B = 555.0; % absorbing boundary (nm)

best_hummer_analysis_simplified('extension', x_t, x_A, x_B)


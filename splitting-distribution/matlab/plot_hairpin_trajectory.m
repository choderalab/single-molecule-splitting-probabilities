% Plot the trajectory

load hairpin_analysis.mat;

times = [0:(T-1)] * dt; % sampling times for trajectory

clf;
figure_handle = gcf;

% PROPERTIES

fmax = 6.0; % maximum free energy

fontsize = 7;
fontname = 'Arial';
fontweight = 'bold';

figurewidth = 7.5; % inches
figureheight = 2.0; % inches

linewidth = 1.5;
edgecolor = 'k';
edgewidth = 0.1;

set(figure_handle, 'PaperUnits', 'inches', 'PaperSize', [figurewidth figureheight], 'PaperPosition', [0 0 figurewidth figureheight]);
set(figure_handle, 'Units', 'inches', 'Position', [0 0 figurewidth figureheight]);

axes_handles = [];
label_handles = [];

% plot histogram
ax = subplot('Position', [0.07 0.15 0.15 0.75]);
axes_handles = [axes_handles ax]; % store axis handle
hold on;
fill([(Px-2*dPx) fliplr(Px+2*dPx)], [xbins fliplr(xbins)], lightgray, 'EdgeColor', edgecolor, 'LineWidth', edgewidth); % 2 sigma
fill([(Px-dPx) fliplr(Px+dPx)], [xbins fliplr(xbins)], darkgray, 'EdgeColor', edgecolor, 'LineWidth', edgewidth); % 1 sigma
plot(Px, xbins, 'k-', 'LineWidth', linewidth); % MLE
axis([0 max(Px+2*dPx) xmin xmax]);
set(ax, 'XTick', [0]);
ylabel('extension / nm');
h = xlabel('p(x)'); label_handles = [label_handles h];
set(ax, 'XDir', 'reverse'); % flip direction of X plot

% plot trajectory
ax = subplot('Position', [0.22 0.15 0.43 0.75]);
axes_handles = [axes_handles ax]; % store axis handle
%plot(times, x_t, 'k.', 'MarkerSize', 0.25);
%semilogx(times, x_t, 'k.', 'MarkerSize', 2.0);
markersize = 0.5;
plot(times, x_t, 'k.', 'markersize', markersize);
h = xlabel('time / s'); label_handles = [label_handles h];
axes_limits = axis;
axis([0 T*dt xmin xmax]);
set(ax, 'YTick', []);

% plot free energy
ax = subplot('Position', [0.65 0.15 0.15 0.75]);
axes_handles = [axes_handles ax]; % store axis handle
hold on;
fill([(fx-2*dfx) fliplr(fx+2*dfx)], [xbins fliplr(xbins)], lightgray, 'EdgeColor', edgecolor, 'LineWidth', edgewidth); % 2 sigma
fill([(fx-dfx) fliplr(fx+dfx)], [xbins fliplr(xbins)], darkgray, 'EdgeColor', edgecolor, 'LineWidth', edgewidth); % 1 sigma
plot(fx, xbins, 'k-', 'LineWidth', linewidth); % MLE
axis([0 fmax xmin xmax]);
h = xlabel('pmf / kT'); label_handles = [label_handles h];
set(ax, 'YTick', []);
set(ax, 'XTick', [0 2 4 6 8]);

% plot splitting probabilities
plot_height = (x_B - x_A) / (xmax - xmin) * 0.75;
ax1 = subplot('Position', [0.8 0.15+(x_A - xmin)/(xmax - xmin)*0.75 0.15 plot_height]);
hold on;
% plot pAx

% store axis handles
axes_handles = [axes_handles ax1];

fill([(pAx-2*dpAx) fliplr(pAx+2*dpAx)], [xbins fliplr(xbins)], lightgray, 'EdgeColor', edgecolor, 'LineWidth', edgewidth,  'Parent', ax1); % 2 sigma 

fill([(pAx-dpAx) fliplr(pAx+dpAx)], [xbins fliplr(xbins)], darkgray, 'EdgeColor', edgecolor, 'LineWidth', edgewidth,  'Parent', ax1); % 1 sigma

line(pAx, xbins, 'Color', 'k', 'LineWidth', linewidth, 'LineStyle', '-', 'Parent', ax1); % MLE

% plot analytical x_A
line(pAx_pmf, xbins, 'Color', 'r', 'LineStyle', '--', 'LineWidth', linewidth, 'Parent', ax1); % MLE

axis(ax1, [0 1 x_A x_B]);

% label axes
h = xlabel(ax1, 'p_A'); label_handles = [label_handles h];

set(ax1, 'YTick', []);

set(ax1, 'XTick', [0 0.5]);

set(ax1, 'XDir', 'reverse'); % flip direction of X plot

% plot x_A and x_B
text(-0.05, x_A - (x_B-x_A)*0.025, 'x_A', 'FontSize', fontsize, 'FontName', fontname, 'FontWeight', fontweight);
text(-0.05, x_B - (x_B-x_A)*0.025, 'x_B', 'FontSize', fontsize, 'FontName', fontname, 'FontWeight', fontweight);

%format_ticks(ax2, {}, {'$x_A$', '$x_B$'}, [], [x_A, x_B]);

% Set font properties for all axes.
for ax = axes_handles
  set(ax, 'FontSize', fontsize);
  set(ax, 'FontName', fontname);
  set(ax, 'FontWeight', fontweight);
  set(ax, 'Box', 'on'); % draw plot box
end

for ax = label_handles
  set(ax, 'FontSize', fontsize);
  set(ax, 'FontName', fontname);
  set(ax, 'FontWeight', fontweight);
end

% draw x_A and x_B lines
ax = axes('Position', [0.1 0.15 0.7 0.75], 'Visible', 'off');
linewidth = 1;
hold on;
plot([0 1], x_A * [1 1], 'k-', 'LineWidth', linewidth);
plot([0 1], x_B * [1 1], 'k-', 'LineWidth', linewidth);
axis(ax, [0 1 xmin, xmax])

% write PDF
figurename = 'dna-hairpin';
filename = sprintf('%s.eps', figurename);
%print('-depsc', sprintf('%s.eps', figurename)); % create EPS figure
%unix(sprintf('epstopdf %s.eps', figurename)); % convert to PDF
%unix(sprintf('ps2pdf13 -dPDFSETTINGS=/prepress %s.pdf %s-embedded.pdf', figurename, figurename)); % embed fonts in PDF
addpath exportfig
addpath epscombine
exportfig(gcf, 'zbuffer', 'color', 'cmyk', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'zbuffer', 'resolution', 300);
exportfig(gcf, 'painter', 'color', 'cmyk', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'painters');
epscombine('painter_t.eps','zbuffer.eps',filename);
system(sprintf('epstopdf %s', filename));

% bw
figurename = 'dna-hairpin-bw';
exportfig(gcf, 'zbuffer', 'color', 'bw', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'zbuffer', 'resolution', 300);
exportfig(gcf, 'painter', 'color', 'bw', 'width', 10, 'height', 2.5, 'separatetext', 1, 'renderer', 'painters');
epscombine('painter_t.eps','zbuffer.eps',filename);
system(sprintf('epstopdf %s', filename));



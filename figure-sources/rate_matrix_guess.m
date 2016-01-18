function Kij = rate_matrix_guess(Tij, tau)
% Guess the neighbor-only rate matrix that produced the given transition matrix for specified lag time.
% 
% Kij = rate_matrix_guess(Tij, tau)

disp('Computing guess of rate matrix...');

% Compute guess of forward and reverse rates.

% Naive.
kf = diag(Tij,1) / tau;
kr = diag(Tij,-1) / tau;

% For q, tau=1
%kf = 5*diag(Tij,1) / tau;
%kr = 5*diag(Tij,-1) / tau;

% For q, tau=2
%kf = 10*diag(Tij,1) / tau;
%kr = 10*diag(Tij,-1) / tau;

% Others?
%D = diag(Tij,0);
%kf = 3*(1-D(1:end-1)) / tau;
%kr = 3*(1-D(2:end)) / tau;

% Define objective function.
% objective = |Tij - expm(Kij)|
rate_matrix = @(X) diag(X(:,1),1)+diag(X(:,2),-1)-diag(sum(diag(X(:,1),1)+diag(X(:,2),-1),2));
objective = @(X) norm(Tij-expm(rate_matrix(X)*tau));

% DEBUG: Starting objective
objective([kf kr])

% DEBUG: Plot transition matrix before optimization.
clf;
subplot(2,2,1);
imagesc(Tij);
colorbar;
axis square
subplot(2,2,2);
imagesc(expm(rate_matrix([kf kr])*tau))
colorbar
axis square

% Optimize rate matrix.
options = optimset('fminsearch');
options.Display = 'notify';
options.MaxFunEvals = 1000;
Xopt = fminsearch(objective, [kf kr], options);

% Compute rate matrix.
Kij = rate_matrix(Xopt);

% DEBUG: Plot transition matrix after optimization.
subplot(2,2,3);
imagesc(Tij);
colorbar;
axis square
subplot(2,2,4);
imagesc(expm(Kij*tau))
colorbar
axis square
print -djpeg -r300 rate-matrix-guess.jpg

% DEBUG: Final objective
objective(Xopt)

return
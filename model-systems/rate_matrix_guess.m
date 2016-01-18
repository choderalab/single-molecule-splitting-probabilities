function Kij = rate_matrix_guess(Tij, tau)
% Guess the neighbor-only rate matrix that produced the given transition matrix for specified lag time.
% 
% Kij = rate_matrix_guess(Tij, tau)

disp('Computing guess of rate matrix...');

% Compute guess of forward and reverse rates.
kf = diag(Tij,1) / tau;
kr = diag(Tij,-1) / tau;

% Define objective function.
% objective = |Tij - expm(Kij)|
rate_matrix = @(X) diag(X(:,1),1)+diag(X(:,2),-1)-diag(sum(diag(X(:,1),1)+diag(X(:,2),-1),2));
objective = @(X) norm(Tij-expm(rate_matrix(X)*tau));

% Optimize rate matrix.
options = optimset('fminsearch');
options.Display = 'notify';
options.MaxFunEvals = 1000;
Xopt = fminsearch(objective, [kf kr], options);

% Compute rate matrix.
Kij = rate_matrix(Xopt);

return
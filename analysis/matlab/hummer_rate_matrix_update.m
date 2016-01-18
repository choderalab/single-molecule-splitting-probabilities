function Kij = hummer_rate_matrix_update(Kij, Nij, tau)
% Produce a (correlated) sample of the transition matrix corresponding to a true rate matrix that satisfies detailed balance using scheme of Gerhard Hummer.
%
% Kij = hummer_rate_matrix_update(Kij, Nij, tau)
%
% ARGUMENTS
%  Kij (MxM matrix) - current rate matrix sample
%    This transition matrix must correspond to a valid rate matrix, such that Tij = expm(Kij * tau)
%  Nij (MxM matrix) - transition conditional count matrix
%    Nij(i,j) is number of observed counts terminating in j at time tau given that they initially started in i at time zero.
%  tau - evolution time to relate rate matrix to transition matrix through Tij = expm(Kij * tau)
%
% RETURN VALUES
%  Kij_new (MxM matrix) - updated sample of rate matrix, may be correlated with previous sample
%
% NOTES
%  The algorithm described in [1] and modified in [2] (especially with regard to the prior) is implemented.
%  Transitions between all states are allowed -- no constraints on the sparisty of Kij are imposed.
%
% REFERENCES
%  [1] Sriraman S, Kevrekidis IG, and Hummer G. Coarse master equation from Bayesian analysis of replica molecular dynamics simulations. JPC B 109(14):6479-6484, 2005.
%  [2] Buchete N-V and Hummer G. Coarse master equations for peptide folding dynamics. JPC B 112(19):6057-6069, 2008.
%
% TODO
%  * Check to be sure we don't need a correction for asymmetric proposal probability.
%  * Automatically adjust delta to provide ~ 50% acceptance before fixing delta.
%  * Automatically determine 'nsteps' to provide an uncorrelated sample.
%  * Speed up sampling with Java methods?
%  * Consider starting sampling from maximum-likelihood estimate if we're going to produce a completely decorrelated estimate.

% PARAMETERS
nsteps = 1000; % number of update steps (set to give approximately uncorrelated new sample)
delta = 1.0; % move width in logarithmic parameterization (set to give acceptance probability of ~ 50%)

% Check that Kij is real.
if ~isreal(Kij)
  Kij
  error('Kij is imaginary');
end

% Determine number of states.
M = size(Kij,1);

% Promote Nij to double.
Nij = double(Nij);

% Compute initial log-likelihood.
logL = log_likelihood(Kij, Nij, tau);

% Update sample using Hummer scheme.
naccept = 0; % number of accepted moves
for step = 1:nsteps
  % Choose an element to perturb.
  if (rand <= 0.5)
    i = randi(M-1);
    j = i + 1;
  else
    i = 1+randi(M-1);
    j = i - 1;
  end

  % Perturb this element in logarithm.
  Kij_new = zeros(M,M) + diag(diag(Kij,+1),+1) + diag(diag(Kij,-1),-1);
  Kij_new(i,j) = exp(log(Kij(i,j)) + (2.0*rand-1.0)*delta);
  Kij_new = Kij_new - diag(sum(Kij_new,2),0);
  
  if ~isreal(Kij_new)
    disp(sprintf('i = %d, j = %d', i, j));
    Kij(i,j)
    log(Kij(i,j))
    Kij_new(i,j)
    error('Proposed move is imaginary Kij');
  end

  % Compute new log-likelihood function.
  logL_new = log_likelihood(Kij_new, Nij, tau);
    
  % Compute change.
  delta_logL = logL_new - logL;
    
  % DEBUG
  %disp(sprintf('step %6d : logL_new = %16.1f : logL = %16.1f : delta_logL = %16.1f', step, logL_new, logL, delta_logL));
    
  % Accept or reject by Metropolis-Hastings procedure.
  if (rand < exp(delta_logL))
    % Store new parameters.
    Kij = Kij_new;
    logL = logL_new;
    naccept = naccept + 1;
  end

end

% Show summary statistics.
disp(sprintf('accepted %d / %d move attempts (%.1f %%)', naccept, nsteps, naccept / nsteps * 100));

return

%%%%%%%%%%%%%%

function logL = log_likelihood(Kij, Nij, tau)
% Compute log-likelihood.
%
% ARGUMENTS
%  Kij (MxM) - row-stochastic rate matrix - Kij(i,j) is rate from i -> j
%  Nij (MxM) - row transition count matrix - Nij(i,j) is number of transitions observed terminating in j a time tau after the system was initially in state i
%
% RETURN VALUE
%  logL - log-likelihood

Tij = expm(Kij * tau);
logTij = log(Tij);

% Compute log-likelihood
logL = sum(sum(Nij .* logTij));

return

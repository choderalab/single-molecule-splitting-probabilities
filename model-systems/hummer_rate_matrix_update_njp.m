function Kij = hummer_rate_matrix_update_njp(Kij, Nij, tau)
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
pdelta = 0.1; % move width in logarithm of probability
kdelta = 1.0; % move width in rate

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

% Compute stationary probability of this sample.  
Tij = expm(Kij * tau);
options = struct('disp', 0);
[V,D] = eigs(Tij',1,'LR', options);
p_i = (V / sum(V))';

% Update sample using Hummer scheme.
naccept = 0; % number of accepted moves
for step = 1:nsteps
  % Choose an element to perturb.
  if (rand <= 0.5)
    % Perturb one element of g_i, i = 1...(M-1)
    i = randi(M-1);
    p_i_new = p_i;
    p_i_new(i) = exp(log(p_i_new(i)) + (2.0*rand-1.0)*pdelta);
    % Reject if we exceed normalization.
    if(sum(p_i_new(1:(M-1))) > 1.0)
      % REJECT
      %disp('normalization failure')
      continue      
    end
    % Renormalize.
    p_i_new(M) = 1.0 - sum(p_i_new(1:(M-1)));
    % Check realness.
    if ~isreal(p_i_new)
      step
      disp('p_i_new is not real');
      p_i_new
      return
    end
    % Update Kij.
    Kij_new = Kij;
    if i > 1
      Kij_new(i,i-1) = Kij_new(i-1,i) * sqrt(p_i_new(i-1) / p_i_new(i));
    end
    Kij_new(M,M-1) = Kij_new(M-1,M) * sqrt(p_i_new(M-1) / p_i_new(M));
    
    % Correct diagonal elements.
    Kij_new = Kij_new - diag(diag(Kij_new));
    Kij_new = Kij_new - diag(sum(Kij_new,2),0);

  else
    % Choose an element to perturb.
    i = randi(M-1); % choose i from {1,2,...,M-1}
    Kij_new = Kij;
    p_i_new = p_i;
    Kij_new(i,i+1) = Kij_new(i,i+1) + kdelta * (2.0*rand-1.0);
    Kij_new(i+1,i) = Kij_new(i,i+1) * sqrt(p_i_new(i) / p_i_new(i+1));

    % Correct diagonal elements.
    Kij_new = Kij_new - diag(diag(Kij_new));
    Kij_new = Kij_new - diag(sum(Kij_new,2),0);     

    if any((Kij_new - diag(diag(Kij_new))) < 0.0)
      disp('Kij_new has negative elements.')
      continue
    end
  end

  if ~isreal(Kij_new)
    % reject
    disp('Kij is not real');
    continue;
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
    p_i = p_i_new;    
    logL = logL_new;
    naccept = naccept + 1;

%    % Update p_i
%    Tij = expm(Kij * tau);
%    options = struct('disp', 0);
%    [V,D] = eigs(Tij',1,'LR', options);
%    p_i = (V / sum(V))';
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

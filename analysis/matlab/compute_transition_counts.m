function Nij = compute_transition_counts(s_t, tau)
% Compute transition counts at specified lag time.
%
% Nij = compute_transition_counts(s_t, tau)
%
% ARGUMENTS
%   s_t - discretized trajectory on 1...S where S = max(s_t)
%   tau - lag time
% RETURNS
%   Nij - S x S matrix of transition counts where Nij(i,j) is number of observed counts from i to j for lag time tau

S = max(s_t); % number of states
T = length(s_t);

% Allocate storage for transition counts.
Nij = zeros(S,S);

% Accumulate counts.
for i = 1:S
  for j = 1:S
    Nij(i,j) = sum( s_t(1:T-tau)==i & s_t(1+tau:T)==j );
  end
end

return




function [potential, gradient] = muller(q)
% Muller potential
%
% [potential, gradient] = muller_brown(q)
%
% Compute the potential energy and gradient vector for the Müller potential.
%
% ARGUMENTS
%   q - 2-vector specifying the coordinates at which the potential is to be evaluated
% 
% RETURN VALUES
%   potential - the potential energy at q
%   gradient - 2-vector containing the gradient of the potential
%
%
% REFERENCE
%   Huo S and Straub JE. The MaxFlux algorithm for calculating variationally optimized reaction
%   paths for conformational transitions in many body systems at finite temperature.
%   J. Chem. Phys. 07(13):5000.  (Eq. 15)

coeff(1:4){'A'} = [-200 -100 -170 15];
coeff(1:4){'a'} = [-1 -1 -6.5 .7];
coeff(1:4){'b'} = [0 0 11 0.6];
coeff(1:4){'c'} = [-10 -10 -6.5 0.7];
coeff(1:4){'x0'} = [1 0 -0.5 -1];
coeff(1:4){'y0'} = [0 0.5 1.5 1];

% Compute contributions to potential from all sets of terms.
potential = 0.0;
gradient = zeros(size(q));
for i = 1:4
  % retrieve coefficients

  % compute contributions to potential and gradient
  term = A*exp(a*(q(1)-x0)^2 + b*(q(1)-x0)*(q(2)-y0) + c*(q(2)-y0)^2);  
  potential = potential + term;
  gradient(1) = gradient(1) + term*(2*a*(q(1)-x0) + b*(q(2)-y0));
  gradient(2) = gradient(2) + term*(b*(q(1)-x0) + 2*c*(q(2)-y0));    
end

return;

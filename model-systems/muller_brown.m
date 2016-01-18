function [potential, gradient] = muller_brown(q)
% Muller-Brown potential
%
% [potential, gradient] = muller_brown(q)
%
% Compute the potential energy and gradient vector for the Muller-Brown potential.
%
% ARGUMENTS
%   q - (2xN)-vector specifying the coordinates at which the potential is to be evaluated
% 
% RETURN VALUES
%   potential - N-vector the potential energy at qs
%   gradient - (2xN)-vector containing the gradient of the potential
%
% REFERENCE
%   Müller K and Brown LD. Location of saddle points and minimum energy paths by a constrained
%   simplex optimization procedure. Theoret. Chim. Acta 53:75, 1979.  (Footnote 7)
%
% Footnote 7 is reproduced below:
%
%   The potential consists in a sum of four terms of the form A.exp [a(x - x0) 2 + b(x - x0) 
% (y - y0) + e(y - y0)2], in which the constants take the following values: (A) = (-200/-100/ 
% -170/15), (a) = (-1/-1/-6.5/0.7), (b) = (0/0/11/0.6), (c) = (-10/-10/-6.5/0.7), (x0) = 
% (1/0/-0.5/-1), (y0) = (0/0.5/1.5/1). 

% Set coefficients.
coeff = struct();
coeff.A = [-200 -100 -170 15];
coeff.a = [-1 -1 -6.5 .7];
coeff.b = [0 0 11 0.6];
coeff.c = [-10 -10 -6.5 0.7];
coeff.x0 = [1 0 -0.5 -1];
coeff.y0 = [0 0.5 1.5 1];

% Determine number of elements to evaluate.
N = size(q,2);

% Compute contributions to potential from all sets of terms.
potential = zeros(1,N);
gradient = zeros(2,N);
for i = 1:4
  % Retrieve coefficients for term i.
  A = coeff.A(i);
  a = coeff.a(i);
  b = coeff.b(i);
  c = coeff.c(i);
  x0 = coeff.x0(i);
  y0 = coeff.y0(i);
  
  % Compute contributions to potential and gradient.
  term = A*exp(a*(q(1,:)-x0).^2 + b*(q(1,:)-x0).*(q(2,:)-y0) + c*(q(2,:)-y0).^2);  
  potential = potential + term;
  gradient(1,:) = gradient(1,:) + term.*(2*a*(q(1,:)-x0) + b*(q(2,:)-y0));
  gradient(2,:) = gradient(2,:) + term.*(b*(q(1,:)-x0) + 2*c*(q(2,:)-y0));    
end

return;

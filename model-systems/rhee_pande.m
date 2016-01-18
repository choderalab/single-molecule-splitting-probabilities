function [potential, gradient] = rhee_pande(q)
% Rhee and Pande potential.
%
% [potential, gradient] = rhee_pande(q)
%
% Compute the potential energy and gradient vector for the Rhee and Pande potential described 
% in Eq. 41 of Ref. [1].
%
% ARGUMENTS
%   q - 2xN vector of [x;y] vectors
% 
% RETURN VALUES
%   potential - the potential energy at q
%   dUdx - x-gradient
%   dUdy - y-gradient
%
% REFERENCE
%
% [1] Rhee YM and Pande VS. One-dimensional reaction coordinae and the corresponding potential of mean
% force from commitment probability distribution. JPC B 109:6780, 2005.

[M,N] = size(q);

x = q(1,:);
y = q(2,:);

potential = (1 - 0.5*tanh(y-x)).*(x+y-5).^2 + 0.2*(((y-x).^2 - 9).^2 + 3*(y-x)) + 15*exp(-(x-2.5).^2 - (y-2.5).^2) - 20*exp(-(x-4).^2 - (y-4).^2);

gradient = zeros(2,N);
gradient(1,:) = (1 - 0.5*tanh(y-x)).*2.*(x+y-5) + (-0.5*sech(y-x).^2)*(-1).*(x+y-5).^2 + 0.2*(2*((y-x).^2 - 9).*2.*(y-x).*(-1) + 3*(-1)) + 15*exp(-(x-2.5).^2 - (y-2.5).^2).*(-2).*(x-2.5) - 20*exp(-(x-4).^2 - (y-4).^2)*(-2).*(x-4);
gradient(2,:) = (1 - 0.5*tanh(y-x)).*2.*(x+y-5) + (-0.5*sech(y-x).^2)*(+1).*(x+y-5).^2 + 0.2*(2*((y-x).^2 - 9).*2.*(y-x).*(+1) + 3*(+1)) + 15*exp(-(x-2.5).^2 - (y-2.5).^2).*(-2).*(y-2.5) - 20*exp(-(x-4).^2 - (y-4).^2)*(-2).*(y-4);

return;

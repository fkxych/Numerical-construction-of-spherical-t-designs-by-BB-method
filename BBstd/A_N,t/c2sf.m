function S = c2sf(X)
% S = c2sf(X)
% Convert set of m cartesion points X on unit sphere S^2 in R^3
% points are stored in columns of X (3 by m) 
% into spherical coordinates theta in [0, pi), phi in [0, 2*pi], r >= 0
% S is 3 by m matrix of elements [theta; phi; r]

% For points on unit sphere r == ones(1,m)
r = sqrt(sum(X.*X));

% Avoid division by 0 if points may be zero
%r(r==0) = eps;

% Polar angle theta
theta = acos(X(3,:)./r);

% Azimuthal angle phi on [0, 2*pi)
phi = atan2(X(2,:), X(1,:)) + (2*pi)*(X(2,:)<0);;

% Matrix of spherical coordinates
S = [theta; phi; r];

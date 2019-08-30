function s = c2sn(X)
% s = c2sn(X)
% Convert set of m normalized symmetric cartesion points X on unit sphere S^2 in R^3
% Points are stored as columns of X (3 by m) where X(:,1) is the north pole
% and X(:,2) is on the prime meridian,
% into vector of spherical coordinates theta in [0, pi), phi in [0, 2*pi]
% with theta(1) = 0 (phi(1) irrelavant) and phi(2) = 0 omitted.
% s is vector of 2*m-3 elements [theta(2:m) phi(3:m)]

m  = size(X,2);

% For points on unit sphere r == ones(1,m)
r = sqrt(sum(X.*X));

% Avoid division by 0 if points may be zero
%r(r==0) = eps;

% Polar angle theta
theta = acos(X(3,:)./r);

% Azimuthal angle phi on [0, 2*pi)
phi = atan2(X(2,:),X(1,:)) + (2*pi)*(X(2,:)<0);;

% Omit theta(1) and phi(1:2) assuming first point is noth pole
% and second point is on primer meridian
s = [theta(2:end) phi(3:end)];

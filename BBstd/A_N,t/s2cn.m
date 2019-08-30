function X = s2cn(s);
% X = s2cn(s);
% Convert the spherical representation s of a fundamental system
% on the unit sphere S^2 in R^3 into cartesian coordinates X
% s is a row vector containing the angles theta and phi
% Because of the rotational invariance of the fundamental system
% first point is always (x,y,z) = (0, 0, 1) <=> theta = 0 (and phi = 0)
% second point always has phi = 0.
% Thus there are m-1 variables theta and m-2 variables phi
% X is a 3 by m matrix containing the m cartesian points on the sphere

ls = length(s);
m = (ls+3) / 2;

s = s(:)';

theta = [0, s(1:m-1)];
phi = [0, 0, s(m:ls)];

X = zeros(3,m);
X(1,:) = sin(theta) .* cos(phi);
X(2,:) = sin(theta) .* sin(phi);
X(3,:) = cos(theta);


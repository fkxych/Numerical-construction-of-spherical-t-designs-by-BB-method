function [A, DAdt, DAdp]=inmds(X, n, ipr)
% [A, DAdt, DAdp] = inmds(X, n, ipr)
%   Return the interpolation matrix A based on the normalized
%   real spherical harmonics Y_{L,K} of degree L at the points x,
%   for -L <= K <= L and 0 <= L <= n,n is t
%   where the default for n has size(X,2) == (n+1)^2
%   Optionally, return DAdt and DAdp.
%   DAdt is the partial derivative of A with respect to theta,
%   DAdp is the partial derivative of A with respect to phi,
%   where theta, phi are the spherical polar coordinates of x.
%
%   See http://mathworld.wolfram.com/SphericalHarmonic.html
%   for definition.
%   See also http://mathworld.wolfram.com/LegendrePolynomial.html
%

if nargin < 3
    ipr = 0;
end;

t0 = cputime;
m = size(X,2);
if nargin < 2
    n = floor(sqrt(m))-1;
end;

CalcDAdt = nargout > 1;
CalcDAdp = nargout > 2;

dn = (n+1)^2;
A    = zeros(dn, m);
if CalcDAdt
    DAdt = zeros(dn, m);
end;
if CalcDAdp
    DAdp = zeros(dn, m);
end;

S       = c2sf(X);
z       = X(3,:)./S(3,:);
phi     = S(2,:);

% Do for each degree
for L = 0:n
    slice = [L^2+1:L^2+2*L+1];
    if CalcDAdp
        [A(slice,:), DAdt(slice,:), DAdp(slice,:)] = spharmrds(L,z,phi);
    elseif CalcDAdt
        [A(slice,:), DAdt(slice,:)] = spharmrds(L,z,phi);
    else
        A(slice,:) = spharmrds(L,z,phi);
    end
end;

if ipr > 0
    fprintf('INMDS: Time taken for degree n = %d and m = %d points: %.2f secs\n', n, m, cputime-t0);
end;
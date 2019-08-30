function [G, DG, DDG] = gramxddL(X, L, ipr)
% [G, DG] = gramxddL(X, L, ipr)
% Calculate the Gram matrix G for degree L and
% its deriviative DG w.r.t. Z = X'*X
% for the fundamental system X on the Sphere S(2) in R^3
% X is an 3 by m matrix contaning the fundamental system
% G is a m by m symmetric positive semidefinite matrix
% G is nonsingular and hence positive definite if X is a fundametnal system
% DG is a m by m symmetric matirx of derivatives of G w.r.t Z = X'*X
% Derivative of G w.r.t. X(i,j) is  ej*ej'*DG*Di + Di*DG*ej*ej'
% where Di = diag(X(i,:)) and ej is j th unit vector in R^m.
if nargin < 3
    ipr = 0;
end;

CalcDeriv = nargout - 1;

t0 = cputime;

% For the sphere in R^3 the number of interpolation points N = (L+1)^2
% L is the degree of the interpolating polynomials on the sphere
N = size(X,2);
if nargin < 2
    L = sqrt(N) - 1;
end;

% Area of unit sphere S^2 in R^3
% Not used in this version, as in phiL
S2area = 4*pi;

% Evaluate using symmetry of G, DG and Gegenbauer polynomials
% Should initliaze diagonal elements of G and DG
% and then only use z = X(:,i)'*X(:,i+1:m) to get G(i,i+1:m) etc
% Explicity using the symmetry roughly halves the number of flops
% to caclulate G, but introduces a slow for loop
% This is faster for n >= 20 or so
Gdiag = (L+1)^2 - 1;                 % (L+1)^2 - 1
DGdiag = L*(L+2)*Gdiag / 4;          %  L * (L+1)^2 * (L+2) / (4); 
DDGdiag = (L-1)*(L+3)*DGdiag / 6;    % (L-1) * L * (L+1)^2 * (L+2) * (L+3) / (24);
G = diag(Gdiag * ones(N,1));
if CalcDeriv > 0
    DG = diag(DGdiag * ones(N,1));
    if CalcDeriv > 1
        DDG = diag(DDGdiag * ones(N,1));
    end;
end;

% Loop over strict upper traingle as diagonal specified
% Strict lower triangle given by symmetry
for i = 1:N-1
    
    I = [i+1:N];
    xi = X(:,i);
    XI = X(:,I);
    z = xi'*XI;
    
    % Make sure elements of A are in [-1, 1]
    z = max(z, -1);
    z = min(z, 1);
    
    if CalcDeriv == 2
        [p, Dp, DDp] = phiL(L, z);
        G(i,I) = p;
        G(I,i) = p';
        DG(i,I) = Dp;
        DG(I,i) = Dp';
        DDG(i,I) = DDp;
        DDG(I,i) = DDp';
    elseif CalcDeriv == 1
        [p, Dp] = phiL(L, z);
        G(i,I) = p;
        G(I,i) = p';
        DG(i,I) = Dp;
        DG(I,i) = Dp';
    else
        p = phiL(L, z);
        G(i,I) = p;
        G(I,i) = p';
    end;
    
end;

tc = cputime  - t0;

if ipr > 0
    fprintf('Calculating Phi_L matrix and %d derivatives for m = %d', CalcDeriv, N);
    fprintf(': Time = %.2f secs\n', tc);
end;

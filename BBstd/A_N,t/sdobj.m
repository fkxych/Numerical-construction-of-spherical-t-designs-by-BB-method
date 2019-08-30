function [f, g, H] = sdobj(arg1, L, ipr)
% [f, g, H] = sdobj(arg1, L, ipr)
% Calculate the spherical L-design objective
% f =  A_N,L = (4*pi/N^2) sum(sum(Phi))
% where the gram matrix Phi depends on the N points in arg1
% arg1 can be a 3 by N array of cartesion points
% or a 2N-3 vector (either row or column) of the normalized
% spherical parametrization of the points (see c2sn and s2cn)
% If two output arguments are given also calculate gradient w.r.t.
% spherical parametrization s
% If ipr > 0 (default ipr = 0) print summary

t0 = cputime;

if min(size(arg1)) == 1
   % arg1 a vector of normalized sperical paramaetrization
   s = arg1;
   X = s2cn(s);
else
   % arg1 should be a 3 by N array of cartesian points in S^2
   X = arg1;
   s = c2sn(X);
end;
% Number of points
N = size(X,2);

% If degree L is not specified, assume N = (L+1)^2
if nargin < 2, L = sqrt(N)-1; end;
if isempty(L), L = sqrt(N)-1; end;

% Default to no output
if nargin < 3, ipr = 0; end;

% Decide whether to calculate gradient
CalcDeriv = nargout > 1;

if CalcDeriv
    [Phi, DPhi] = gramxddL(X, L);
else
     Phi = gramxddL(X, L);
end;

% Shperical design objective 
%wavg = 4*pi/N;
% Moved contant |S^2| = 4*pi into calculation of Phi
const = 1/N^2;
%f = const*sum(sum(Phi));
%Phi(i,i) = (L+1)^2 - 1;
f = const*sum(sum(Phi));
clear Phi

if CalcDeriv
        
    s = s(:)';
    ns = length(s);
    theta = [0, s(1:N-1)];
    phi = [0, 0, s(N:end)];
    ct = cos(theta);
    st = sin(theta);
    cp = cos(phi);
    sp = sin(phi);
    
    g = zeros(size(s));
    
    e = ones(N,1);
    XDG = X*DPhi;
    XDG1 = (e*X(1,:)) .* DPhi;
    XDG2 = (e*X(2,:)) .* DPhi;
    XDG3 = (e*X(3,:)) .* DPhi;

    for j = 1:N
                
        w1 = XDG1(:,j)'; w1(j) = w1(j) + XDG(1,j);
        w2 = XDG2(:,j)'; w2(j) = w2(j) + XDG(2,j);
        w3 = XDG3(:,j)'; w3(j) = w3(j) + XDG(3,j);
        
        % Derivatives w.r.t spherical parametrization
        Dtheta = w1.*(ct.*cp) + w2.*(ct.*sp) - w3.*st;
        Dphi = -w1.*(st.*sp) + w2.*(st.*cp);
        
        % Derivative of c = wavg*G*e - e w.r.t. variables s
        g = g + [Dtheta(2:N) Dphi(3:N)];
        
    end;
    
    g = const * g';
    
    % Finite difference approximation of Hessian
    if nargout > 2
        
        step = eps^(1/3);
        ns = max(size(s));
        H = zeros(ns,ns);
        for j = 1:ns
            %j
            sj = s;
            sj(j) = sj(j)+step;
            [fj, gjp] = sdobj(sj, L);
            sj(j) = sj(j)-step;
            [fj, gjm] = sdobj(sj, L);
            H(:,j) = (gjp - gjm)/(2*step);
        end;
        H = (H+H')/2;
        
    end;
    
end;

tc = cputime - t0;

if ipr > 0
    fprintf('SDOBJ: Spherical %d-design objective A_{L,N} with %d points\n', L, N);
    fprintf('f = %.6e', f);
    if CalcDeriv
        fprintf(', ||g||_inf = %.4e', norm(g,inf));
    end;
    fprintf(', Time = %.2f secs\n', tc);
end;

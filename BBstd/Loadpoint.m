function [X0]=Loadpoint(t,n)
%Author: Yuchen Xiao
%Time: Nov 18, 2018
% Load md point data
% The form of data is Cartisian framework [X,Y,Z]

tstr=sprintf('%.2f',t/1e2);
tstr=tstr(3:end);
nstr=sprintf('%.4f',n/1e4);
nstr=nstr(3:end);
X0=['load md' tstr '.' nstr];
eval(X0);
eval(['X0=md' tstr '(:,1:3);']);
X0=X0';
end
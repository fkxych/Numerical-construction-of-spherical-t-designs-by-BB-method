function [X0]=Loadpoint2(t,n)
%Author: Yuchen Xiao
%Time: Nov 18, 2018
% Load md point data
% The form of data is Cartisian framework [X,Y,Z]

prompt = 'Please input the file name ';
str = input(prompt,'s');
X0=load(str);
X0=X0(:,1:3);
X0=X0';
end
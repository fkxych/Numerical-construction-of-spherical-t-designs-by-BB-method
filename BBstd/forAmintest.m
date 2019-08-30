%Computing A_{N,t} with t up to 191 with N = t^2 + 2t + O(1) points;
%Initial points sets-->Extremal (maximum determiant) points. 
%Sources of points: https://web.maths.unsw.edu.au/~rsw/
%Author: Yuchen Xiao
%Time: Nov 18, 2018
%%
  clc,clear
  close all
  F=[]; T=[]; MY=[]; MY1=[];
  currentFolder = pwd;
  addpath(genpath(currentFolder));
  profile on
for t= 50:50
           N = (t+1)^2;    
    if t>80
          X0 = Loadpoint2(t,N);          
    else
       X0 = Loadpoint(t,N);
    end
%   tic
%     [f,XX] = Amintest(X0,t);
    [f,XX,minY,minY1] = Amintest(X0,t);
%   toc
    F(end+1,:) = f;
    T(end+1,:) = t;
    MY(end+1,:) = minY;
    MY1(end+1,:) = minY1;
end

  figure(11),plot(T,log(F),'*'),grid on,xlabel('t');title('The behavior of final log {\it A_{N,t}} ','fontSize',12)
  figure(12),plot(T,F,'*'),grid minor,xlabel('t');title('The behavior of final {\it A_{N,t}} ','fontSize',12)
  figure(13),plot(T,MY,'*'),grid on,xlabel('t');title('The behavior of minimal singular vaule of Y_{t}','fontSize',12)
  figure(14),plot(T,MY1,'*'),grid on,xlabel('t');title('The behavior of minimal singular vaule of Y_{t+1}','fontSize',12)
  
  profile viewer
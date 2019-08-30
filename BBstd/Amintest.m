function [f,XX,minY,minY1,ht,at,rs] = Amintest(X0,t)
%Minimize spherical design objective A_{N,t} from the following paper:
%I. H. Sloan, R. S. Womersley, A variational characterisation of spherical designs, 
%Journal of Approximation Theory 159 (2) (2009) 308-318.
%Original Code Author: Rob Womersley and Ian Hugh Sloan
%Updater and Modifier: Yuchen Xiao and Congpei An
%Time: Nov 18, 2018
%%
%Preprocess
   format long
    % Loaded initial point set
    s0 = c2sn(X0);
    fprintf('No perturbation of intial points\n');

%%
%Select optimization method
%Minimize objective A_{L,N}
%By using minFunc
    options = [];
    options.display = 'iter';

    %qn
%     options.Method = 'qnewton';
%     options.qnUpdate = 0;%BFGS
    %bb
    options.Method = 'bb';
    options.bbType = 1;

    options.MaxIter = 10000;
    options.optTol = 1e-16;
    options.progTol = 1e-16;
%      options.LS_type = 0;
%     options.LS_init = 2;
%     options.Fref = 20;
%     options.useMex = 1;
    ipr=0;%in sdobj decide whether to show the result(at line 111 sdobj.m)
    varargin=[t,ipr];
    [sss,fs,flag,output] = minFunc(@sdobj,s0,options,varargin);
    fprintf('FMINUNC: flag = %d\n', flag);
    [f, g] = sdobj(sss, t, 3);
    XX = s2cn(sss);

%%
%Output
if nargout>2
fprintf('\n ---- terminate point set -----\n');
%Fundamental systems
  Y=inmds(XX,t);
  svdY=svd(Y);
  minY=min(svd(Y));
  maxY=max(svd(Y));

  Y1=inmds(XX,t+1);
  svdY1=svd(Y1);
  minY1=min(svd(Y1));
  maxY1=max(svd(Y1));
else
    
end

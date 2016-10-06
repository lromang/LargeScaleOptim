function[x, iter] = LBFGS ( fgname, maxiter, tol, m, prtlevel )
%
% ...This is an implementation of the L-BFGS method for solving large-scale
% unconstrained optimization problems. The routine requires the objective
% function and its gradient (f, g).
% This version makes use of automatic differentiation, provided by AMPL, 
% to compute f and g.  The interface AMPL-Matlab is through spamfunc.
%
% routine loop_LBFGS provided by javier urena carrion
%
% j-l morales
% itam 2015
%--------------------------------------------------------------------------

[ x0, ~, ~, ~, ~, ~ ] = spamfunc( fgname );   
n  = length(x0);
sn = sqrt(n);

output = strcat(fgname,'.outL');
fout = fopen(output,'w');

fprintf(fout, ' Name of the problem  % 10s \n', fgname);
fprintf(fout, ' Number of variables  % 5i  \n', n);

x    = x0;  
fev  = 1;
ind  = 1;
iter = 0;
%
% More-Thuente line search translated to Matlab  by DP O'Leary ------------
%
c1     = 1.0e-4; 
c2     = 0.9;
xtol   = 0.1e0;
stpmin = eps;
stpmax = 1.0e0;
maxfev = 20;
%
%--------------------------------------------------------------------------
%
tic
S = zeros(n,m);
Y = zeros(n,m);
[f, ~] = spamfunc(x, 0);
[g, ~] = spamfunc(x, 1);
norm_g = norm(g,inf);
normg0 = norm_g;

fprintf(fout,' %3i  %21.15e   %8.2e   %5.3f    %8s  \n', ...
                            iter, f,  norm_g/(normg0), [], []);

fprintf(fout,' Iter            f              ||g||    alpha        s*y \n');
while  norm_g > tol*(1.0 + normg0)  &&  iter < maxiter   
    % Add this restriction to C code
    if iter >= 1
        p = - loop_LBFGS(S, Y, g, gamma_k, ind);
    else 
        p = - g;
    end
    g0 = g;

    alpha = 1;
    [ x, f, g, alpha, info, nfev ] = cvsrch( fgname, n, x, f, g, ...
               p, alpha, c1, c2, xtol, stpmin, stpmax, maxfev);

    fev = fev + nfev;
    % 
    % ... store pair (s,y)
    %
    s = alpha * p;
    y = g - g0;   
    sTy    = s'*y;
    
    
    if sTy <= eps
        fprintf(fout, ' Warning, small sTy < eps  % 8.2e \n', sTy);
    else        
        a = mod(iter,m) + 1;
        ind(a)= iter + 1;
        S(:,a) = s;
        Y(:,a) = y;
 
        gamma_k  = (sTy / (y' * y));
        norm_g = norm(g,inf);
        iter   = iter + 1;
    end
   
    %if prtlevel > 0 & mod(iter,0) == 0
       fprintf(fout, ' %3i  %21.15e   %8.2e   %5.3f    %8s  \n', iter, f, norm_g/(1.0+normg0), alpha, sTy );
    %end
end
fprintf(fout, ' %3i  %21.15e   %8.2e   %5.3f    %8s  \n', ...
                           iter, f,  norm_g/(1.0+normg0), alpha, s'*y );
fprintf(fout, ' Number of function evaluations %4i \n', fev );

toc

% tic
% semilogy(normg);
% toc
fclose(fout);

function [r]= loop_LBFGS( s, y, g, gamma, indice);
%--------------------------------------------------------------------------
%
%
% ... Doble loop para encontrar producto H_k*grad_k en el metodo BFGS con
% memoria limitada
%
% ENTRADA:
%      s   matriz con vectores s
%      y   matriz con vectores y
%      g   gradiente en el punto actual
%     H0   aproximacion incial a H_k
%  indice  vector de indices que indique el orden de 's',  'y' 

% SALIDA:
%      r   Aproximacion a H_k * grad_k
%    
%      Javier Urena Carrion    
%      CU: 125319           ITAM 2015
%
%--------------------------------------------------------------------------    
%indice
m      = length(indice);
[~,b]  = min(indice);
[~,bm] = max(indice);
a = b;
q = g;

ro = zeros(m,1);
for i = 1:m
   ro(a) = 1 / (y(:,a)' * s(:,a));
   a = mod(a,m) + 1;
end
a = bm;
alfa=zeros(m,1);
%
% Two loops start here
%
for i=1:m
    alfa(a)= ro(a) * s(:,a)' * q;
    q = q - alfa(a) * y(:,a);
    a = mod(a - 2,m) + 1;
end
r = gamma * q;
a = b;
for i=1:m
    beta = ro(a) * y(:,a)' * r;
    r = r + (alfa(a) - beta) * s(:,a);
    a = mod(a,m) + 1;
end
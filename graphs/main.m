clear
clc

N = 4;

xfix = [ 1:N ]'
% 
% cond = 16;
% 
% D = diag(linspace(1,cond,N))
% 
% W = [ .5 -.5 .5 .5 ]'
% 
% Q = eye(N) - 2 * (W * W')
% 
% A = Q*D*Q'
% 
% b = A*xfix

A = [
    4 1 0 0
    1 4 1 0
    0 1 4 1 
    0 0 1 4
]

b = A*xfix

x2 = [1.0001 1.9998 3.0005 3.9999]'

norm(x2 - xfix, 1)
norm(x2 - xfix, 1)/norm(xfix, 1)
norm(A*x2 - b, 1)
clear
clc

N = 4;

xfix = [ 1:N ]'

cond = 16;

D = diag(linspace(1,cond,N))

W = [ .5 -.5 -.5 .5 ]'

Q = eye(N) - 2 * (W * W')

A = Q*D*Q'

b = A*xfix

[L,U, P] = lu(A)

sum(sum(P == eye(N))) / N^2

L*U

% L(4,3) = 1
% 
% U(4,4) = (A(4,4) - L(4,1:4-1)*U(1:4-1,4)) / U(4,4)

n = N;
% Ly = b
y = zeros(n,1);
y(1) = b(1);
for i = 2:n
    y(i) = b(i) - L(i,1:i-1)*y(1:i-1);
end
y

% Ux = y
x = zeros(n,1);
x(n) = y(n)/U(n,n);
for i = n-1:-1:1
    x(i) = (y(i) - U(i,i+1:n)*x(i+1:n))/U(i,i);
end
x

x2 = [-0.3498, 1.1773, 2.1753, 4.7272]'

norm(x2 - xfix, 1)
norm(xfix, 1)
norm(A*x2 - b, 1)
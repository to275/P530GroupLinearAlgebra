function [x,er,V,H] = mygmres(l,b,x0,n,M,A)
% [x,er,V,H] = mygmres(l,b,x0,n,M,A)
% function that runs the GMRES algorithm for the problem Ax = b
% Inputs: l: number of iterations, b: soultion vector, x0: initial guess,
% n is the dimension of
% the problem, A is an n x n matrix, M is an n x n matrix used to define
% the inner product: <u,v> = u'Mv
% output: x is the approximate solution to the problem Ax=b, er: error, V:
% set of basis vectors, H: Hessenberg matrix

% define constants from the Fan paper.
r0 = b-A*x0;
beta = norm(r0);
v1 = r0/beta; % first basis vector
% initialize the matrices (don't do this, it causes issues when we don't go
% through the entire thing
% V = zeros([length(x0),n]); % V subspace
V(:,1) = v1;
% H = zeros([length(x0),n]);
% W = zeros([length(x0),n]);
for j = 1:l
    W(:,j) = A*V(:,j);
    for i = 1:j
        H(i,j) = W(:,j)'*M*V(:,i);
        W(:,j) = W(:,j)-H(i,j)*V(:,i);
    end % i = 1:j
    H(j+1,j) = norm(W(:,j));
    if H(j+1,j) <= 1e-8 % close enough to zero
        disp('break')
        break
    end
    V(:,j+1) = W(:,j)/H(j+1,j);
end % j = 1:n

e1 = eye(size(H));
e1 = e1(:,1);
y = H\(beta*e1); % <-- this is running into a conditioning error
x = x0 + V(:,1:length(y))*y;
er = norm(A*x-b)/n;
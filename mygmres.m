function [x,V,H] = mygmres(I,b,x0,n,M,A)
% function [V,H] = mygmres(I,b,x0,n,M,A)
% function that runs the GMRES algorithm for the problem Ax = b
% Inputs: I: number of iterations, x0: initial guess, n is the dimension of
% the problem, A is an n x n matrix, M is an n x n matrix used to define
% the inner product: <u,v> = u'Mv
% output: x is the approximate solution to the problem Ax=b

% define constants from the Fan paper.
r0 = b-A*x0;
beta = norm(r0);
v1 = r0/beta; % first basis vector
% initialize the matrices
V = zeros([length(x0),n]); % V subspace
V(:,1) = v1;
H = zeros([length(x0),n]);
W = zeros([length(x0),n]);
for j = 1:n
    W(:,j) = A*V(:,j);
    for i = 1:j
        H(i,j) = W(:,j)'*M*V(:,i);
        W(:,j) = W(:,j)-H(i,j)*V(:,j);
    end % i = 1:j
    H(j+1,j) = norm(W(:,j));
    if H(j+1,j) == 0
        m=j;
        break
    end
    V(:,j+1) = W(:,j)/H(j+1,j);
end % j = 1:n

e1 = eye(size(H));
e1 = e1(:,1);
y = H\(beta*e1);
x = x0 + V*y;
function [A,b] = CalcAandb(n,gamma)
%This script will calculate the matrix A and b for the GMRES problem
% Input: n - the number of basis functions, gamma the value of the
% potential: V(x) = gamma.
% Output: A, a tridiagonal matrix with elements as defined in part b
% b - the source vector.

% NEED TO UPDATE TO CREATE A SPRASE REPRESENTATION

% regular representation
dx = n+1;
A1 = zeros(n);
A2 = zeros(n);
b = dx*ones([n,1]);

for i = 2:n % assign off diagonal elements
    A1(i,i-1) = -1/dx;
    A2(i,i-1) = -gamma/2;
end
A1 = A1+A1';
A2 = A2-A2';
A1 = A1+eye(n)*2/dx;

A = A1+A2;

end


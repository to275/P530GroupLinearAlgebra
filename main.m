%% Linear Algebra Group Assignment
% Authors: Tyce Olaveson and Levi Moats 
% Class: Physics 530 
% Professor: Dr Trastrum 
% Date: 02/26/2024 
% 
% Consider the one dimenstional boundary value problem that arises in fluid
% dynamics
%
% $$-u''(x) + V(x)u'(x) = f(x), x\in[0,1]$$
%
% where
%
% $$u(0) = u(1) = 0$$
%
% where we will take $V(x),f(x)$ to be constants: $V(x) = \gamma$ and $f(x)
% = 1$
%% Problem Formulation
% a) Using a test function $\phi$, the above boundary value problem can be
% rewritten in the weak form as the following:
%
% $$A_1[\phi,u] + A_2[\phi,u] = F[\phi]$$
%
% where 
%
% $$A_1[\phi,u] = \int_{0}^{1} \phi' \cdot u' dx$$
%
% $$A_2[\phi,u] = \int_{0}^{1} v(x) \cdot u' \cdot \phi dx$$
%
% $$F[\phi] = \int_{0}^{1} f \cdot \phi dx$$
%
% b) The "hat" functions from class, $\phi_i$, can be used to aproximate
% $u$ using linear combinations: $u(x) = \sum_{i} u_i \phi_i(x)$. This
% allows the weak form of the boundary value problem to be written in the
% form:
%
% $$Ax = b$$
%
% where
%
% $$A = A_1 + A_2$$
% 
% $$b = F_i$$
%
% $$x = u(x)$$
%
% From above
%
% It can be shown that $A_1$ is a symetric tridiagonal matrix with
% $\frac{2}{\Delta x}$ on the diagonal and $\frac{-1}{\Delta x}$ on the off
% diagonals. $A_2$ is a skew-symmetric tridiagonal matrix with 0 on the
% diagonal, $\frac{\gamma}{2}$ on the upper off diagonal, and 
% $\frac{-\gamma}{2}$ on the lower off diagonal.
%
% It can be shown that $F_i$ is a vector is an n by 1 vector where all the
% values are $\Delta x$.
%
% c) A function CalcAandb.m has been written and is included in this
% repository. It accepts $n$ and $\gamma$ as inputs and returns a sparse A
% matrix and 
%% Implement the GMRES algorithm
%
%% Solving the FEM problem
%
%% Preconditioning the GMRES
%
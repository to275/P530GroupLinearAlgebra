%% Linear Algebra Group Assignment
% Authors: Tyce Olaveson and Levi Moats 
% Class: Physics 530 
% Professor: Dr Transtrum 
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
% Examples from the text:
% A = [1 4 7; 2 9 7; 5 8 3];
% b = [1;8;2]; x = [-2.18; 1.84; -0.6]
% b = [2;3;9]; x = [-2.1; -0.22; 0.11]
% b = [5;-3;8]; x = [4.8; -2.6; 1.5]
A = [1 4 7; 2 9 7; 5 8 3];

b = [1;8;2]; x = [-2.18; 1.84; -0.6]; x0 = zeros([3,1]); M = eye(3);
[xfit,er,V,H] = mygmres(3,b,x0,3,M,A);
disp(round(xfit,2))
disp(norm(er))

b = [2;3;9]; x = [-2.1; -0.22; 0.11]; x0 = zeros([3,1]); M = eye(3);
[xfit,er,V,H] = mygmres(3,b,x0,3,M,A);
disp(round(xfit,2))

b = [5;-3;8]; x = [4.8; -2.6; 1.5]; x0 = zeros([3,1]); M = eye(3);
[xfit,er,V,H] = mygmres(3,b,x0,3,M,A);
disp(round(xfit,2))
disp(' ')
%% Solving the FEM problem
% Use your GMRES function to solve the finite-element formulation of the
% variational problem for the cases $V(x) = 1$ and $V(x) = n+1$ using $M$
% as the identity matrix. For each case, run with $n=16,32,64,128$ and
% $l=2,4,8,16,32,64,...,$ increasing $l$ until the error (i.e., the norm of
% the residual divided by $n$) is below $10^{-6}$. Plot your most accurate
% solution (as a function of x) as well as teh error versus functions of
% $n$ and $l$.

% define arrays
ns = [16,32,64,128]; % number of basis functions to use
ls = 2.^(1:7); % number of iterations

Vs = {@(n) 1, @(n) n+1}; % set of different (constant) functions for

errors = zeros([length(ls),length(ns),length(Vs)]); % empty matrix to store the errors
solutions = cell(size(errors)); % empty cell array to store the solutions after they have been calculated
xplot = cell([length(ns),1]);
for i = 1:length(ns)
    xplot{i} = 0:1/(ns(i)+1):1;
end

% loop through each case
for i = 1:length(Vs)
    V = Vs{i};
    for j = 1:length(ns)
        n = ns(j);
        for k = 1:length(ls)
            l = ls(k);
            % create the input matrices
            [A,b] = CalcAandb(n,V(n));
            [x,errors(k,j,i),~,~] = mygmres(l,b,zeros([n,1]),n,eye(n),A);
            solutions{k,j,i} = [0;x;0];
        end  % k = 1:length(ls)
    end % j = 1:length(ns)
end % i = 1:length(Vs)

%%%
% Produce the plot for $V(x) = 1$
sq=squeeze(errors(:,:,1));

figure
[X,Y] = meshgrid(ns,ls);
surf(X,Y,sq)
ax=gca;
ax.ZScale = 'log';
title("V(x) = 1","FontSize",16)
ylabel("l","FontSize",14)
xlabel("n","FontSize",14)
zlabel("Error","FontSize",14)
view([143.7,25.8])
%zlim([1e-6 max(sq,[],"all")+max(sq,[],"all")*0.05])

hold on
Hplane = surf(X,Y,1e-6*ones(size(X)));
Hplane.EdgeColor = 'none';
Hplane.FaceAlpha = 0.5;
snapnow;

figure()
for i = 1:length(ns)
    plot(xplot{i},solutions{7,i,1},'.-')
    hold on
end
legend(string(ns))
hold off
xlabel('x')
ylabel('u')
%%%
% Produce the plot for $V(x) = n+1$
sq=squeeze(errors(:,:,2));

figure
surf(X,Y,sq)
ax=gca;
ax.ZScale = 'log';
title("V(x) = n+1","FontSize",16)
ylabel("l","FontSize",14)
xlabel("n","FontSize",14)
zlabel("Error","FontSize",14)
view([143.7,25.8])
%zlim([1e-6 max(sq,[],"all")+max(sq,[],"all")*0.05])

hold on
Hplane = surf(X,Y,1e-6*ones(size(X)));
Hplane.EdgeColor = 'none';
Hplane.FaceAlpha = 0.5;
snapnow;

figure()
for i = 1:length(ns)
    plot(xplot{i},solutions{7,i,2},'.-')
    hold on
end
legend(string(ns))
hold off
xlabel('x')
ylabel('u')
%% Preconditioning the GMRES
% a) Simple derivation here
% Let $x$ be a solution to $Ax=b$. We now left multiple both sides by
% $A_1^{-1}$ to get $A_1^{-1}Ax=A_1^{-1}b$ which is equivalent to
% $\tilde{A}x = \tilde{b}$ by definition of $\tilde{A}$ and $\tilde{b}$.
% Thus $x$ is a solution to the preconditioned problem.
%
% Now suppose that $x$ is a solution to the preconditioned problem : $\tilde{A}x =
% \tilde{b}$. We multiply both sides by $A_1$ and expand to get:
% $A_1 A_1^{-1} A x = A_1 A_1^{-1} b \rightarrow Ax = b$.
% Thus the solution to the preconditioned problem is also a solution to the
% original problem.
%  
% 
% b) $\tilde{A}$ and $\tilde{b}$ are easy to calculate be cause they are
% both derived using $A_1$. $A_1$ is a symetric tridiagonal matrix making it very
% easy to find the invers reguardless of size using a variety of
% algorithms. 
%
% c) Repeat problem 3 using the preconditioned matrix. 
%

% define arrays
ns = [16,32,64,128]; % number of basis functions to use
ls = 2.^(1:7); % number of iterations

Vs = {@(n) 1, @(n) n+1}; % set of different (constant) functions for

errors = zeros([length(ls),length(ns),length(Vs)]); % empty matrix to store the errors
solutions = cell(size(errors)); % empty cell array to store the solutions after they have been calculated
xplot = cell([length(ns),1]);
for i = 1:length(ns)
    xplot{i} = 0:1/(ns(i)+1):1;
end

% loop through each case
for i = 1:length(Vs)
    V = Vs{i};
    for j = 1:length(ns)
        n = ns(j);
        for k = 1:length(ls)
            l = ls(k);
            % create the input matrices
            [A,b,M] = CalcAandb(n,V(n));
            % convert to the pre-conditioned problem.
            Atilde = M\A; % This solves the problem M Atilde = A
            btilde = M\b;
            disp([cond(A), cond(Atilde)])
            [x,errors(k,j,i)] = mygmres(l,btilde,zeros([n,1]),n,eye(n),Atilde);
            solutions{k,j,i} = [0;x;0];
        end  % k = 1:length(ls)
    end % j = 1:length(ns)
end % i = 1:length(Vs)

%%%
% Produce the plot for $V(x) = 1$
sq=squeeze(errors(:,:,1));

figure
[X,Y] = meshgrid(ns,ls);
surf(X,Y,sq)
ax=gca;
ax.ZScale = 'log';
title("V(x) = 1","FontSize",16)
ylabel("l","FontSize",14)
xlabel("n","FontSize",14)
zlabel("Error","FontSize",14)
view([143.7,25.8])
%zlim([1e-6 max(sq,[],"all")+max(sq,[],"all")*0.05])

hold on
Hplane = surf(X,Y,1e-6*ones(size(X)));
Hplane.EdgeColor = 'none';
Hplane.FaceAlpha = 0.5;
snapnow;

figure()
for i = 1:length(ns)
    [~,mind] = min(errors(:,i,2));
    plot(xplot{i},solutions{mind,i,1},'.-')
    hold on
end
legend(string(ns))
hold off
xlabel('x')
ylabel('u')
title('V(x) = 1')

%%%
% Produce the plot for $V(x) = n+1$
sq=squeeze(errors(:,:,2));

figure
surf(X,Y,sq)
ax=gca;
ax.ZScale = 'log';
title("V(x) = n+1","FontSize",16)
ylabel("l","FontSize",14)
xlabel("n","FontSize",14)
zlabel("Error","FontSize",14)
view([143.7,25.8])
%zlim([1e-6 max(sq,[],"all")+max(sq,[],"all")*0.05])

hold on
Hplane = surf(X,Y,1e-6*ones(size(X)));
Hplane.EdgeColor = 'none';
Hplane.FaceAlpha = 0.5;
snapnow;

figure()
for i = 1:length(ns)
    [~,mind] = min(errors(:,i,2));
    plot(xplot{i},solutions{mind,i,2},'.-')
    hold on
end
legend(string(ns))
hold off
xlabel('x')
ylabel('u')
title('V(x) = n+1')

% c) condition numbers
cond(A)
cond(Atilde)
%%%
% The rate of convergence for the precodnitioned GMRES problem was much
% faster than the previous one. The first $A$ matrix has a condition
% number of 324 while $tilde{A}$ has a conditon number of 220. This means
% that the preconditioned matrix is farther from being singular, causing it
% to converge faster. 

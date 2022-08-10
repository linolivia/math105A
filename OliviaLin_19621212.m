%% Math 105LA Final Project
% Olivia Lin, 19621212

%% Question 1
% Algorithm: Inverse Matrices using Gaussian Elimination
%
% Goal: Row reduce the matrix $[A | I]$ to $[I | A^-1]$, where A is any
% invertible $n\times n$ matrix.
% 

A = [1 -1 2 -1;2 -2 3 -3;1 1 1 0;1 -1 4 3];
tol = 10^-5;

inverse(A)

% Checking that my function produced the correct inverse of matrix A
abs(inv(A)-inverse(A))<tol;
%%
% The absolute error is less than the tolerance, so these two matrices are equal. Our function does, in fact, return the inverse of the matrix.

%% Question 2
% Algorithm: Cholesky Factorization
%
% Goal: Factor a positive definite $n \times n$ matrix $A$ into $LL^T$, where $L$ is lower
% triangular.


% Given Matrix
A = [4 2 2;2 6 2;2 2 5];
b = [0 1 0]';

L = cholesky(A)

% Solving the linear system using Cholesky Factorization
x = cholsolve(L,b)
%%
% * An advantage of Cholesky factorization is that it uses a less number of calculations. When doing Gaussain elimination, one must perform the 
% elementary matrix operations multiple times to reach echelon form.
% However, with Cholesky factorization, one can simply find the lower
% triangular matrix and use that for the forwards and backwards
% substitution, which use no elementary matrix operations.
% 
% * On the other hand, the type of matrix is limited when using Cholesky factorization. In order to use this factorization,
% you must have a symmetric, positive definite matrix, whereas normal
% Gaussain Elimination accepts any matrix.
%% Question 3
% Algortihm: Successive Over-Relaxation (SOR) Method
%
% Goal: solve linear systems that occur in the numerical solution of
% certain partial differential equations


A = [4 -1 0;-1 4 -1;0 -1 4]; Ag = A;
b = [-1 4 -5]';

% D-L-U Decomposition
D = diag(diag(A));
L = D-tril(A);
U = D-triu(A);

% Jacobi Matrix Calculation
Tj = inv(D)*(L+U);
c = inv(D)*b;


tol = 10^-5;
x = zeros(3,1); xg = x;
N = 1000;
rho = (max(abs(eigs(Tj)))); % spectral radius of Tj
w = 2/(1+(sqrt(1-rho^2))); % optimal value of the parameter w* for the SOR method

[x,k] = sor(A,b,x,w,tol,N);

% Jacobi Method
Ab = [4 -1 0 -1;-1 4 -1 4;0 -1 4 -5];

Jacobi_mat(Ab,[0 0 0]',10^-5,1000)

% Gauss-Siedel Method

sor(Ag,b,xg,1,tol,N);

%%
% The SOR method converges after 14 iterations, while the Gauss-Siedel
% method converges after 13 methods.
%
% Like the Gauss-Siedel method, the Jacobi method converges after 13 iterations.
% 
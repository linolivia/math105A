# math105A Final Project
This is a final project that includes codes on iterative methods such as the successive over-relaxation method, Cholesky factorization, and Gaussian elimination
to attain an inverse matrix given a matrix A.

Iterative methods use initial estimations to form sequences that approximate function values, eventually converging to a solution that is below a given tolerance, implying a convergence. In my project, I was assigned to code three different iterative methods to solve different problems.

## Successive Over-Relaxation Method (SOR)
The Successive Over-Relaxation Method, also written as the SOR method, is used to solve a system of linear combinations. This can be in the form of a matrix and vector, written as such:
$$A \vec{x} = \vec{b},$$
where $A$ is a matrix and both $\vec{x}$ and $\vec{b}$ are vectors. In an iterative method of solving this system of equations, we first decompose our matrix into three distinct matrices:
$$A: D-L-U,$$ where $D$ is the main diagonal of matrix $A$, $L$ is the lower triangular matrix, and $U$ is the corresponding upper diagonal matrix. We then set an arbitrary initial vector $\vec{x}^{(0)}$ whose vectors are selected by the mathematician. Then, the new vector, which is, in our case, $\vec{x}^{(1)}$ is obtained in the following equation:
$$\vec{x}^{(1)} = \vec{x}^{(0)} + \omega D^{-1}\vec{r}^{(1)}$$

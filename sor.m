% Algortihm: Successive Over-Relaxation (SOR) Method
% Goal: solve linear systems that occur in the numerical solution of
% certain partial differential equations

function [x,k] = sor(A,b,x,w,tol,N)
k = 1;
XO = x;
[n,~] = size(b);
    while k <= N
        for i = 1:n
            l = 0;
            u = 0;
            for j = 1:i-1
                l = l + A(i,j)*XO(j);
            end
            for j = i+1:n
                u = u + A(i,j)*x(j);
            end
            x(i,:) = (1-w)*x(i) + (w/A(i,i))*(b(i)-l-u);
        end
        x
        if norm(x-XO) < tol
            disp("number of iterations used: "+ k)
            return
        end
        k = k+1;
        XO = x;
    end
 if k>N   
    disp("Maximum number of iterations exceeded")
 end
end
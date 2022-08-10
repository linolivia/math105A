% Solving a linear system with Cholesky Factorization

function x = cholsolve (L,b)
[n,~] = size(L);
y = zeros(n,1);
x = zeros(n,1);
w = 0;
% Forwards Substitution
    for i = 1:n
        for j = 1:i-1
            w = w + L(i,j)*y(j);
        end
        y(i) = (b(i)-w)/L(i,i);
    end
% Backwards Substitution
L = L';
    for i = n:-1:1
        x(i) = y(i);
        for j = i+1:n
            x(i) = x(i) - L(i,j)*x(j);
        end
        x(i) = x(i)/L(i,i);
    end
end
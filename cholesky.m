% Algorithm: Cholesky Factorization
% Goal: Factor a positive definite nxn matrix A into LL', where L is lower
% triangular.

function L = cholesky(A)
    % Step 1
    [n,~] = size(A);
    L = zeros(n);
    L(1,1) = sqrt(A(1,1));
    x = 0;
    y = 0;
    z = 0;
    % Step 2
    for j = 2:n
        L(j,1) = A(j,1)/L(1,1);
    end
    % Step 3
    for i = 2:n-1
        % Step 4
        for k = 1:i-1
            x = x + L(i,k)^2;
        end
        L(i,i) = sqrt(A(i,i)-x);
        % Step 5
        for j = i+1:n
            for k = 1:i-1
                y = y + L(j,k)*L(i,k);
            end
            L(j,i) = (A(j,i)-y)/L(i,i);
        end
    end
    % Step 6
    for k = 1:n-1
        z = z + (L(n,k))^2;
    end
    L(n,n) = sqrt(A(n,n)-z);
end

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
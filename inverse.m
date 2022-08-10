% Algorithm: Inverse Matrix
% Goal: Obtain the inverse of a matrix using Gaussian Elimination


function A = inverse(A)
    [n,~] = size(A); % n is the number of equations or unkonwns
    A = [A eye(n)];
    % Step 1 
    for i = 1:n-1
        % Step 2
        p = find(A(i:n,i),1) + (i-1); %find the first nonzero row index
        if isempty(p)  
            x = 'no unique solution exists';
            return
        end
        % Step 3
        if p~=i
            A([i p],:) = A([p i],:); % swap row i and row p
        end
        % Step 4
        for j = i+1:n
            % Step 5
            m = A(j,i)/A(i,i);
            % Step 6
            A(j,:) = A(j,:) - m*A(i,:);
        end
    end 
    % Step 7
    if A(n,n)==0
        x = 'no unique solution exists';
        return
    end
    % Step 8

     for j = n:-1:1
        for i = j+1:2*n
            A(j,i) = A(j,i)/A(j,j);
        end
        A(j,j) = A(j,j)/A(j,j);
     end
for k = n-1:-1:1
    for j = k:-1:1
        A(j,:) = A(j,:)-A(j,k+1)*A(k+1,:);
    end
end
A = A(:,n+1:2*n);
end
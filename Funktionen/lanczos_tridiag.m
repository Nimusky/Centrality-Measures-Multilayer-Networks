function [T, Q] = lanczos_tridiag(A,b,m)
% Input:    A: symmetric nxn matrix
%           b: starting vector of length n
%           m: maximal number of Lanczos iterations (the method 
%               terminates before if stopping criterion is met)
% Output:   tridiagonal matrix T and the matrix Q with orthonormal columns

b = b/norm(b);
Q(:,1) = b;
alpha = [];
beta = [];
T = [];
for k = 1:m
    if k == 1
        v = A*Q(:,k);
    else
        v = A*Q(:,k)-beta(k-1)*Q(:,k-1);
    end
    alpha(k) = v'*Q(:,k);
    v = v-alpha(k)*Q(:,k);
    beta(k) = norm(v);
    Q(:,k+1) = v/beta(k);
    T(k,k) = alpha(k);
    T(k,k+1) = beta(k);
    T(k+1,k) = beta(k);
    if abs(beta(k))<1e-12
        break
    end
end
T = T(1:end-1, 1:end-1);
Q = Q(:,1:end-1);
end
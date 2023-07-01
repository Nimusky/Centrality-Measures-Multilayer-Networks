function [H, V] = Tensor_Global_Arnoldi(A, V_init, m)

% Input:    A:          4D tensor of size NxLxNxL
%           V_init:     Initialmatrix of size NxL
%           m:          maximal number of iterations
%
% Output:   V           Orthonormal basis
%           H           Hessenbergmatrix 

[N, L, ~, ~] = size(A);
V = zeros(N, L, m+1);

V(:,:,1) = V_init/norm(V_init, 'fro');
H = [];


for j = 1:m
    W = Tensor4D_Matrix_Prod_2N(A, V(:,:,j));

    for i = 1:j
        H(i,j) = trace(V(:,:,i)' * W);
        W = W - H(i,j) * V(:,:,i);
    end

    H(j+1, j) = norm(W, 'fro');
    if abs(H(j+1, j))<1e-12
        break
    end
    V(:,:, j+1) = W / H(j+1, j);
    
end
H = H(1:end-1, 1:end);
V = V(:,:,1:end-1);
end
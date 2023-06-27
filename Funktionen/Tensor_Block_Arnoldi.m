function [H_m1m, V_basis] = Tensor_Block_Arnoldi(A, V_init, m)

% Input:    A:          4D tensor of size NxLxNxL
%           V_init:     Initialtensor of size NxLxR
%           m:          maximal number of iterations
%
% Output:   V_basis     Orthonormal basis
%           H_block     Block-hessenbergmatrix 

[N, L, ~, ~] = size(A);
R = size(V_init, 3);

V_basis = zeros(N, L, R*(m+1));
H_m1m = zeros(R*(m+1), R*m);

[V_basis(:,:, 1:R), ~] = Tensor3D_Matrix_QR_Factor(V_init);

for j = 1:m
    W = Tensor4D_Tensor3D_Prod_2N(A, V_basis(:,:, ((j-1)*R + 1): (j*R)));
    for i = 1:j
        H_m1m( ((i-1)*R + 1) : (R*i), ((j-1)*R + 1) : (R*j)) = Tensor3D_Tensor3D_Prod_2N(permute(V_basis(:,:, ((i-1)*R + 1):(i*R)), [3,1,2]) , W);
        W = W - Tensor3D_Matrix_Prod_1N(V_basis(:,:, ((i-1)*R + 1):(i*R)), H_m1m( ((i-1)*R+1) : (R*i), ((j-1)*R+1) : (R*j)));
    end

    [V_basis(:,:, ((j)*R + 1): ((j+1)*R)), H_m1m( j*R+1 : (j+1)*R , (j-1)*R+1 : R*j)] = Tensor3D_Matrix_QR_Factor(W);
    
end
end
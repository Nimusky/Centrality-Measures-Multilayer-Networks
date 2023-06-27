% This function computes the approx. RSC centralities
% using the block Arnoldi process

function [RSC_approx] = Block_Arnoldi_RSC(AdjTensor, V, alpha, m)

[~, ~, R] = size(V);

[H_m1m, V_basis] = Tensor_Block_Arnoldi(AdjTensor, V, m);

V_m_basis = V_basis(:,:,1:R*m);
H_m = H_m1m(1:R*m,:);

rechte_seite = Tensor3D_Tensor3D_Prod_2N(permute(V_m_basis, [3,1,2]), V);
linke_seite = Tensor3D_Tensor3D_Prod_2N(permute(V, [3,1,2]), V_m_basis);

f_Hm_RSC = inv(eye(size(H_m)) - alpha*H_m) - eye(size(H_m));

mitte_RSC = Tensor2D_Matrix_Prod_1N(f_Hm_RSC, rechte_seite);
RSC_approx = Tensor2D_Matrix_Prod_1N(linke_seite, mitte_RSC);
RSC_approx = diag(RSC_approx);

end
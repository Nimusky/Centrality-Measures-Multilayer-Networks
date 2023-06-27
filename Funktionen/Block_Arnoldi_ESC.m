% This function computes the approx. ESC centralities
% using the block Arnoldi process

function [ESC_approx] = Block_Arnoldi_ESC(AdjTensor, V, beta, m)

[~, ~, R] = size(V);

[H_m1m, V_basis] = Tensor_Block_Arnoldi(AdjTensor, V, m);

V_m_basis = V_basis(:,:,1:R*m);
H_m = H_m1m(1:R*m,:);

rechte_seite = Tensor3D_Tensor3D_Prod_2N(permute(V_m_basis, [3,1,2]), V);
linke_seite = Tensor3D_Tensor3D_Prod_2N(permute(V, [3,1,2]), V_m_basis);

f_Hm_ESC = expm(H_m * beta) - eye(size(H_m));

mitte_ESC = Tensor2D_Matrix_Prod_1N(f_Hm_ESC, rechte_seite);
ESC_approx = Tensor2D_Matrix_Prod_1N(linke_seite, mitte_ESC);
ESC_approx = diag(ESC_approx);

end
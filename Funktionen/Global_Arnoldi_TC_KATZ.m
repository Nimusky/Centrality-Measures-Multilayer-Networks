% This function computes the approx. TC and Katz centralities
% using the global Arnoldi process

function [TC_approx,KATZ_approx] = Global_Arnoldi_TC_KATZ(AdjTensor, alpha, beta, m)

[nodes, Layers, ~, ~] = size(AdjTensor);

V = ones(nodes, Layers);
[H, V_m] = Tensor_Global_Arnoldi(AdjTensor, V, m);

H_TC = exmpp_0(H, beta);
H_KATZ = inv(eye(size(H)) - alpha*H ) - eye(size(H));

rechte_seite = zeros(nodes*Layers, 1);
rechte_seite(1,1) = norm(V, 'fro');

linke_seite_TC = Tensor3D_Matrix_Prod_1N(V_m, H_TC);

TC_approx = Tensor3D_Vektor_Prod_1N(linke_seite_TC, rechte_seite);
TC_approx = TC_approx(:);

linke_seite_KATZ = Tensor3D_Matrix_Prod_1N(V_m, H_KATZ);

KATZ_approx = Tensor3D_Vektor_Prod_1N(linke_seite_KATZ, rechte_seite);
KATZ_approx = KATZ_approx(:);
end
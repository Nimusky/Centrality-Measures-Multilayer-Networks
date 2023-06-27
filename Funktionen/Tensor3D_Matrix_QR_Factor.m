% This function performs a QR-Factorization for a given 3D Tensor
% and returns a tensor Q and matrix R, that form the QR-Factorization.

function [Q_tensor,R] = Tensor3D_Matrix_QR_Factor(V_tensor)

[N, L, K] = size(V_tensor);
V_matrix = reshape(V_tensor, [], K);
[Q,R] = qr(V_matrix,0);
Q_tensor = reshape(Q, N, L,[]);

end

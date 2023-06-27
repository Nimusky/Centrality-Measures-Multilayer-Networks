function [ESC] = Tensor_MNC_ESC(AdjTensor, beta)
% This function computes the ESC_MNC centralities of all nodes for the 
% adjacency tensor A of a multilayer-network with L layers and n nodes in each layer.
%
% Output: (nx1)-vector, where (jx1) corresponds to ESC_MNC(j) of node x_j.

[n_nodes, L_layers, ~, ~] = size(AdjTensor);
AdjMatrix = convert_4dTensor_To_Matrix(AdjTensor);

ESC = diag(expm(AdjMatrix * beta) - eye(size(AdjMatrix)) );
ESC(ESC < 0) = 0;

ESC = sum(reshape(ESC,[n_nodes,L_layers]),2);
end
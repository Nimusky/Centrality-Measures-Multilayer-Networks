function [RSC] = Tensor_JC_RSC(AdjTensor, alpha)
% This function computes the RSC_JC centralities of all nodes for the 
% adjacency tensor A of a multilayer-network with L layers and n nodes in each layer.
%
% Output: (nx1)-vector, where (jx1) corresponds to RSC_JC(j) of node x_j.

[n_nodes, L_layers, ~, ~] = size(AdjTensor);
AdjMatrix = convert_4dTensor_To_Matrix(AdjTensor);

RSC = diag( inv(eye(n_nodes*L_layers) - alpha*AdjMatrix) - eye(n_nodes*L_layers));
end
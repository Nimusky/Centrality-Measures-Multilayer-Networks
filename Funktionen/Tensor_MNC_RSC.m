function [RSC] = Tensor_MNC_RSC(AdjTensor, alpha)
% This function computes the RSC_MNC centralities of all nodes for the 
% adjacency tensor A of a multilayer-network with L layers and n nodes in each layer.
%
% Output: (nx1)-vector, where (jx1) corresponds to RSC_MNC(j) of node x_j.

[n_nodes, L_layers, ~, ~] = size(AdjTensor);
AdjMatrix = convert_4dTensor_To_Matrix(AdjTensor);

matRSC = diag( inv(eye(n_nodes*L_layers) - alpha*AdjMatrix) - eye(n_nodes*L_layers));
    
RSC = sum(reshape(matRSC,[n_nodes,L_layers]),2);
end
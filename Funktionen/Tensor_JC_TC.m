function [TC] = Tensor_JC_TC(AdjTensor, beta)
% This function computes the TC_JC centralities of all nodes for the 
% adjacency tensor A of a multilayer-network with L layers and n nodes in each layer.
%
% Output: (nx1)-vector, where (jx1) corresponds to TC_JC(j) of node x_j.

[n_nodes, L_layers, ~, ~] = size(AdjTensor);
E = ones(n_nodes, L_layers);
AdjMatrix = convert_4dTensor_To_Matrix(AdjTensor);

matTC = Tensor4D_Matrix_Prod_2N(convert_Matrix_To_4dTensor (expm(AdjMatrix * beta) - eye(size(AdjMatrix)), n_nodes, L_layers), E);

TC = matTC(:);
TC(TC < 0) = 0;
end
function [KC] = Tensor_JC_KATZ(AdjTensor, alpha)
% This function computes the KC_JC centralities of all nodes for the 
% adjacency tensor A of a multilayer-network with L layers and n nodes in each layer.
%
% Output: (nx1)-vector, where (jx1) corresponds to KC_JC(j) of node x_j.

[n_nodes, L_layers, ~, ~] = size(AdjTensor);
E = ones(n_nodes, L_layers);
AdjMatrix = convert_4dTensor_To_Matrix(AdjTensor);

matKC = Tensor4D_Matrix_Prod_2N( ...
    convert_Matrix_To_4dTensor ((inv(eye(n_nodes*L_layers) - alpha*AdjMatrix ) - eye(n_nodes*L_layers) ), n_nodes, L_layers ) ...
    , E);

KC = matKC(:);
end
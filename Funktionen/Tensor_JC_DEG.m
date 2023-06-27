function [DEG] = Tensor_JC_DEG(AdjTensor)
% This function computes the DEG_JC centralities of all nodes for the 
% symmetric adjacency tensor A of a multilayer-network with L layers and n nodes in each layer.
% In the special case, that A is NOT a symmetric tensor, this function
% computes the JC_DEG_OUT centralities.

% Output: (nx1)-vector, where (jx1) corresponds to DEG_JC(j) of node x_j.

[n_nodes, L_layers, ~, ~] = size(AdjTensor);

E = ones(n_nodes, L_layers);
matJC = Tensor4D_Matrix_Prod_2N(AdjTensor, E);

DEG = matJC(:);
end


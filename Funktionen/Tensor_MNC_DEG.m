function [DEG] = Tensor_MNC_DEG(AdjTensor)
% This function computes the DEG_MNC centralities of all nodes for the 
% symmetric adjacency tensor A of a multilayer-network with L layers and n nodes in each layer.
% In the special case, that A is NOT a symmetric tensor, this function
% computes the MNC_DEG_OUT centralities.

% Output: (nx1)-vector, where (jx1) corresponds to DEG_MNC(j) of node x_j.

[n_nodes, L_layers, ~, ~] = size(AdjTensor);

E = ones(n_nodes, L_layers);
matJC = Tensor4D_Matrix_Prod_2N(AdjTensor, E);

DEG = sum(reshape(matJC(:),[n_nodes,L_layers]),2);
end
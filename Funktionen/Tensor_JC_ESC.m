function [ESC] = Tensor_JC_ESC(AdjTensor, beta)
% This function computes the ESC_JC centralities of all nodes for the 
% adjacency tensor A of a multilayer-network with L layers and n nodes in each layer.
%
% Output: (nx1)-vector, where (jx1) corresponds to ESC_JC(j) of node x_j.

AdjMatrix = convert_4dTensor_To_Matrix(AdjTensor);


ESC = diag(expm(AdjMatrix * beta) - eye(size(AdjMatrix)));
ESC(ESC < 0) = 0;
end
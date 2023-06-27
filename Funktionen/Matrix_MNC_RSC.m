function [RSC] = Matrix_MNC_RSC(AdjMatrix, n_nodes, L_layers, alpha)
% This function computes the RSC_MNC centralities of all nodes for the 
% adjacency matrix A of a multilayer-network with L layers and n nodes in each layer.
%
% Output: (nx1)-vector, where (jx1) corresponds to RSC_MNC(j) of node j.

RSC = sum(reshape(diag(inv(eye(size(AdjMatrix))-alpha*AdjMatrix)),[n_nodes,L_layers]),2);
end

function [KC] = Matrix_MNC_KATZ(AdjMatrix, n_nodes, L_layers, alpha)
% This function computes the KC_MNC centralities of all nodes for the 
% adjacency matrix A of a multilayer-network with L layers and n nodes in each layer.
%
% Output: (nx1)-vector, where (jx1) corresponds to KC_MNC(j) of node j.

KC = sum(reshape(sum(inv(eye(size(AdjMatrix))-alpha*AdjMatrix),2),[n_nodes,L_layers]),2);
end

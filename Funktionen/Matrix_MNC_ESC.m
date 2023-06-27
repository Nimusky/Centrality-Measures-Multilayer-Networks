function [ESC] = Matrix_MNC_ESC(AdjMatrix, n_nodes, L_layers, beta)
% This function computes the ESC_MNC centralities of all nodes for the 
% adjacency matrix A of a multilayer-network with L layers and n nodes in each layer.
%
% Output: (nx1)-vector, where (jx1) corresponds to ESC_MNC(j) of node j.

ESC = sum(reshape(diag(expm(beta*AdjMatrix)),[n_nodes,L_layers]),2);
end

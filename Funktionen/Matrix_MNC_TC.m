function [TC] = Matrix_MNC_TC(AdjMatrix, n_nodes, L_layers, beta)
% This function computes the TC_MNC centralities of all nodes for the 
% adjacency matrix A of a multilayer-network with L layers and n nodes in each layer.
%
% Output: (nx1)-vector, where (jx1) corresponds to TC_MNC(j) of node j.

TC = sum(reshape(sum(expm(beta*AdjMatrix),2),[n_nodes,L_layers]),2);
end

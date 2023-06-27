function [DEG] = Matrix_MNC_DEG(AdjMatrix, n_nodes, L_layers)
% This function computes the DEG_MNC centralities of all nodes for the 
% symmetric adjacency matrix A of a multilayer-network with L layers and n nodes in each layer.
% In the special case, that A is NOT a symmetric matrix, this function
% computes the DEG_OUT_MNC centralities, meaning:
%   e_i * A * [1,...,1]^T   for i=1,...,n.
%
% Output: (nx1)-vector, where (jx1) corresponds to DEG_MNC(j) of node x_j.

%   e_i * A * (e_j)^T           <=>  edge from x_i to x_j
%   sum(A, 1) = [1,...,1] * A   <=>  DEG_In  of all x_j for j=1,...,n
%   sum(A, 2) = A * [1,..,1]^T  <=>  DEG_Out of all x_i for i=1,...,n

DEG = sum(reshape(sum(AdjMatrix,2), [n_nodes,L_layers]),2);
end


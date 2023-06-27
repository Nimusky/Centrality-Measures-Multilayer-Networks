function [KC] = Agg_Matrix_KATZ(AggMatrix, alpha_Agg)
% This function computes the KC centralities of all nodes for the 
% aggregated adjacency matrix A of a multilayer-network.
%
% Output: (nx1)-vector, where (jx1) corresponds to KC(j) of node x_j.

KC = sum(inv(eye(size(AggMatrix))-alpha_Agg*AggMatrix),2);
end

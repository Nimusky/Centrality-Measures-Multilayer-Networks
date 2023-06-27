function [RSC] = Agg_Matrix_RSC(AggMatrix, alpha_Agg)
% This function computes the RSC centralities of all nodes for the 
% aggregated adjacency matrix A of a multilayer-network.
%
% Output: (nx1)-vector, where (jx1) corresponds to RSC(j) of node x_j.

RSC = diag(inv(eye(size(AggMatrix))-alpha_Agg*AggMatrix));
end

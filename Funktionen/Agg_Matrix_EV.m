function [EV] = Agg_Matrix_EV(AggMatrix)
% This function computes the EV centralities of all nodes for the 
% aggregated adjacency matrix A of a multilayer-network.
%
% Output: (nx1)-vector, where (jx1) corresponds to EV(j) of node x_j.

[largestEigenvector,~] = eigs(AggMatrix,1);

EV = largestEigenvector;
end

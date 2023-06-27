function [EV] = Matrix_JC_EV(AdjMatrix)
% This function computes the EV_JC centralities of all nodes for the 
% adjacency matrix A of a multilayer-network.
%
% Output: (nx1)-vector, where (jx1) corresponds to EV_JC(j) of node j.

[largestEigenvector,~] = eigs(AdjMatrix,1);

EV = largestEigenvector;
end


function [KC] = Matrix_JC_KATZ(AdjMatrix, alpha)
% This function computes the KC_JC centralities of all nodes for the 
% adjacency matrix A of a multilayer-network.
%
% Output: (nx1)-vector, where (jx1) corresponds to KC_JC(j) of node j.


KC = sum(inv(eye(size(AdjMatrix))-alpha*AdjMatrix),2);
end

function [ESC] = Matrix_JC_ESC(AdjMatrix, beta)
% This function computes the ESC_JC centralities of all nodes for the 
% adjacency matrix A of a multilayer-network.
%
% Output: (nx1)-vector, where (jx1) corresponds to ESC_JC(j) of node j.

ESC = diag(expm(beta*AdjMatrix));
end

function [TC] = Matrix_JC_TC(AdjMatrix, beta)
% This function computes the TC_JC centralities of all nodes for the 
% adjacency matrix A of a multilayer-network.
%
% Output: (nx1)-vector, where (jx1) corresponds to TC_JC(j) of node j.

TC = sum(expm(beta*AdjMatrix),2);
end

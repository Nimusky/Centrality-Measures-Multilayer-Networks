function [RSC] = Matrix_JC_RSC(AdjMatrix, alpha)
% This function computes the RSC_JC centralities of all nodes for the 
% adjacency matrix A of a multilayer-network.
%
% Output: (nx1)-vector, where (jx1) corresponds to RSC_JC(j) of node j.


RSC = diag(inv(eye(size(AdjMatrix)) - alpha*AdjMatrix));
end

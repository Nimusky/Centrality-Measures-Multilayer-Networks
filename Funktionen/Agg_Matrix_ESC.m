function [ESC] = Agg_Matrix_ESC(AggMatrix, beta_Agg)
% This function computes the ESC centralities of all nodes for the 
% aggregated adjacency matrix A of a multilayer-network.
%
% Output: (nx1)-vector, where (jx1) corresponds to ESC(j) of node x_j.


ESC = diag(expm(beta_Agg*AggMatrix));
end

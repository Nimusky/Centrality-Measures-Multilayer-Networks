function [TC] = Agg_Matrix_TC(AggMatrix, beta_Agg)
% This function computes the TC centralities of all nodes for the 
% aggregated adjacency matrix AggMatrix of a multilayer-network.
%
% Output: (nx1)-vector, where (jx1) corresponds to TC(j) of node x_j.


TC = sum(expm(beta_Agg*AggMatrix),2);
end

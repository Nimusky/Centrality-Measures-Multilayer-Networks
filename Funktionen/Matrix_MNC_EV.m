function [EV] = Matrix_MNC_EV(AdjMatrix, n_nodes, L_layers)
% This function computes the EV_MNC centralities of all nodes for the 
% adjacency matrix A of a multilayer-network with L layers and n nodes in each layer.
%
% Output: (nx1)-vector, where (jx1) corresponds to EV_MNC(j) of node j.

% Compute the eigenvalues and eigenvectors of the matrix A
[eigen_vectors, eigen_values] = eig(AdjMatrix);

% Find the index of the largest eigenvalue
[~, index] = max(diag(eigen_values));

% Find the eigenvector corresponding to the largest eigenvalue
eigenvector_max = eigen_vectors(:, index);

EV = sum(reshape(eigenvector_max, [n_nodes,L_layers]), 2);
end


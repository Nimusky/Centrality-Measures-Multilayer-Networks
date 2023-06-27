function A_agg = convert_AdjMatrix_to_AggMatrix(A, n_nodes, L_layers)
% This function creates an aggregated adjacency matrix for a given
% supraadjacency matrix A with n_nodes and L_layers

A_agg = zeros(n_nodes);
    
for i = 1:L_layers
    for j = 1:L_layers
        block = A(1+(i-1)*n_nodes : i*n_nodes, 1+(j-1)*n_nodes : j*n_nodes);
        A_agg = A_agg + block;
    end
end

end

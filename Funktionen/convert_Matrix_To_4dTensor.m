function [Tensor] = convert_Matrix_To_4dTensor(Matrix, n_nodes,L_layers)
% This function converts a (nL x nL) Matrix to a (n x L x n x L) Tensor

% create a null 4D Tensor of size n_node x L_layers x n_nodes x L_layers
Tensor = zeros(n_nodes, L_layers, n_nodes, L_layers);

% Iteration from left to right:
for l_2 = 1:L_layers
    for n_2 = 1:n_nodes
        for l_1 = 1:L_layers
            for n_1 = 1:n_nodes
                Tensor(n_1, l_1, n_2, l_2) ...
                    = Matrix( (n_1 + (l_1 -1)*n_nodes), (n_2 + (l_2 -1)*n_nodes) );
            end
        end
    end
end
%   AdjTensor(1, 1, 1, 1)   = AdjMatrix(1, 1)
%   ...                     = ...
%   AdjTensor(n, 1, 1, 1)   = AdjMatrix(n, 1)

%   AdjTensor(1, 2, 1, 1)   = AdjMatrix(n + 1, 1)
%   AdjTensor(2, 2, 1, 1)   = AdjMatrix(n + 2, 1)
%   ...                     = ...
%   AdjTensor(n, 2, 1, 1)   = AdjMatrix(2n, 1)

%   ...                     = ...
%   AdjTensor(n, L, 1, 1)   = AdjMatrix(nL, 1)

%   ...                     = ...
%   AdjTensor(n, L, n, L)   = AdjMatrix(nL, nL)
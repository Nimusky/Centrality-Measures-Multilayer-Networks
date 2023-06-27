function [Matrix] = convert_4dTensor_To_Matrix(Tensor)
% This function converts a (n x L x n x L) Tensor to a (nL x nL) Matrix

[n_nodes, L_layers, ~, ~] = size(Tensor);

% create a null Matrix of size (nL x nL)
Matrix = zeros(n_nodes*L_layers, n_nodes*L_layers);

% Iteration from left to right:
for l_2 = 1:L_layers
    for n_2 = 1:n_nodes
        for l_1 = 1:L_layers
            for n_1 = 1:n_nodes
                Matrix( (n_1 + (l_1 -1)*n_nodes), (n_2 + (l_2 -1)*n_nodes) ) ...
                    = Tensor(n_1, l_1, n_2, l_2);
            end
        end
    end
end

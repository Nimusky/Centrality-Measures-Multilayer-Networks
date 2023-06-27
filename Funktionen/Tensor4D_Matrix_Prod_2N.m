% This function performs a 4DTensor-Matrix-Multiplication for a 
% tensor of size (N_1 x L_1 x N_2 x L_2) and matrix of size (N_2 x L_2)
% and returns a tensor/matrix of size (N_1 x L_1)
%
% (N_1 x L_1 x N_2 x L_2) *_2 (N_2 x L_2) = (N_1 x L_1)

function result = Tensor4D_Matrix_Prod_2N(tensor, matrix)

    [N_1, L_1, N_2, L_2] = size(tensor);
    result = zeros(N_1, L_2);
    
    for n_1 = 1:N_1
        for l_1 = 1:L_1

            for n_2 = 1:N_2
                for l_2 = 1:L_2
                    result(n_1, l_1) = result(n_1, l_1) + tensor(n_1, l_1, n_2, l_2) * matrix(n_2, l_2);
                end
            end
        end
    end
end

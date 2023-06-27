% This function performs a 3DTensor-3DTensor-Multiplication for a 
% tensor of size (N_1 x L_1 x R) and tensor of size (R x N_2 x L_2)
% and returns a tensor of size (N_1 x L_1 x N_2 x L_2)
%
% (N_1 x L_1 x R) *_1 (R x N_2 x L_2) = (N_1 x L_1 x N_2 x L_2)

function result = Tensor3D_Tensor3D_Prod_1N(tensor_1, tensor_2)
    [N_1, L_1, R] = size(tensor_1);
    [~, N_2, L_2] = size(tensor_2);
    result = zeros(N_1, L_1, N_2, L_2);
    
    for n_1 = 1:N_1
        for l_1 = 1:L_1
            for n_2 = 1:N_2
                for l_2 = 1:L_2

                    for r = 1:R
                        result(n_1, l_1, n_2, l_2) = result(n_1, l_1, n_2, l_2) + tensor_1(n_1, l_1, r) * tensor_2(r, n_2, l_2);
                    end
                end
            end
        end
    end
end
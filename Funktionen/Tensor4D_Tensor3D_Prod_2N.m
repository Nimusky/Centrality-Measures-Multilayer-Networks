% This function performs a 4DTensor-3DTensor-Multiplication for a 
% tensor of size (N_2 x L_2 x N_1 x L_1) and tensor of size (N_1 x L_1 x R)
% and returns a tensor of size (N_2 x L_2 x R)
%
% (N_2 x L_2 x N_1 x L_1) *_2 (N_1 x L_1 x R) = (N_2 x L_2 x R)

function result = Tensor4D_Tensor3D_Prod_2N(tensor4D, tensor3D)
    [N_2, L_2, ~, ~] = size(tensor4D);
    [N_1, L_1, R] = size(tensor3D);

    result = zeros(N_2, L_2, R);
    
    for n_2 = 1:N_2
        for l_2 = 1:L_2
            for r = 1:R

                for n_1 = 1:N_1
                    for l_1 = 1:L_1

                        result(n_2, l_2,r) = result(n_2, l_2,r) + tensor4D(n_2, l_2, n_1, l_1) * tensor3D(n_1, l_1,r);
                    end
                end
            end
        end
    end
end
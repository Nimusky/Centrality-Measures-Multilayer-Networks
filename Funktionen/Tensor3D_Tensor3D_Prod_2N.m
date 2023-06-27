% This function performs a 3DTensor-3DTensor-Multiplication for a 
% tensor of size (R_1 x N x L) and matrix of size (N x L x R_2)
% and returns a tensor/matrix of size (R_1 x R_2)
%
% (R_1 x N x L) *_2 (N x L x R_2) = (R_1 x R_2)

function result = Tensor3D_Tensor3D_Prod_2N(tensor_1, tensor_2)

    [R_1, N, L] = size(tensor_1);
    [~, ~, R_2] = size(tensor_2);
    result = zeros(R_1, R_2);
    
    for r_1 = 1:R_1
        for r_2 = 1:R_2

            for n = 1:N
                for l = 1:L
                    result(r_1, r_2) = result(r_1, r_2) + tensor_1(r_1, n, l) * tensor_2(n, l, r_2);
                end
            end

        end
    end
end
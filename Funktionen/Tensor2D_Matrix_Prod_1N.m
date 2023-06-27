% This function performs a 2DTensor-Matrix-Multiplication for a 
% tensor of size (N x L) and matrix of size (L x K)
% and returns a matrix/tensor of size (N x K)
%
% (N x L) *_1 (L x K) = (N x K)

function result = Tensor2D_Matrix_Prod_1N(tensor, matrix)

    [N, L] = size(tensor);
    [~, K] = size(matrix);
    result = zeros(N, K);
    
    for n = 1:N
        for k = 1:K

            for l = 1:L
                result(n,k) = result(n,k) + tensor(n, l) * matrix(l, k);
            end
        end
    end
end
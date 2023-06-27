% This function performs a 3DTensor-Matrix-Multiplication for a 
% tensor of size (M x N x L) and matrix of size (N x L)
% and returns a vector of size (M x 1)
%
% (M x N x L) *_2 (N x L) = (M x 1)

function result = Tensor3D_Matrix_Prod_2N(tensor, matrix)

    [M, N, L] = size(tensor);
    result = zeros(M,1);
    
    for m = 1:M
        
        for n = 1:N
            for l = 1:L
                    result(m, 1) = result(m, 1) + tensor(m, n, l) * matrix(n, l);
            end
        end
    end
end

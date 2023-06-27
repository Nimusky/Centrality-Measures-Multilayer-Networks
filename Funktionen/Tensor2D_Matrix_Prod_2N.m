% This function performs a 2DTensor-Matrix-Multiplication for a 
% tensor of size (N x L) and matrix of size (N x L)
% and returns a scalar m
%
% (N L) *_2 (N L) = m

function result = Tensor2D_Matrix_Prod_2N(tensor, matrix)

    [N, L] = size(tensor);
    result = 0;
    
    for n = 1:N
        for l = 1:L
            result = result + tensor(n, l) * matrix(n, l);
        end
    end
end
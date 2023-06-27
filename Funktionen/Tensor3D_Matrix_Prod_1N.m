% This function performs a 3DTensor-Matrix-Multiplication for a 
% tensor of size (N x L x M_1) and matrix of size (M_1 x M_2)
% and returns a tensor of size (N x L x M_2)
%
% (N x L x M_1) *_1 (M_1 x M_2) = (N x L x M_2)
 
function result = Tensor3D_Matrix_Prod_1N(tensor, matrix)

    [N, L, M_1] = size(tensor);
    [~,M_2] = size(matrix);
    result = zeros(N,L,M_2);
    
    for n = 1:N
        for l = 1:L
            for m_2 = 1:M_2

                for m_1 = 1:M_1
                    result(n,l,m_2) = result(n,l,m_2) + tensor(n, l, m_1) * matrix(m_1, m_2);
                end
            end
        end
    end
end

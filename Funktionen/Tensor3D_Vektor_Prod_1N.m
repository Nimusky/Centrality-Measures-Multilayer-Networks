% This function performs a 3DTensor-Vector-Multiplication for a 
% tensor of size (N x L x M) and vector of size (M x 1)
% and returns a tensor/matrix of size (N x L)
%
% (N x L x M) *_1 (M) = (N x L)

function result = Tensor3D_Vektor_Prod_1N(tensor, vector)


    [N, L, M] = size(tensor);
    result = zeros(N,L);
    
    for n = 1:N
        for l = 1:L
            
            for m = 1:M
                    result(n,l) = result(n,l) + tensor(n, l, m) * vector(m);
            end
        end
    end
end


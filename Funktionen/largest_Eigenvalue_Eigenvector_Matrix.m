function [largestEigenvalue,largestEigenvector] = largest_Eigenvalue_Eigenvector_Matrix(A)
% This function computes the largest eigenvalue and 
% the corresponding eigenvector of the matrix A.

[largestEigenvector,largestEigenvalue] = eigs(A,1);
end

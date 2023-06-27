%% Acknowledgements:
%  The network used for this test is from:
%
%  "Navigability of interconnected networks under random failures",
%   Manlio De Domenico, Albert Solé-Ribalta, Sergio Gómez, and Alex 
%   Arenas PNAS 2014 111 (23) 8351-8356. 
%   Available under: http://www.pnas.org/content/111/23/8351.abstract

addpath('../Funktionen')


%% Convert the data from the network to a adjacency matrix, adjacency tensor and agg. matrix
nodes = 369;
Layers = 3;
AdjMatrix = convert_LnnW_Nodes0_Data_To_SymmAdjMatrix(nodes, Layers, 'london_transport_multiplex.edges');

% Choose Interlayer-Coupling
interlayer_adjacency_matrix = [0,1,1; 1,0,1; 1,1,0];
omega = 0;

AdjMatrix = AdjMatrix + omega*kron(interlayer_adjacency_matrix, eye(nodes));
AdjTensor = convert_Matrix_To_4dTensor(AdjMatrix, nodes, Layers);
Agg = convert_AdjMatrix_to_AggMatrix(AdjMatrix, nodes, Layers);


%% Compute the largest eigenvalue of Agg. and some alpha / beta
[lmax_Agg, ~] = largest_Eigenvalue_Eigenvector_Matrix(Agg); 
alpha_Agg = 0.5 / lmax_Agg;
beta_Agg = 0.5;


%% Compute the aggregated centralities
Agg_KATZ = Agg_Matrix_KATZ(Agg, alpha_Agg);
Agg_RSC = Agg_Matrix_RSC(Agg, alpha_Agg);
Agg_ESC = Agg_Matrix_ESC(Agg, beta_Agg);
Agg_TC = Agg_Matrix_TC(Agg, beta_Agg);


%% Compute the largest eigenvalue of AdjMatrix. and some alpha / beta
[lmax_AdjMatrix, ~] = largest_Eigenvalue_Eigenvector_Matrix(AdjMatrix); 
alpha_AdjMatrix = 0.5 / lmax_AdjMatrix;
beta_AdjMatrix = 0.5;


%% Compute the exact matrix-based joint centralities
JC_KATZ_Matrix = Matrix_JC_KATZ(AdjMatrix, alpha_AdjMatrix);
JC_RSC_Matrix = Matrix_JC_RSC(AdjMatrix, alpha_AdjMatrix);
JC_ESC_Matrix = Matrix_JC_ESC(AdjMatrix, beta_AdjMatrix);
JC_TC_Matrix = Matrix_JC_TC(AdjMatrix, beta_AdjMatrix);


%% Compute the exact tensor-based joint centralities
JC_KATZ_Tensor = Tensor_JC_KATZ(AdjTensor, alpha_AdjMatrix);
JC_RSC_Tensor = Tensor_JC_RSC(AdjTensor, alpha_AdjMatrix);
JC_ESC_Tensor = Tensor_JC_ESC(AdjTensor, beta_AdjMatrix);
JC_TC_Tensor = Tensor_JC_TC(AdjTensor, beta_AdjMatrix);


%% Create a tolerance for the exponential-based centralities
tolerance = 1e-12;
closeToZeroMask = abs(JC_ESC_Tensor) < tolerance & JC_ESC_Tensor ~= 0;
JC_ESC_Tensor(closeToZeroMask) = 0;
closeToZeroMask = abs(JC_TC_Tensor) < tolerance & JC_TC_Tensor ~= 0;
JC_TC_Tensor(closeToZeroMask) = 0;
closeToOneIndices = abs(JC_ESC_Matrix - 1.0) < tolerance;
JC_ESC_Matrix(closeToOneIndices) = 1.0;
closeToOneIndices = abs(JC_TC_Matrix - 1.0) < tolerance;
JC_TC_Matrix(closeToOneIndices) = 1.0;


%% Compute the matrix-based marginal-node centralities
MNC_KATZ_Matrix = sum(reshape(JC_KATZ_Matrix, [nodes,Layers]),2);
MNC_RSC_Matrix = sum(reshape(JC_RSC_Matrix, [nodes,Layers]),2);
MNC_ESC_Matrix = sum(reshape(JC_ESC_Matrix, [nodes,Layers]),2);
MNC_TC_Matrix = sum(reshape(JC_TC_Matrix, [nodes,Layers]),2);


%% Compute the tensor-based marginal-node centralities
MNC_KATZ_Tensor = sum(reshape(JC_KATZ_Tensor, [nodes,Layers]),2);
MNC_RSC_Tensor = sum(reshape(JC_RSC_Tensor, [nodes,Layers]),2);
MNC_ESC_Tensor = sum(reshape(JC_ESC_Tensor, [nodes,Layers]),2);
MNC_TC_Tensor = sum(reshape(JC_TC_Tensor, [nodes,Layers]),2);


%% Compute the Kendall-coefficients and the correlations
CorrMatTen = Corr_CKS(JC_ESC_Matrix, JC_RSC_Matrix, JC_TC_Matrix, JC_KATZ_Matrix, ...
    JC_ESC_Tensor, JC_RSC_Tensor, JC_TC_Tensor, JC_KATZ_Tensor);

CorrMatAgg = Corr_CKS(MNC_ESC_Matrix, MNC_RSC_Matrix, MNC_TC_Matrix, MNC_KATZ_Matrix, ...
    Agg_ESC, Agg_RSC, Agg_TC, Agg_KATZ);

CorrTenAgg = Corr_CKS(MNC_ESC_Tensor, MNC_RSC_Tensor, MNC_TC_Tensor, MNC_KATZ_Tensor, ...
    Agg_ESC, Agg_RSC, Agg_TC, Agg_KATZ);


figure;
scatter(JC_TC_Matrix, JC_TC_Tensor, 'filled');
xlabel('Matrixbasierte Joint-TC');
ylabel('Tensorbasierte Joint-TC');

figure;
scatter(MNC_RSC_Matrix, Agg_RSC, 'filled');
xlabel('Matrixbasierte Marginal-Node-RSC');
ylabel('Agg. RSC');

figure;
scatter(MNC_TC_Matrix, Agg_TC, 'filled');
xlabel('Matrixbasierte Marginal-Node-TC');
ylabel('Agg. TC');



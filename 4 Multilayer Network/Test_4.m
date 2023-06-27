%%  The network used for this test is a self made synthetic network.

addpath('../Funktionen')


%% Convert the data from the network to a adjacency matrix, adjacency tensor and agg. matrix
load('matlab.mat');
nodes = 200;
Layers = 35;
nL = nodes * Layers;

AdjMatrix = zeros(nL);

table_size = size(Multilayer_Netzwerk);
rows = table_size(1);
for row = 1:rows
    AdjMatrix( (Multilayer_Netzwerk{row, 2} - 1)*nodes + Multilayer_Netzwerk{row, 1}, ...
        (Multilayer_Netzwerk{row, 4} - 1)*nodes + Multilayer_Netzwerk{row, 3}) = Multilayer_Netzwerk{row, 5};
end

AdjTensor = convert_Matrix_To_4dTensor(AdjMatrix, nodes, Layers);
Agg = convert_AdjMatrix_to_AggMatrix(AdjMatrix, nodes, Layers);

clear Multilayer_Netzwerk; clear rows; clear row; clear table_size;


%% Compute the largest eigenvalue of agg. and some alpha / beta
[lmax_Agg, ~] = largest_Eigenvalue_Eigenvector_Matrix(Agg); 
alpha_Agg = 0.9 / lmax_Agg;
beta_Agg = 0.01;


%% Compute the agg. centralities
Agg_KATZ = Agg_Matrix_KATZ(Agg, alpha_Agg);
Agg_RSC = Agg_Matrix_RSC(Agg, alpha_Agg);
Agg_ESC = Agg_Matrix_ESC(Agg, beta_Agg);
Agg_TC = Agg_Matrix_TC(Agg, beta_Agg);


%% Compute the largest eigenvalue of AdjMatrix and some alpha / beta
[lmax_AdjMatrix, ~] = largest_Eigenvalue_Eigenvector_Matrix(AdjMatrix); 
alpha_AdjMatrix = 0.98 / lmax_AdjMatrix;
beta_AdjMatrix = 0.6;


%% Compute the exact tensor-based joint centralities
JC_KATZ_Tensor_exact = Tensor_JC_KATZ(AdjTensor, alpha_AdjMatrix);
JC_RSC_Tensor_exact = Tensor_JC_RSC(AdjTensor, alpha_AdjMatrix);
JC_ESC_Tensor_exact = Tensor_JC_ESC(AdjTensor, beta_AdjMatrix);
JC_TC_Tensor_exact = Tensor_JC_TC(AdjTensor, beta_AdjMatrix);


%% Tolerance for the exponential-based centralities wasn't necessary for this network
%tolerance = 1e-12;
%closeToZeroMask = abs(JC_ESC_Tensor_exact) < tolerance & JC_ESC_Tensor_exact ~= 0;
%JC_ESC_Tensor_exact(closeToZeroMask) = 0;
%closeToZeroMask = abs(JC_TC_Tensor_exact) < tolerance & JC_TC_Tensor_exact ~= 0;
%JC_TC_Tensor_exact(closeToZeroMask) = 0;


%% Compute the exact tensor-based marginal-node centralities
MNC_KATZ_Tensor_exact = sum(reshape(JC_KATZ_Tensor_exact, [nodes,Layers]),2);
MNC_RSC_Tensor_exact = sum(reshape(JC_RSC_Tensor_exact, [nodes,Layers]),2);
MNC_ESC_Tensor_exact = sum(reshape(JC_ESC_Tensor_exact, [nodes,Layers]),2);
MNC_TC_Tensor_exact = sum(reshape(JC_TC_Tensor_exact, [nodes,Layers]),2);


%% Compute the approximated tensor-based joint centralities TC / Katz with Arnoldi
maxit=5:12;

JC_TC_Tensor_approx=zeros(nL,length(maxit));
JC_KATZ_Tensor_approx=zeros(nL,length(maxit));

bar=waitbar(0,'Computing the approx. TC and Katz centralities');
for j=1:length(maxit)
    [JC_TC_Tensor_approx(:,j),JC_KATZ_Tensor_approx(:,j)] = Global_Arnoldi_TC_KATZ(AdjTensor, alpha_AdjMatrix, beta_AdjMatrix, maxit(j));

    waitbar(j/length(maxit),bar);
end
close(bar)

maxit=5:12;
figure(1)
semilogy(maxit,corr(JC_KATZ_Tensor_exact, JC_KATZ_Tensor_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)

figure(2)
semilogy(maxit,corr(JC_TC_Tensor_exact, JC_TC_Tensor_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)


%% Compute the approximated tensor-based joint centralities ESC with Arnoldi for the R most important nodes
[Top_Values_JC_ESC, Top_Indices_JC_ESC] = sort(JC_ESC_Tensor_exact, 'descend');
R = 50;
V = zeros(nodes, Layers, R);

for k = 1:R
    x = floor(Top_Indices_JC_ESC(k) / nodes);
    y = Top_Indices_JC_ESC(k) - nodes * x;

    % V(N,L,R):   N + (L-1)*nodes = ID of indices from exact
    V(y, x+1, k) = 1;
    
end


maxit=2:9;

JC_ESC_Tensor_approx=zeros(R,length(maxit));

bar=waitbar(0,'Computing the approx. ESC centrality for the R most important nodes');
for j=1:length(maxit)
    [JC_ESC_Tensor_approx(:,j)] = Block_Arnoldi_ESC(AdjTensor, V, beta_AdjMatrix, maxit(j));
    waitbar(j/length(maxit),bar);
end
close(bar)

Top_Values_JC_ESC = Top_Values_JC_ESC(1:R);

maxit=2:9;
figure(3)
semilogy(maxit,corr(Top_Values_JC_ESC, JC_ESC_Tensor_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")


%% Compute the approximated tensor-based joint centralities RSC with Arnoldi for the R most important nodes
[Top_Values_JC_RSC, Top_Indices_JC_RSC] = sort(JC_RSC_Tensor_exact, 'descend');
R = 50;
V = zeros(nodes, Layers, R);

for k = 1:R
    x = floor(Top_Indices_JC_RSC(k) / nodes);
    y = Top_Indices_JC_RSC(k) - nodes * x;

    % V(N,L,R):   N + (L-1)*nodes = ID of indices from exact
    V(y, x+1, k) = 1;
end

maxit=2:9;

JC_RSC_Tensor_approx=zeros(R,length(maxit));

bar=waitbar(0,'Computing the approx. RSC centrality for the R most important nodes');
for j=1:length(maxit)
    [JC_RSC_Tensor_approx(:,j)] = Block_Arnoldi_RSC(AdjTensor, V, alpha_AdjMatrix, maxit(j));
    waitbar(j/length(maxit),bar);
end
close(bar)

Top_Values_JC_RSC = Top_Values_JC_RSC(1:R);

maxit=2:9;
figure(4)
semilogy(maxit,corr(Top_Values_JC_RSC, JC_RSC_Tensor_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")


%% Compute the Kendall-coefficients and the correlations
CorrTenAgg = Corr_CKS(MNC_ESC_Tensor_exact, MNC_RSC_Tensor_exact, MNC_TC_Tensor_exact, MNC_KATZ_Tensor_exact, ...
    Agg_ESC, Agg_RSC, Agg_TC, Agg_KATZ);


%% Acknowledgements:
%  The network used for this test is from:
%
%  "Multiplex social ecological network analysis reveals how social 
%   changes affect community robustness more than resource depletion",
%   J.A. Baggio, S.B. Burnsilver, A. Arenas, J.S. Magdanz, G.P. Kofinas, and M. De Domenico. 
%   Available under: https://github.com/manlius/Alaska


addpath('../Funktionen')
addpath('../gauss')


%% Convert the data from the network to a adjacency matrix / agg. matrix 
nodes = 218;
Layers = 36;
nL = nodes * Layers;
AdjMatrix = convert_nLnLW_Data_To_SymmAdjMatrix(nodes, Layers, 'Wainwright.edges');
Agg = convert_AdjMatrix_to_AggMatrix(AdjMatrix, nodes, Layers);


%% Compute the largest eigenvalue of Agg and some alpha / beta
[lmax_Agg, ~] = largest_Eigenvalue_Eigenvector_Matrix(Agg); 

alpha_Agg = 0.95 / lmax_Agg;
beta_Agg = 0.00003;


%% Compute the aggregated centralities
Agg_KATZ = Agg_Matrix_KATZ(Agg, alpha_Agg);
Agg_RSC = Agg_Matrix_RSC(Agg, alpha_Agg);
Agg_ESC = Agg_Matrix_ESC(Agg, beta_Agg);
Agg_TC = Agg_Matrix_TC(Agg, beta_Agg);


%% Compute the largest/smallest eigenvalue of AdjMatrix and some alpha / beta
lmin_AdjMatrix=eigs(AdjMatrix,1,'smallestreal');
lmax_AdjMatrix=eigs(AdjMatrix,1,'largestreal');
alpha_AdjMatrix = 0.95 / lmax_AdjMatrix;
beta_AdjMatrix = 0.0001;


%% Compute the exact matrix-based joint centralities
JC_KATZ_Matrix_exact = Matrix_JC_KATZ(AdjMatrix, alpha_AdjMatrix);
JC_RSC_Matrix_exact = Matrix_JC_RSC(AdjMatrix, alpha_AdjMatrix);
JC_ESC_Matrix_exact = Matrix_JC_ESC(AdjMatrix, beta_AdjMatrix);
JC_TC_Matrix_exact = Matrix_JC_TC(AdjMatrix, beta_AdjMatrix);


%% Create a tolerance for the exponential-based centralities
tolerance = 1e-12;
closeToOneIndices = abs(JC_ESC_Matrix_exact - 1.0) < tolerance;
JC_ESC_Matrix_exact(closeToOneIndices) = 1.0;
closeToOneIndices = abs(JC_TC_Matrix_exact - 1.0) < tolerance;
JC_TC_Matrix_exact(closeToOneIndices) = 1.0;


%% Compute the exact matrix-based marginal-node centralities
MNC_KATZ_Matrix = sum(reshape(JC_KATZ_Matrix_exact, [nodes,Layers]),2);
MNC_RSC_Matrix = sum(reshape(JC_RSC_Matrix_exact, [nodes,Layers]),2);
MNC_ESC_Matrix = sum(reshape(JC_ESC_Matrix_exact, [nodes,Layers]),2);
MNC_TC_Matrix = sum(reshape(JC_TC_Matrix_exact, [nodes,Layers]),2);


%% Compute the approximated matrix-based joint centralities TC / Katz with Lanczos
b = ones(nL, 1);
maxit=1:12;

JC_KATZ_Matrix_approx=zeros(nL,length(maxit));
JC_TC_Matrix_approx=zeros(nL,length(maxit));

bar=waitbar(0,'Computing the approx. matrix-based TC and Katz centralities');
for j=1:length(maxit)
    [T, Q] = lanczos_tridiag(AdjMatrix, b, maxit(j));
    [S, D] = eig(T);
    
    JC_KATZ_Matrix_approx(:,j) = Q * S * (inv(eye(size(D))-alpha_AdjMatrix*D)) * S' * Q' * b;
    JC_TC_Matrix_approx(:,j) = Q * S * expm(beta_AdjMatrix*D) * S' * Q' * b;
    waitbar(j/length(maxit),bar);
end
close(bar)

maxit=1:12;
figure(1)
semilogy(maxit,corr(JC_KATZ_Matrix_exact, JC_KATZ_Matrix_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)

figure(2)
semilogy(maxit,corr(JC_TC_Matrix_exact, JC_TC_Matrix_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)


%% Compute the approximated matrix-based joint centralities ESC / RSC with Gauss
maxit = 2:4;
JC_RSC_Matrix_approx=zeros(nL,length(maxit));
JC_ESC_Matrix_approx=zeros(nL,length(maxit));

bar=waitbar(0,'Computing the approx. matrix-based ESC and RSC centralities');
for j=1:length(maxit)
    for i=1:nL
        ei=zeros(nL,1); 
        ei(i)=1;
        T_jp1 = lanczos_tridiag_Gauss(AdjMatrix,ei,maxit(j));

        JC_RSC_Matrix_approx(i,j) = gauss_resolvent(T_jp1,alpha_AdjMatrix);
        JC_ESC_Matrix_approx(i,j) = gauss_subgraph(T_jp1,beta_AdjMatrix);
    end
    waitbar(j/length(maxit),bar);
end
close(bar)

maxit = 2:4;
figure(3)
semilogy(maxit,corr(JC_RSC_Matrix_exact, JC_RSC_Matrix_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)

figure(4)
semilogy(maxit,corr(JC_ESC_Matrix_exact, JC_ESC_Matrix_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)

%% Compute the Kendall-coefficients and the correlations
CorrMatAgg = Corr_CKS(MNC_ESC_Matrix, MNC_RSC_Matrix, MNC_TC_Matrix, MNC_KATZ_Matrix, ...
    Agg_ESC, Agg_RSC, Agg_TC, Agg_KATZ);


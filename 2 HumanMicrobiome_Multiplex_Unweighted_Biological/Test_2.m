%% Acknowledgements:
%  The network used for this test is from:
%
%  "Dynamics and associations of microbial community types across the human body",
%   T. Ding and P. Schloss Nature 2014 509, 357–360.
%   Available under: https://www.nature.com/articles/nature13178
%
%  "Spectral Entropies as Information-Theoretic Tools for Complex Network
%  Comparison", M. De Domenico and J. Biamonte Physical Review X 2016 6, 041062.
%  Available under: https://journals.aps.org/prx/abstract/10.1103/PhysRevX.6.041062


addpath('../Funktionen')
addpath('../gauss')


%% Convert the data from the network to a adjacency matrix / adjacency tensor 
nodes = 305;
Layers = 18;
nL = nodes * Layers;
AdjMatrix = convert_LnnW_Nodes1_Data_To_SymmAdjMatrix(nodes, Layers, 'HumanMicrobiome_multiplex.edges');
AdjTensor = convert_Matrix_To_4dTensor(AdjMatrix, nodes, Layers);


%% Compute the largest/smallest eigenvalue of AdjMatrix and some alpha / beta
lmin_AdjMatrix=eigs(AdjMatrix,1,'smallestreal');
lmax_AdjMatrix=eigs(AdjMatrix,1,'largestreal');
alpha_AdjMatrix = 0.9 / lmax_AdjMatrix;
beta_AdjMatrix = 0.3;


%% Compute the exact matrix-based joint centralities
JC_KATZ_Matrix_exact = Matrix_JC_KATZ(AdjMatrix, alpha_AdjMatrix);
JC_RSC_Matrix_exact = Matrix_JC_RSC(AdjMatrix, alpha_AdjMatrix);
JC_ESC_Matrix_exact = Matrix_JC_ESC(AdjMatrix, beta_AdjMatrix);
JC_TC_Matrix_exact = Matrix_JC_TC(AdjMatrix, beta_AdjMatrix);


%% Compute the exact tensor-based joint centralities
JC_KATZ_Tensor_exact = Tensor_JC_KATZ(AdjTensor, alpha_AdjMatrix);
JC_RSC_Tensor_exact = Tensor_JC_RSC(AdjTensor, alpha_AdjMatrix);
JC_ESC_Tensor_exact = Tensor_JC_ESC(AdjTensor, beta_AdjMatrix);
JC_TC_Tensor_exact = Tensor_JC_TC(AdjTensor, beta_AdjMatrix);


%% Create a tolerance for the exponential-based centralities
tolerance = 1e-12;
closeToZeroMask = abs(JC_ESC_Tensor_exact) < tolerance & JC_ESC_Tensor_exact ~= 0;
JC_ESC_Tensor_exact(closeToZeroMask) = 0;
closeToZeroMask = abs(JC_TC_Tensor_exact) < tolerance & JC_TC_Tensor_exact ~= 0;
JC_TC_Tensor_exact(closeToZeroMask) = 0;
closeToOneIndices = abs(JC_ESC_Matrix_exact - 1.0) < tolerance;
JC_ESC_Matrix_exact(closeToOneIndices) = 1.0;
closeToOneIndices = abs(JC_TC_Matrix_exact - 1.0) < tolerance;
JC_TC_Matrix_exact(closeToOneIndices) = 1.0;


%% Compute the approximated matrix-based joint centralities TC / Katz with Lanczos
b = ones(nL, 1);
maxit=5:10;

JC_KATZ_Matrix_approx=zeros(nL,length(maxit));
JC_TC_Matrix_approx=zeros(nL,length(maxit));

bar=waitbar(0,'Computing the approx. matrix-based TC and Katz centralities');
for j=1:length(maxit)
    [T, Q] = lanczos_tridiag(AdjMatrix, b, maxit(j));
    [S, D] = eig(T);

    JC_TC_Matrix_approx(:,j) = Q * S * expm(beta_AdjMatrix*D) * S' * Q' * b;
    JC_KATZ_Matrix_approx(:,j) = Q * S * (inv(eye(size(D))-alpha_AdjMatrix*D)) * S' * Q' * b;
    waitbar(j/length(maxit),bar);
end
close(bar)

maxit=5:10;
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
maxit=2:5;
JC_RSC_Matrix_approx=zeros(nL,length(maxit));
JC_RSC_Matrix_approx_lower_radau=zeros(nL,length(maxit));
JC_RSC_Matrix_approx_upper_radau=zeros(nL,length(maxit));

JC_ESC_Matrix_approx=zeros(nL,length(maxit));
JC_ESC_Matrix_approx_lower_radau=zeros(nL,length(maxit));
JC_ESC_Matrix_approx_upper_radau=zeros(nL,length(maxit));

bar=waitbar(0,'Computing the approx. matrix-based ESC and RSC centralities');
for j=1:length(maxit)
    for i=1:nL
        ei=zeros(nL,1); 
        ei(i)=1;
        T_jp1 = lanczos_tridiag_Gauss(AdjMatrix,ei,maxit(j));
    
        JC_RSC_Matrix_approx(i,j)=gauss_resolvent(T_jp1,alpha_AdjMatrix);
        JC_RSC_Matrix_approx_lower_radau(i,j)=gauss_radau_resolvent(T_jp1,alpha_AdjMatrix,lmin_AdjMatrix);
        JC_RSC_Matrix_approx_upper_radau(i,j)=gauss_radau_resolvent(T_jp1,alpha_AdjMatrix,lmax_AdjMatrix);
        
        JC_ESC_Matrix_approx(i,j)=gauss_subgraph(T_jp1,beta_AdjMatrix);
        JC_ESC_Matrix_approx_lower_radau(i,j)=gauss_radau_subgraph(T_jp1,beta_AdjMatrix,lmin_AdjMatrix);
        JC_ESC_Matrix_approx_upper_radau(i,j)=gauss_radau_subgraph(T_jp1,beta_AdjMatrix,lmax_AdjMatrix);
    end
    waitbar(j/length(maxit),bar);
end
close(bar)

maxit=2:5;
figure(3)
semilogy(maxit,corr(JC_RSC_Matrix_exact, JC_RSC_Matrix_approx, 'type', 'Kendall'),'LineWidth',2)
hold on
semilogy(maxit,corr(JC_RSC_Matrix_exact, JC_RSC_Matrix_approx_lower_radau, 'type', 'Kendall'),'LineWidth',2)
semilogy(maxit,corr(JC_RSC_Matrix_exact, JC_RSC_Matrix_approx_upper_radau, 'type', 'Kendall'),'LineWidth',2)
hold off
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)
legend('Gauß-Regel','Untere Gauß-Radau-Regel','Obere Gauß-Radau-Regel')

figure(4)
semilogy(maxit,corr(JC_ESC_Matrix_exact, JC_ESC_Matrix_approx, 'type', 'Kendall'),'LineWidth',2)
hold on
semilogy(maxit,corr(JC_ESC_Matrix_exact, JC_ESC_Matrix_approx_lower_radau, 'type', 'Kendall'),'LineWidth',2)
semilogy(maxit,corr(JC_ESC_Matrix_exact, JC_ESC_Matrix_approx_upper_radau, 'type', 'Kendall'),'LineWidth',2)
hold off
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)
legend('Gauß-Regel','Untere Gauß-Radau-Regel','Obere Gauß-Radau-Regel')


%% Compute the approximated tensor-based joint centralities TC / Katz with Arnoldi
maxit=1:8;

JC_TC_Tensor_approx=zeros(nL,length(maxit));
JC_KATZ_Tensor_approx=zeros(nL,length(maxit));

bar=waitbar(0,'Computing the approx. tensor-based TC and Katz centralities');
for j=1:length(maxit)
    [JC_TC_Tensor_approx(:,j),JC_KATZ_Tensor_approx(:,j)] = Global_Arnoldi_TC_KATZ(AdjTensor, alpha_AdjMatrix, beta_AdjMatrix, maxit(j));
    waitbar(j/length(maxit),bar);
end
close(bar)

maxit=1:8;
figure(5)
semilogy(maxit,corr(JC_KATZ_Tensor_exact, JC_KATZ_Tensor_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)
figure(6)
semilogy(maxit,corr(JC_TC_Tensor_exact, JC_TC_Tensor_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)


%% Compute the approximated tensor-based joint centralities ESC with Arnoldi for the R most important nodes
[Top_Values_JC_ESC, Top_Indices_JC_ESC] = sort(JC_ESC_Tensor_exact, 'descend');
R = 20;
V = zeros(nodes, Layers, R);

for k = 1:R
    x = floor(Top_Indices_JC_ESC(k) / nodes);
    y = Top_Indices_JC_ESC(k) - nodes * x;

    % V(N,L,R):   N + (L-1)*nodes = ID of indices from exact
    V(y, x+1, k) = 1;
end

maxit=2:4;

JC_ESC_Tensor_approx=zeros(R,length(maxit));

bar=waitbar(0,'Computing the approx. tensor-based ESC centrality');
for j=1:length(maxit)
    [JC_ESC_Tensor_approx(:,j)] = Block_Arnoldi_ESC(AdjTensor, V, beta_AdjMatrix, maxit(j));
    waitbar(j/length(maxit),bar);
end
close(bar)

Top_Values_JC_ESC = Top_Values_JC_ESC(1:R);

maxit=2:4;
figure(7)
semilogy(maxit,corr(Top_Values_JC_ESC, JC_ESC_Tensor_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)


%% Compute the approximated tensor-based joint centralities RSC with Arnoldi for the R most important nodes
[Top_Values_JC_RSC, Top_Indices_JC_RSC] = sort(JC_RSC_Tensor_exact, 'descend');
R = 20;
V = zeros(nodes, Layers, R);

for k = 1:R
    x = floor(Top_Indices_JC_RSC(k) / nodes);
    y = Top_Indices_JC_RSC(k) - nodes * x;

    % V(N,L,R):   N + (L-1)*nodes = ID of indices from exact
    V(y, x+1, k) = 1;
end

maxit=2:4;

JC_RSC_Tensor_approx=zeros(R,length(maxit));

bar=waitbar(0,'Computing the approx. tensor-based RSC centrality');
for j=1:length(maxit)
    [JC_RSC_Tensor_approx(:,j)] = Block_Arnoldi_RSC(AdjTensor, V, alpha_AdjMatrix, maxit(j));
    waitbar(j/length(maxit),bar);
end
close(bar)

Top_Values_JC_RSC = Top_Values_JC_RSC(1:R);

maxit=2:4;
figure(8)
semilogy(maxit,corr(Top_Values_JC_RSC, JC_RSC_Tensor_approx, 'type', 'Kendall'),'LineWidth',2)
xlabel('m Iterationen')
ylabel("Kendall'scher Koeffizient")
xticks(maxit)


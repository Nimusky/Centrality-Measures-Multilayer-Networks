function [AdjMatrix] = convert_LnnW_Nodes0_Data_To_SymmAdjMatrix(n_nodes,L_layers, DataName)
% This function returns a symmetric adjacency matrix A of size (nL x nL),
% that corresponds to a multiplex undirected network (with n nodes and L layers) 
% from 'DataName', where the node IDs go from 0 to n-1 
% and edges are only allowed between nodes of the same layer.
%
% Data Format: layerID nodeID nodeID weight

% create a null matrix of size (nL x nL)
AdjMatrix = zeros(n_nodes * L_layers);

% create a DataTable with the Informations from the File 'DataName',
%   so that we can use data directly from .edges files
opts = detectImportOptions(DataName, "FileType","text");
DataTable = readtable(DataName, opts);
clear opts;

% measure the number of rows in DataTable
table_size = size(DataTable);
rows = table_size(1);

% iterate all rows from DataTable and save the data in the Adjacency Matrix A
for row = 1:rows

    % Since the IDs of the nodes from the given network data start with 0, 
    % we add '+1' to the indices i,j of A(i,j) and have:
    %    A(i + 1, j + 1) = weight of edge (x_i, x_j)

    % Formula: A( ((l-1)*n + i + 1), ((l-1)*n + j + 1) ) = 
    %                                   = weight of edge (x^{l}_{i}, x^{l}_{j})

    % save weight of edge (x_i, x_j) in A
    AdjMatrix( (DataTable{row, 1} - 1)*n_nodes + DataTable{row, 2} + 1, ...
        (DataTable{row, 1} - 1)*n_nodes + DataTable{row, 3} + 1 ) = DataTable{row, 4};

    % the network is NOT directed, hence
    % save weight of edge (x_j, x_i) in A
    AdjMatrix( (DataTable{row, 1} - 1)*n_nodes + DataTable{row, 3} + 1, ...
        (DataTable{row, 1} - 1)*n_nodes + DataTable{row, 2} + 1 ) = DataTable{row, 4};

end
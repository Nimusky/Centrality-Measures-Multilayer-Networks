function [AdjMatrix] = convert_nLnLW_Data_To_SymmAdjMatrix(n_nodes,L_layers, DataName)
% This function returns an adjacency matrix A of size (nL x nL),
% that corresponds to a multilayer network (with n nodes and L layers) 
% from 'DataName', where the node IDs go from 1 to n.
%
% Data Format: nodeFrom layerFrom nodeTo layerTo weight

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

    % Formula: A( ((l_1-1)*n + i), ((l_2-1)*n + j) ) = 
    %                                   = weight of edge (x^{l_1}_{i}, x^{l_2}_{j})

    % save weight of edge (x^{l_1}_{i}, x^{l_2}_{j}) in A
    AdjMatrix( (DataTable{row, 2} - 1)*n_nodes + DataTable{row, 1}, ...
        (DataTable{row, 4} - 1)*n_nodes + DataTable{row, 3}) = DataTable{row, 5};


    AdjMatrix( (DataTable{row, 4} - 1)*n_nodes + DataTable{row, 3}, ...
        (DataTable{row, 2} - 1)*n_nodes + DataTable{row, 1}) = DataTable{row, 5};


end
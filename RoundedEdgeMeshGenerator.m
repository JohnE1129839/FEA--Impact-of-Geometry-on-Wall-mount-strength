clc; clear; close all;

%% --- User settings ---
stlFile = 'roundededgehole.STL';  % your STL file
Hmax = 2;                     % max mesh element size
geoOrder = 'linear';              % 'linear' or 'quadratic'

%% --- Step 1: Create PDE model and import STL ---
model = createpde();
importGeometry(model, stlFile);  % no scaling

%% --- Step 2: Generate mesh ---
msh = generateMesh(model, 'Hmax', Hmax, 'GeometricOrder', geoOrder);

%% --- Step 3: Access nodes and elements ---
nodes = msh.Nodes;       % 3 x N array of node coordinates
elements = msh.Elements; % 3 x M array of triangular element indices

%% --- Step 4: Visualize ---
figure;
pdegplot(model, 'EdgeLabels', 'on');
axis equal
title('Imported STL Geometry');

figure;
pdemesh(msh);
axis equal
title('Mesh from STL');

% Optional: display number of nodes/elements
fprintf('Number of nodes: %d\n', size(nodes,2));
fprintf('Number of elements: %d\n', size(elements,2));

nodes = msh.Nodes * 1e-3;        % 2xN or 3xN array
elements = msh.Elements;  % element connectivity

edgeIDs = [1,2,7,8,12,11,6,5];  % array of edges you want
edgeNodes = findNodes(msh, 'region', 'Edge', edgeIDs);

writematrix(nodes, 's1r_nodes.txt');
writematrix(elements, 's1r_elements.txt');
writematrix(edgeNodes, "s1r_edgenodes.txt");
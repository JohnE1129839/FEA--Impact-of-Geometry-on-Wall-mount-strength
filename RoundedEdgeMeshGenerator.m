clc; clear; close all;

filename = "c3_off4";
stlFile = "./STLs/" + filename + ".STL";  % your STL file
Hmax = 2;                     % max mesh element size
geoOrder = 'linear';              % 'linear' or 'quadratic'

model = createpde();
importGeometry(model, stlFile);  % no scaling


msh = generateMesh(model, 'Hmax', Hmax, 'GeometricOrder', geoOrder);


nodes = msh.Nodes;       % 3 x N array of node coordinates
elements = msh.Elements; % 3 x M array of triangular element indices


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

edgeIDs = [5];  % array of edges you want
edgeNodes = findNodes(msh, 'region', 'Edge', edgeIDs);

writematrix(nodes, "./Meshes/"+ filename + "_nodes.txt");
writematrix(elements, "./Meshes/" + filename + "_elements.txt");
writematrix(edgeNodes, "./Meshes/"+filename+ "_edgenodes.txt");
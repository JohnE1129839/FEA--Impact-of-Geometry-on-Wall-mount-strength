%Material Properties
E = 2.1e11; %N/m2
h = 5e-3; %5mm plate
v = 0.28;

nodes = readmatrix("./Meshes/s1r_nodes.txt");
elements = readmatrix("./Meshes/s1r_elements.txt");
semiNodes = load('./Meshes/s1r_edgeNodes.txt');  
%% 

n_nodes = size(nodes,2);
n_elems = size(elements,2);
globalStiffness = zeros(2*n_nodes, 2*n_nodes);
globalForce = zeros(2*n_nodes, 1);  % if not already defined

for j = 1:n_elems
    id1 = elements(1, j);
    id2 = elements(2, j);
    id3 = elements(3, j);
    n1 = nodes(:, id1);
    n2 = nodes(:, id2);
    n3 = nodes(:, id3);
    ke = stiffnessMatrixUniformBar(E, h, v, n1, n2, n3);
    idx = [id1*2-1 id1*2 id2*2-1 id2*2 id3*2-1 id3*2];
    globalStiffness(idx, idx) = globalStiffness(idx, idx) + ke;
end


tol = 1e-11;  % tolerance

% Node coordinates
x = nodes(1,:);
y = nodes(2,:);

% Degrees of freedom (both x and y) for Dirichlet BC
semiDOFs = reshape([2*semiNodes-1; 2*semiNodes], [], 1);

% Apply Dirichlet BCs
globalStiffness(semiDOFs,:) = 0;
globalStiffness(:,semiDOFs) = 0;
globalStiffness(semiDOFs,semiDOFs) = eye(length(semiDOFs));
globalForce(semiDOFs) = 0;


q_y = -10;   % N/m (downward uniform load)
q_x = 0;       % no horizontal traction
tol = 1e-11;    % tolerance for detecting y=0 boundary



% loop over all elements, find edges on y=0
for e = 1:n_elems
    elemNodes = elements(:, e);
    nodeCoords = nodes(:, elemNodes);
    yCoords = nodeCoords(2,:);
    
    % find which nodes lie on y = 0
    onBottom = abs(yCoords - 0) < tol;
    
    if sum(onBottom) == 2
        % this element has one bottom edge
        edgeNodes = elemNodes(onBottom);
        n1 = edgeNodes(1); n2 = edgeNodes(2);

        % get coordinates
        x1 = nodes(1, n1); y1 = nodes(2, n1);
        x2 = nodes(1, n2); y2 = nodes(2, n2);
        L = sqrt((x2 - x1)^2 + (y2 - y1)^2);

        % equivalent nodal forces for this edge
        f_edge = L/2 * [q_x; q_y; q_x; q_y];

        % assemble into global force vector
        globalForce([2*n1-1, 2*n1, 2*n2-1, 2*n2]) = ...
            globalForce([2*n1-1, 2*n1, 2*n2-1, 2*n2]) + f_edge;
    end
end

U = globalStiffness \ globalForce;
Ux = U(1:2:end);   % x-displacement for each node
Uy = U(2:2:end);   % y-displacement for each node
Umag = sqrt(Ux.^2 + Uy.^2);  % displacement magnitude

scale = 1e7;  % exaggerate deformation for visualization
p = nodes;
t = elements;

% Compute deformed coordinates
p_def = p + scale * [Ux'; Uy'];

figure; hold on;
triplot(t', p(1,:), p(2,:), 'k--', 'LineWidth', 0.5);      % undeformed
triplot(t', p_def(1,:), p_def(2,:), 'b', 'LineWidth', 1.2); % deformed
axis equal tight
title('Undeformed (black dashed) vs Deformed (blue)');
legend('Undeformed','Deformed');

%% --- Material properties (plane stress)
C = stressStrainMatrix(E,v);   % 3x3

%% --- Preallocate
n_elems = size(elements, 2);
elemVM = zeros(1, n_elems);   % von Mises stress per element
elemSigma = zeros(3, n_elems); % [sigma_x; sigma_y; tau_xy] per element

%% --- Loop over elements
for e = 1:n_elems
    id = elements(:, e);          % node indices
    x1 = nodes(1,id(1)); y1 = nodes(2,id(1));
    x2 = nodes(1,id(2)); y2 = nodes(2,id(2));
    x3 = nodes(1,id(3)); y3 = nodes(2,id(3));

    % Triangle area
    A = 0.5 * det([1 x1 y1; 1 x2 y2; 1 x3 y3]);

    % B matrix
    b1 = y2 - y3; c1 = x3 - x2;
    b2 = y3 - y1; c2 = x1 - x3;
    b3 = y1 - y2; c3 = x2 - x1;
    B = 1/(2*A) * [b1 0 b2 0 b3 0;
                    0 c1 0 c2 0 c3;
                    c1 b1 c2 b2 c3 b3];

    % Element displacement vector using your Ux and Uy
    ue = [Ux(id(1)); Uy(id(1)); Ux(id(2)); Uy(id(2)); Ux(id(3)); Uy(id(3))];

    % Strain and stress
    eps = B * ue;       % [epsilon_x; epsilon_y; gamma_xy]
    sigma = C * eps;    % [sigma_x; sigma_y; tau_xy]
    elemSigma(:, e) = sigma;

    % Von Mises stress
    sxx = sigma(1); syy = sigma(2); txy = sigma(3);
    elemVM(e) = sqrt(sxx^2 - sxx*syy + syy^2 + 3*txy^2);
end

%% --- Nodal averaging
n_nodes = size(nodes, 2);
nodalVM = zeros(1, n_nodes);
nodeCount = zeros(1, n_nodes);

for e = 1:n_elems
    id = elements(:, e);
    for k = 1:3
        nodalVM(id(k)) = nodalVM(id(k)) + elemVM(e);
        nodeCount(id(k)) = nodeCount(id(k)) + 1;
    end
end
nodalVM = nodalVM ./ nodeCount;

%% --- Plot nodal-averaged von Mises stress
figure;
pdeplot(msh, 'XYData', nodalVM, 'ColorMap', 'jet', ...
        'Deformation', [Ux Uy], 'DeformationScaleFactor', 100);
axis equal tight
title('Nodal-averaged von Mises Stress (Deformed Shape)');
colorbar

max(nodalVM)

%Uncomment this section if you do not have the pde addon
% figure;
% trisurf(elements', nodes(1,:), nodes(2,:), nodalVM, 'EdgeColor','none');
% view(2);              % top-down view for 2D
% axis equal;           % equal scaling
% colorbar;
% colormap('jet');
% title('Von Mises Stress');
% xlabel('X [m]');
% ylabel('Y [m]');
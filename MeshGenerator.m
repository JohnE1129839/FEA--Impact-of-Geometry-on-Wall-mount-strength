model = createpde();
%Creating mesh geometries
R1 = [3,4,0,0.1,0.1,0,0,0,0.1,0.1,0,0]';  % General 10cm x 10 cm rectangle

%Hole dimensions
xc = 0.05; %5cm
yc = 0.025;%2.5cm
r = 0.02; %2cm

C1 = [1,xc,yc,r,0,0,0,0,0,0,0,0]';  % Circular hole dimension

x1 = xc - r;  % left
x2 = xc + r;  % right
y1 = yc - r;  % bottom
y2 = yc + r;  % top
S1 = [3,4,x1,x2,x2,x1,y1,y1,y2,y2,0,0]'; %Square hole dimension

% Combine
gd = [R1 S1];
ns = char('R1','S1')';
sf = 'R1 - S1';   % subtract the circle from the rectangle
g = decsg(gd, sf, ns);

geometryFromEdges(model, g);
msh = generateMesh(model,'Hmax',0.002, 'GeometricOrder','linear'); % smaller Hmax = finer mesh

nodes = msh.Nodes;        % 2xN or 3xN array
elements = msh.Elements;  % element connectivity

writematrix(nodes, 's1_nodes.txt');
writematrix(elements, 's1_elements.txt');


%% Solve laplace's equation on a square region with a circle cut out
% Requires PDE Toolbox
clear all;

%% Specify the problem geometry
model = createpde(1);
R1 = [3, 4, -1, 1, 1, -1, 1, 1, -1, -1];
C1 = [1, 0, 0, 0.5]';
C1 = [C1;zeros(length(R1) - length(C1),1)];
R1 = R1';
gm = [R1, C1];
sf = 'R1-C1';
ns = char('R1', 'C1')';
g = decsg(gm, sf, ns);
geometryFromEdges(model, g);
figure(1);
pdegplot(model, 'EdgeLabels', 'on');
axis equal;

%% Specify boundary conditions
% Volatage of 1 on top edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 4, 'r', 1);
% Voltage of 2 on bottom edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 2, 'r', -1);
% Normal dervative must be 0 at all other boundaries
applyBoundaryCondition(model, 'neumann', 'Edge', [1,3,5,6,7,8], 'q', 0, 'g', 0);

%% Specify equation coefficients
% Equation is (del)^2 u = 0
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);

%% Generate mesh and solve equation
generateMesh(model);
result = solvepde(model);
u = result.NodalSolution;

%% Plot Result
figure(2);
pdeplot(model,'XYData',u,'ZData',u);
view(2);
colormap('hsv');
axis equal;
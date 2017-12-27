function [x, y, phi] = calc_potential(V0, res)
%calc_potential Solve Laplace's equation on a z-shaped conductor

% Specify geometry
R1 = [3 4 -3 0 0 -3 2.5 2.5 2.3 2.3]';
R2 = [3 4 -0.25 0.25 0.25 -0.25 2.5 2.5 -2.5 -2.5]';
R3 = [3 4 3 0 0 3 -2.5 -2.5 -2.3 -2.3]';
gm = [R1, R2, R3];
sf = 'R1+R2+R3';
ns = char('R1', 'R2', 'R3')';

model = createpde(1);
g = decsg(gm, sf, ns);
geometryFromEdges(model, g);

% Equation is (del)^2 phi = 0 (Laplace's equation)
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);

% Specify boundary conditions
% Volatage of V0 on top edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 1, 'r', 0);
% Voltage of -V0 on bottom edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 2, 'r', V0);
% Normal dervative must be 0 at all other boundaries
applyBoundaryCondition(model, 'neumann', 'Edge', 3:17, 'q', 0, 'g', 0);

% Solve equation
generateMesh(model, 'Hmax', 0.1);
result = solvepde(model);

% Interpolate solution
x = linspace(-5, 5, 10*res);
y = linspace(-5, 5, 10*res);
[x, y] = meshgrid(x, y);

phi = interpolateSolution(result, x, y);
phi = reshape(phi, size(x));

end


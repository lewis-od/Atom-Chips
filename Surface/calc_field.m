function [phi, xq, yq, Bx, By] = calc_field(params, gm, ns, sf)
%calc_field Calculates the magnetic field given the problem params and
%geometry

mu_0 = pi*4e-7; % Permeability of free space
sigma = params(1);
V0 = params(2);
d = params(3);
z = params(4);

% Set up problem geometry
model = createpde(1);
g = decsg(gm, sf, ns);
geometryFromEdges(model, g);
dim = size(g);
edges = 1:dim(2); % All edges
% Remove edges 2 and 4 as these have dirichlet boundary conditions
edges = edges(edges~=2);
edges = edges(edges~=4);

% Equation is (del)^2 phi = 0 (Laplace's equation)
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);

% Specify boundary conditions
% Volatage of V0 on top edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 4, 'r', V0);
% Voltage of -V0 on bottom edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 2, 'r', -V0);
% Normal dervative must be 0 at all other boundaries
applyBoundaryCondition(model, 'neumann', 'Edge', edges, 'q', 0, 'g', 0);

% Generate mesh and solve PDE
generateMesh(model);
result = solvepde(model);

% Assume geometry is given by rectangle then circles. Length of an array
% specifying a rectangle is 10, 2nd coord says how many points are given
n = gm(2);
rect_x = gm(3:2+n);
rect_y = gm(3+n:2+2*n);
max_x = max(rect_x);
max_y = max(rect_y);
min_x = min(rect_x);
min_y = min(rect_y);

% Interpolate solution from mesh onto linearly spaced grid
res = 100; % How many points per unit
xq = linspace(min_x, max_x, (max_x-min_x)*res);
yq = linspace(min_y, max_y, (max_y-min_y)*res);
[xq, yq] = meshgrid(xq, yq);

phi = interpolateSolution(result, xq, yq);
phi = reshape(phi, size(xq));

% Calculate electric field
[Ex, Ey] = gradient(phi);
Ex = -1.*Ex;
Ey = -1.*Ey;

% Remove infinities
Ex(isnan(Ex)) = 0.0;
Ey(isnan(Ey)) = 0.0;

% Ohm's law
jx = sigma.*Ex;
jy = sigma.*Ey;

% Calculate f(kx, ky, z, d) in Fourier space
kx = (0:1:length(jx));
ky = (0:1:length(jy));

[kx, ky] = meshgrid(kx, ky);
k = sqrt(kx.^2 + ky.^2);
f_hat = (mu_0/2)*d .* ((1-exp(-d*k))./(k*d)) .* exp(-k.*z);

f_hat(isnan(f_hat)) = 0.0; % Remove infinities for inverse transform
f = real(ifft2(f_hat)); % Calculate f in real space

% Calculate B using the convolution theorem
Bx = conv2(jy, f, 'same');
By = -conv2(jx, f, 'same');

end


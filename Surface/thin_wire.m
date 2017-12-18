%% Calculate the field for a thin 2D wire, and compare to the analytical
%  solution for a thin 1D wire

clear all;

%% Parameters
sigma = 5.96e7; % Conductivity of conductor [S*m^-1]
V0 = 1.0; % Voltage is V0 and -V0 at ends of conductor [V]
d = 0.1; % Thickness of conductor [m]
z = 0.5; % Distance from conductor [m]

%% Specify Geometry
% Specify geometry
R1 = [3, 4, -0.1, 0.1, 0.1, -0.1, 20, 20, -20, -20]';
gm = [R1];
sf = 'R1';
ns = char('R1')';

%% Calculate electric potential
% Set up problem geometry
model = createpde(1);
g = decsg(gm, sf, ns);
geometryFromEdges(model, g);

% Equation is (del)^2 phi = 0 (Laplace's equation)
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);

% Specify boundary conditions
% Volatage of V0 on top edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 1, 'r', V0);
% Voltage of -V0 on bottom edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 3, 'r', -V0);
% Normal dervative must be 0 at all other boundaries
applyBoundaryCondition(model, 'neumann', 'Edge', [2 4], 'q', 0, 'g', 0);

% Generate mesh and solve PDE
generateMesh(model, 'Hmax', 0.01);
result = solvepde(model);

% Assume geometry is given by rectangle. Length of an array specifying a
% rectangle is 10, 2nd coord says how many points are given
n = gm(2);
rect_x = gm(3:2+n);
rect_y = gm(3+n:2+2*n);
max_x = max(rect_x);
max_y = max(rect_y);
min_x = min(rect_x);
min_y = min(rect_y);

% Interpolate solution from mesh onto linearly spaced grid
res = 50; % How many points per unit
xq = linspace(min_x, max_x, (max_x-min_x)*res+1);
yq = linspace(min_y, max_y, (max_y-min_y)*res+1);
[xq, yq] = meshgrid(xq, yq);

phi = interpolateSolution(result, xq, yq);
phi = reshape(phi, size(xq));

xq = xq.*1e-6;
yq = yq.*1e-6;
z = z.*1e-6;
res = res.*1e6;

dims = size(xq);
x0 = ceil(dims(2)/2);
y0 = ceil(dims(1)/2);

% Extend x and y, then pad phi with zeros in the regions surrounding the
% conductor
% xq = linspace(-1, 1, 2*res+1);
% yq = linspace(-1.5, 1.5, 3*res+1);
% [xq, yq] = meshgrid(xq, yq);
% diff = size(xq) - size(phi);
% phi = padarray(phi, diff./2, 0, 'both');

%% Calculate magnetic field
[Bx, By] = calc_field(phi, sigma, res, d, z);

% Trim wire between y=-1 and y=1 so that the effects of the ends of the
% wire are negligible and the field can be compared to that of an infinite
% wire
n = res*1e-6; % Num points in interval of 1 micron
xq = xq(y0-n:y0+n, :);
yq = yq(y0-n:y0+n, :);
Bx = Bx(y0-n:y0+n, :);
By = By(y0-n:y0+n, :);

%% Calculate magnetic field analytically - 2D wire
mu_0 = 4e-7 * pi;
[Ex, Ey] = gradient(phi);
Ex = -Ex;
Ey = -Ey;
W = 0.2e-6;

J = sigma.*Ey; % Jx should be negligible
J = mean(mean(J)); % J should be uniform

Bx_analytic = ((mu_0*J)/(2*pi)).* (atan((xq+W/2)./z) - atan((xq-W/2)./z));

%% Calculate magnetic field analytically - 1D wire
r = sqrt(xq.^2 + z.^2);
I = J*W;

B_1d = ((mu_0*I)/2*pi)./r;

%% Plot results
subplot(1, 2, 1);
surf(xq, yq, Bx, 'EdgeColor', 'none', 'FaceColor', 'interp');
view(2);
axis equal;
colormap jet;
title('B_x (Computational)');
xlabel('x');
ylabel('y');

subplot(1, 2, 2);
surf(xq, yq, Bx_analytic, 'EdgeColor', 'none', 'FaceColor', 'interp');
view(2);
axis equal;
colormap jet;
title('B_x (Analytic)');
xlabel('x');
ylabel('y');
%% Solve laplace's equation on a square region with a circle cut out
% Requires PDE Toolbox
clear all;

%% Parameters
sigma = 5.96e7; % Conductivity of conductor [S*m^-1]
V0 = 1.0; % Voltage is V0 and -V0 at ends of conductor [V]
mu_0 = pi*4e-7; % Permeability of free space
d = 0.1; % Thickness of conductor [m]
z = 0.5; % Distance from conductor [m]

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
title('Problem Geometry', 'FontSize', 18);
xlabel('x', 'FontSize', 16);
ylabel('y', 'FontSize', 16);
axis equal;

%% Specify boundary conditions
% Volatage of 1 on top edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 4, 'r', V0);
% Voltage of 2 on bottom edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 2, 'r', -V0);
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
plot = pdeplot(model,'XYData',u,'ZData',u);
view(2);
title("Solution", 'FontSize', 18);
colormap('hsv');
xlabel('x', 'FontSize', 16);
ylabel('y', 'FontSize', 16);
plot(2).Label.String = "\phi (x,y) [V]";
plot(2).Label.FontSize = 16;
axis equal;

%% Interpolate from mesh to a square grid
xq = linspace(-1, 1, 200);
yq = xq;
[xq, yq] = meshgrid(xq, yq);
u_interp = interpolateSolution(result, xq, yq);
u_interp = reshape(u_interp, [length(xq), length(yq)]);

figure(3);
surf(xq, yq, u_interp, 'Mesh', 'no', 'FaceColor', 'interp');
view(2);
xlabel('x');
ylabel('y');
colorbar();
axis equal;

%% Calculate electric field
[Ex, Ey] = gradient(u_interp);

% Prepare for Fourier transform
Ex(isnan(Ex)) = 0.0;
Ey(isnan(Ey)) = 0.0;

%% Calculate magnetic field
% Ohm's law
jx = sigma.*Ex;
jy = sigma.*Ey;

% Calculate f(kx, ky, z, d) in Fourier space
kx = (0:1:length(jx));
ky = (0:1:length(jy));

[kx, ky] = meshgrid(kx, ky);
k = sqrt(kx.^2 + ky.^2);
f_hat = (mu_0/2)*d .* ((1-exp(-d*k))./(k*d)) .* exp(-k.*z);

f_hat(isnan(f_hat)) = 0.0; % Prepare for inverse transform
f = real(ifft2(f_hat)); % Calculate f in real space

%% Apply convolution theorem
Bx = conv2(jy, f, 'same');
By = -conv2(jx, f, 'same');

figure(4);
B = sqrt(Bx.^2 + By.^2);
surf(xq, yq, B, 'EdgeColor', 'none');

%% Calculate B field from 2 notches in a conductor
clear all;

%% Parameters
z = 2;
d = 0.25;
sigma = 44.2e6; % Conductivity of gold (S m^-1)

%% Specify Geometry
model = createpde(1);
R1 = [3, 4, -5, 5, 5, -5, 5, 5, -5, -5]';
R2 = [3, 4, -5, -2.5, -2.5, -5, 1.25, 1.25, 0.75, 0.75]';
R3 = [3, 4, 5, 2.5, 2.5, 5, -1.25, -1.25, -0.75, -0.75]';
gm = [R1, R2, R3];
sf = 'R1-R2-R3';
ns = char('R1', 'R2', 'R3')';
g = decsg(gm, sf, ns);
geometryFromEdges(model, g);

neumann_edges = 3:length(g);

%% Specify equation and boundary conditions
% Equation is (del)^2 phi = 0 (Laplace's equation)
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);

% Specify boundary conditions
% Volatage of V0 on top edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 1, 'r', 56.6e-3/2);
% Voltage of -V0 on bottom edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 2, 'r', -56.6e-3/2);
% Normal dervative must be 0 at all other boundaries
applyBoundaryCondition(model, 'neumann', 'Edge', neumann_edges, 'q', 0, 'g', 0);

%% Generate mesh and solve PDE
generateMesh(model, 'Hmax', 0.1);
result = solvepde(model);

%% Interpolate onto grid
res = 501;
xq = linspace(-5, 5, res);
yq = linspace(-5, 5, res);
[xq, yq] = meshgrid(xq, yq);

dims = size(xq);
x0 = ceil(dims(2)/2);
y0 = ceil(dims(1)/2);

phi = interpolateSolution(result, xq, yq);
phi = reshape(phi, size(xq));

figure();
surf(xq, yq, phi, 'EdgeColor', 'none', 'FaceColor', 'interp');
view(2);
c = colorbar();
c.Label.String = '\phi [V]';
c.Label.FontSize = 16;
c.FontSize = 16;
axis equal;
colormap jet;
title('Electric Potential', 'FontSize', 18);
xlabel('x [\mum]', 'FontSize', 16);
ylabel('y [\mum]', 'FontSize', 16);
hold on;
gplot = pdegplot(model);
gplot.Color = [0 0 0];

%% Calc and plot B field
[Bx, By] = calc_field(phi, sigma, d, z);
B = sqrt(Bx.^2 + By.^2);
figure();
surf(xq, yq, B, 'EdgeColor', 'none', 'FaceColor', 'interp');
view(2);
c = colorbar();
c.Label.String = '|B| [T]';
c.Label.FontSize = 16;
c.FontSize = 16;
colormap jet;
xlabel('x [\mum]', 'FontSize', 16);
ylabel('y [\mum]', 'FontSize', 16);
title(['Magnetic Field at z=' num2str(z) '\mum'], 'FontSize', 18);

B_slice_x = B(y0, :);
xq_slice_x = xq(y0, :);
B_slice_y = B(:, x0);
yq_slice_y = yq(:, x0);

figure();
title(['Magnetic Field at z=' num2str(z) '\mum'], 'FontSize', 18);

subplot(1, 2, 1);
plot(xq_slice_x, B_slice_x);
xlabel('x [\mum]', 'FontSize', 16);
ylabel('B [T]', 'FontSize', 16);
title('y = 0\mum', 'FontSize', 16);

subplot(1, 2, 2);
plot(yq_slice_y, B_slice_y);
xlabel('y [\mum]', 'FontSize', 16);
ylabel('B [T]', 'FontSize', 16);
title('x = 0\mum', 'FontSize', 16);

% zq = [2, 5, 10];
% figure();
% for i = 1:length(zq)
%     [Bx, By] = calc_field(phi, sigma, 0.25, zq(i));
%     B = sqrt(Bx.^2 + By.^2);
%     Bx_slice_y = Bx(y0, :);
%     By_slice_y = By(y0, :);
%     x_slice_y = xq(y0, :);
% 
%     subplot(1, 2, 1);
%     hold on;
%     plot(x_slice_y, Bx_slice_y);
%     
%     subplot(1, 2, 2);
%     hold on;
%     plot(x_slice_y, By_slice_y);
% end
% 
% subplot(1, 2, 1);
% xlabel('x [\mum]', 'FontSize', 14);
% ylabel('B_x [T]', 'FontSize', 14);
% legend(cellstr(num2str(zq', 'z=%.2f')));
% hold off;
% 
% subplot(1, 2, 2);
% xlabel('x [\mum]', 'FontSize', 14);
% ylabel('B_y [T]', 'FontSize', 14);
% legend(cellstr(num2str(zq', 'z=%.2f')));
% hold off;
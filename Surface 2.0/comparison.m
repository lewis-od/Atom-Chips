%% Compare the analytic expression for a 2D wire to the 2D FEM
clear all;

%% Parameters

% Number of points in each dimension to evaluate the field at
Ny = 20;
Nz = 20;

% Size of current elements to use
dx = 1e-5;
dy = 50e-7;

L = 1e-4; % Length of wire [m]
W = 50e-6; % Width of wire [m]
J = 100; % Current density [A/m]
mu_0 = 4e-7*pi; % Permeability of free space

%% Points to evaluate field at
xq = 0;
yq = linspace(-100e-6, 100e-6, Ny);
zq = linspace(-100e-6, 100e-6, Nz);
[yq, zq] = meshgrid(yq, zq);

%% Current density for 2D FEM
y = linspace(-L/2, L/2, 200);
x = linspace(-W/2, W/2, 200);
[x, y] = meshgrid(x, y);

Jx = J.*ones(size(x));
Jy = zeros(size(x));

%% Calculate field computationally

[Bx, By, Bz] = calc_field(x, y, Jx, Jy, dx, dy, zeros(size(yq)), yq, zq);

%% Calculate field from analytic expression
[By_analytic, Bz_analytic] = eval_B(xq, yq, zq, W, L, J);

%% Plot results

figure()
hold on;
quiver(yq, zq, By_analytic, Bz_analytic);
quiver(yq, zq, By, Bz);
hold off;
xlabel('y [m]', 'FontSize', 18);
ylabel('z [m]', 'FontSize', 18);
xlim([min(min(yq)) max(max(yq))]);
ylim([min(min(zq)) max(max(zq))]);
legend({'B_{analytic}', 'B_{computational}'}, 'FontSize', 16);

%% Long thin wire
W = 1e-7;
dy = 1e-8;

x = linspace(-L/2, L/2, 200);
y = linspace(-W/2, W/2, 200);
[x, y] = meshgrid(x, y);

Jx = J.*ones(size(x));
Jy = zeros(size(x));

[Bx, By, Bz] = calc_field(x, y, Jx, Jy, dx, dy, zeros(size(yq)), yq, zq);

r = sqrt(yq.^2 + zq.^2);

I = J*W;
B = (mu_0*I)./(2*pi.*r);

theta = atan2(zq, yq);
By_analytic = -B.*sin(theta);
Bz_analytic = B.*cos(theta);

figure();
hold on;
quiver(yq, zq, By, Bz);
quiver(yq, zq, By_analytic, Bz_analytic);
hold off;
xlabel('y [m]', 'FontSize', 18);
ylabel('z [m]', 'FontSize', 18);
xlim([min(min(yq)) max(max(yq))]);
ylim([min(min(zq)) max(max(zq))]);
legend({'B_{analytic}', 'B_{computational}'}, 'FontSize', 16);
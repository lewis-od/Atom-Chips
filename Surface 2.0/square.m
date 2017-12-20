%% Calculate the magnetic field of a 2D conducting sheet using the finite 
% element method

clear all;

%% Parameters and constants
sigma_0 = 5.96e7; % Conductivity of conductor [S*m^-1]
rho_0 = 1./sigma_0;
V0 = 1.0; % Voltage is V0 and -V0 at ends of conductor [V]
z = 0.5; % Distance from conductor [m]
mu_0 = 4e-7 * pi; % Permeability of free space

% Size of current elements to consider
dx = 0.05;
dy = 0.05;

%% Analytic solution for electric potential
res = 100; % Points per unit
x = linspace(-1, 1, 2*res);
y = linspace(-1, 1, 2*res);
a = (max(x) - min(x))/2; % Boudaries at (+/-)a

[x, y] = meshgrid(x, y);

phi = -V0/a .* x;

%% Calculate current density
[Ex, Ey] = gradient(phi);
Ex = -Ex;
Ey = -Ey;

% Gaussian spike of resisitivity at centre
rho = 1e6*exp(-(x.^2+y.^2)./0.005) + rho_0;
sigma = 1./rho;

% y component of current density should be ~0. Only consider Jx
Jx = sigma.*Ex; % Ohm's law

%% Apply finite element method
Nx = (max(x(1,:)) - min(x(1,:)))/dx;
Nx = floor(Nx);
Ny = (max(y(:,1)) - min(y(:,1)))/dy;
Ny = floor(Ny);

% Number of data points per current element
dims = size(x);
dxN = floor(dims(2)/Nx);
dyN = floor(dims(1)/Ny);

resolution = 200;
xq = linspace(-1, 1, resolution);
yq = linspace(-1, 1, resolution);
[xq, yq] = meshgrid(xq, yq);

By = zeros(resolution, resolution);
Bz = zeros(resolution, resolution);
% Loop over all current elements of conductor
for nx = 1:Nx
    for ny = 1:Ny
        % Centre of current element
        cx = (nx-1)*dx + dx/2;
        cx = cx - 1;
        cy = (ny-1)*dy + dy/2;
        cy = cy - 1;
    
        % Indices of current density matrix corresponding to this current
        % element
        xN_min = (nx-1)*dxN + 1;
        xN_max = nx*dxN;
        yN_min = (ny-1)*dyN + 1;
        yN_max = ny*dyN;
        
        % Treat the segment as having uniform current density
        Jn = Jx(yN_min:yN_max, xN_min:xN_max);
        Jn = mean(mean(Jn));
        
        % Calculate the field due to the current segment
        [dBy, dBz] = eval_B(xq-cx, yq-cy, z, dy, dx, Jn);
        By = By + dBy;
        Bz = Bz + dBz;
    end
end

%% Plot results
B = sqrt(By.^2 + Bz.^2);
subplot(1, 2, 1);
surf(xq, yq, B, 'EdgeColor', 'none');
view(2);
xlabel('x');
ylabel('y');
title('B Field');
colorbar();

subplot(1, 2, 2);
surf(x, y, rho, 'EdgeColor', 'none');
view(2);
xlabel('x');
ylabel('y');
title('Resistivity');
colorbar();

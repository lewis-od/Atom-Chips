%% Calculate the magnetic field of a 2D conducting sheet using the finite 
% element method

clear all;

%% Parameters and constants
sigma = 5.96e7; % Conductivity of conductor [S*m^-1]
V0 = 1.0; % Voltage is V0 and -V0 at ends of conductor [V]
z = 0.5; % Distance from conductor [m]
mu_0 = 4e-7 * pi; % Permeability of free space

dx = 0.1;
dy = 0.1;

%% Set up PDE solver
% Specify geometry
R1 = [3, 4, -1, 1, 1, -1, 1, 1, -1, -1]';
gm = [R1];
sf = 'R1';
ns = char('R1')';

model = createpde(1);
g = decsg(gm, sf, ns);
geometryFromEdges(model, g);

% Specify equation [(del)^2 phi = 0]
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);

% Specify boundary conditions
% Volatage of V0 on top edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 2, 'r', V0);
% Voltage of -V0 on bottom edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 4, 'r', -V0);
% Normal dervative must be 0 at all other boundaries
applyBoundaryCondition(model, 'neumann', 'Edge', [1 3], 'q', 0, 'g', 0);

% Generate mesh on which to solve PDE
generateMesh(model, 'Hmax', 0.1);

%% Solve PDE
result = solvepde(model);
res = 100; % Points per unit
x = linspace(-1, 1, 2*res);
y = linspace(-1, 1, 2*res);
[x, y] = meshgrid(x, y);

phi = interpolateSolution(result, x, y);
phi = reshape(phi, size(x));

%% Calculate current density
[Ex, Ey] = gradient(phi);
Ex = -Ex;
Ey = -Ey;

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

xq = linspace(-1, 1, 50);
yq = linspace(-1, 1, 50);

By = zeros(50, 50);
Bz = zeros(50, 50);
% Loop over each point in space
for xn = xq
    for yn = yq
        % Loop over all current elements of conductor
        for nx = 1:Nx
            for ny = 1:Ny
                % Centre of current element
                cx = nx*dx/2;
                cy = ny*dy/2;
                
                xN_min = (nx-1)*dxN + 1;
                xN_max = nx*dxN;
                yN_min = (ny-1)*dyN + 1;
                yN_max = ny*dyN;
                
                Jn = Jx(yN_min:yN_max, xN_min:xN_max);
                Jn = mean(mean(Jn));
                [dBy, dBz] = eval_B(xn-cx, yn-cy, z, dy, dx, Jn);
                By = By + dBy;
                Bz = Bz + dBz;
            end
        end
    end
end

B = sqrt(By.^2 + Bz.^2);
surf(xq, yq, B, 'EdgeColor', 'none');

%% Calculate the field for a z shaped wire on a 2DEG
% Based on the paper "Atom chips with free-standing two-dimensional 
% electron gases: advantages and challenges" by G. A. Sinuco-Leon, P.
% Kruger, and T.M. Fromhold
clear all;

%% Parameters and constants
V0 = 1; % Voltage difference of 2*V0 across wire
E0 = 0; % Electric field applied across 2DEG [V m^-1]
n = 3.3e15; % Mean electron density of 2DEG [m^-2]
mu = 140; % Mobility of 2DEG [m^2 V^-2 s^-1]
mu_0 = 4e-7 * pi; % Permeability of free space
B_bias = 0.081 * 3.4e-10;
z = 0.5; % z position to evaluate field at

% Size of current elements to consider
dx = 0.2;
dy = 0.2;

sigma = n*mu*1.6e-19; % Conductivity of 2DEG [S m^-1]

%% Specify geometry
R1 = [3 4 -3 0 0 -3 2.5 2.5 2.3 2.3]';
R2 = [3 4 -0.25 0.25 0.25 -0.25 2.5 2.5 -2.5 -2.5]';
R3 = [3 4 3 0 0 3 -2.5 -2.5 -2.3 -2.3]';
gm = [R1, R2, R3];
sf = 'R1+R2+R3';
ns = char('R1', 'R2', 'R3')';

model = createpde(1);
g = decsg(gm, sf, ns);
geometryFromEdges(model, g);

pdegplot(model, 'EdgeLabels', 'on');

% Equation is (del)^2 phi = 0 (Laplace's equation)
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);

%% Specify boundary conditions
% Volatage of V0 on top edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 1, 'r', V0);
% Voltage of -V0 on bottom edge
applyBoundaryCondition(model, 'dirichlet', 'Edge', 2, 'r', -V0);
% Normal dervative must be 0 at all other boundaries
applyBoundaryCondition(model, 'neumann', 'Edge', 3:17, 'q', 0, 'g', 0);

%% Solve equation
generateMesh(model, 'Hmax', 0.1);
result = solvepde(model);

%% Interpolate
res = 50;
x = linspace(-5, 5, 10*res);
y = linspace(-5, 5, 10*res);
[x, y] = meshgrid(x, y);

phi = interpolateSolution(result, x, y);
phi = reshape(phi, size(x));

surf(x, y, phi, 'EdgeColor', 'none');
view(2);
axis equal;
colormap cool;

%% Add potential fluctiations

% TODO: Add fluctuations calculated using Thomas-Fermi theory
% Assume delta_n(x,y) is white noise

%% Calculate current desnity
[Ex, Ey] = gradient(phi);
Ex = -Ex;
Ey = -Ey;

% TODO: remove NaNs before noise is added
Ex(isnan(Ex)) = 0.0;
Ey(isnan(Ey)) = 0.0;

Ex = Ex + E0;

Jx = sigma.*Ex;
Jy = sigma.*Ey;

%% Apply finite element method to calculate magnetic field
% Number of current elements in each dimension
Nx = (max(x(1,:)) - min(x(1,:)))/dx;
Nx = floor(Nx);
Ny = (max(y(:,1)) - min(y(:,1)))/dy;
Ny = floor(Ny);

% Number of data points per current element
dims = size(x);
dxN = floor(dims(2)/Nx);
dyN = floor(dims(1)/Ny);

resolution = 101;
xq = linspace(-5, 5, resolution);
yq = linspace(-5, 5, resolution);
[xq, yq] = meshgrid(xq, yq);

Bx = zeros(resolution, resolution);
By = zeros(resolution, resolution);
Bz = zeros(resolution, resolution);
% Loop over all current elements of conductor
for nx = 1:Nx
    for ny = 1:Ny
        % Centre of current element
        cx = (nx-1)*dx + dx/2;
        cx = cx - 5;
        cy = (ny-1)*dy + dy/2;
        cy = cy - 5;
    
        % Indices of current density matrix corresponding to this current
        % element
        xN_min = (nx-1)*dxN + 1;
        xN_max = nx*dxN;
        yN_min = (ny-1)*dyN + 1;
        yN_max = ny*dyN;
        
        % Treat the segment as having uniform current density, given by
        % the average of the density across the element
        Jn_x = Jx(yN_min:yN_max, xN_min:xN_max);
        Jn_y = Jy(yN_min:yN_max, xN_min:xN_max);
        
        Jn = sqrt(Jn_x.^2 + Jn_y.^2);
        Jn = mean(mean(Jn));
        
        if Jn == 0
            % No current => no magnetic field
            continue
        end
        
        % Calculate the field due to the current segment
        dBx = 0;
        dBy = 0;
        dBz = 0;
        if round(mean(Jx > Jy))
            % Current in x direction contributes to By
            [dBy, dBz] = eval_B(xq-cx, yq-cy, z, dy, dx, Jn);
        else
            % Current in y direction contributes to Bx
            [dBx, dBz] = eval_B(xq-cx, yq-cy, z, dy, dx, Jn);
        end
        
        Bx = Bx + dBx;
        By = By + dBy;
        Bz = Bz + dBz;
    end
end

Bx = Bx + ones(resolution, resolution).*B_bias;
B = sqrt(Bx.^2 + By.^2 + Bz.^2);

%% Plot results
surf(xq, yq, B, 'EdgeColor', 'none', 'FaceColor', 'interp');
xlabel('x');
ylabel('y');
colorbar();
colormap jet;

dims = size(xq);
x0 = ceil(dims(2)/2);
y0 = ceil(dims(1)/2);

figure();
hold on;
plot(xq(y0, :), B(y0, :));
plot(yq(:, x0), B(:, x0));
xlabel('x/y');
ylabel('|B|');
legend({'Varying x', 'Varying y'});

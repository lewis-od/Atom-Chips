% Calculate the B field of a 1d wire using the finite element method
% Inspired by https://uk.mathworks.com/matlabcentral/fileexchange/47368-3d-magnetic-field-computation-of-a-straight-wire-of-finite-length-using-biot-savart-s-law
clear all;

%% Define constants and parameters
mu_0 = pi*4e-7;
I = 355e-6;
L = 5e-6; % Actually L/2 - only named L for convinience

% How many points to evaluate the field at in each spatial dimension
Nx = 15;
Ny = 15;
Nz = 15;

% How many line elements to split the wire into
% Currently must be a multiple of 3 as we have 3 wires (150/3 = 50 elements
% per wire)
N = 150;

% x, y, and z coordinates of the wire
wx = zeros(1, N);
wy = zeros(1, N);
wz = zeros(1, N);

% First wire - parallel to y axis
wx(1:(N/3)) = -L;
wy(1:(N/3)) = linspace(0, 2*L, N/3);

% Second wire - parallel to x axis
wx((N/3)+1:(2*N/3)) = linspace(-L, L, N/3);
wy((N/3)+1:(2*N/3)) = 0.0;

% Third wire - parallel to y axis
wx((2*N/3)+1:N) = L;
wy((2*N/3)+1:N) = linspace(-2*L, 0, N/3);

%% Points in space where we will evaluate the field
x = linspace(-15e-6, 15e-6, Nx);
y = linspace(-15e-6, 15e-6, Ny);
z = linspace(0, 15e-6, Nz);
[X, Y, Z] = meshgrid(x, y, z);

% x, y, and z components of the total field at each point
Bx = zeros(Nx, Ny, Nz);
By = zeros(Nx, Ny, Nz);
Bz = zeros(Nx, Ny, Nz);

%% Calculate the magnetic field
% Loop through all points in space
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            % Loop through all line segments
            for n = 1:N-1
                % Vector of length dl pointing in direction of current flow
                dl = [wx(n+1)-wx(n) wy(n+1)-wy(n) wz(n+1)-wz(n)];
                % Mid-point of the line segment
                midpoint = 0.5 .* [wx(n)+wx(n+1) wy(n)+wy(n+1) wz(n)+wz(n+1)];
                % Position vector pointing from midpoint of line segment to
                % the point where we're calculating the field
                r = [X(i,j,k) Y(i,j,k) Z(i,j,k)] - midpoint;
                r_hat = (1/norm(r)) .* r;
                
                % Biot-Savart law
                dB = cross(dl, r_hat);
                dB = dB .* ((mu_0*I)/(4*pi*norm(r)^2));
                % Add the contribution from this line segment to the total
                % B-field at this point
                Bx(i,j,k) = Bx(i,j,k) + dB(1);
                By(i,j,k) = By(i,j,k) + dB(2);
                Bz(i,j,k) = Bz(i,j,k) + dB(3);
            end
        end
    end
end

%% Plot the results
figure(1);
hold on;
quiver3(X,Y,Z, Bx,By,Bz, 'AutoScaleFactor', 3.0);
xlabel('x');
ylabel('y');
zlabel('z');
line(wx(1:(N/3)), wy(1:(N/3)), wz(1:(N/3)), 'color', 'r', 'linewidth', 2);
line(wx((N/3)+1:(2*N/3)), wy((N/3)+1:(2*N/3)), wz((N/3)+1:(2*N/3)), 'color', 'r', 'linewidth', 2);
line(wx((2*N/3):N), wy((2*N/3):N), wz((2*N/3):N), 'color', 'r', 'linewidth', 2);

hold off;


% Calculate the B field of a 1d wire using the finite element method
% Inspired by https://uk.mathworks.com/matlabcentral/fileexchange/47368-3d-magnetic-field-computation-of-a-straight-wire-of-finite-length-using-biot-savart-s-law
clear all;

%% Define constants and parameters
mu_0 = pi*4e-7;
I = 355e-6;
L = 5e-6; % Actually L/2 - only named L for convinience

% How many points to evaluate the field at in each spatial dimension
Nx = 10;
Ny = 10;
Nz = 5;

% How many line elements to split the wire into
N = 50;

% x, y, and z coordinates of the wire
wx = linspace(-L, L, N);
wy = zeros(N);
wz = zeros(N);

%% Points in space where we will evaluate the field
x = linspace(-10e-6, 10e-6, Nx);
y = linspace(-10e-6, 10e-6, Ny);
z = linspace(0, 10e-6, Nz);
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
quiver3(X,Y,Z, Bx,By,Bz);
xlabel('x');
ylabel('y');
zlabel('z');
line([-L L], [0 0], [0 0], 'color', 'r', 'linewidth', 3);
hold off;


% Calculate the B field of a 1d wire using the finite element method
% Inspired by https://uk.mathworks.com/matlabcentral/fileexchange/47368-3d-magnetic-field-computation-of-a-straight-wire-of-finite-length-using-biot-savart-s-law
clear all;

%% Define constants and parameters
mu_0 = pi*4e-7;
% Parameters taken from [Haase et al]
I = 15; % Current through wire (A)
L1 = 6e-3; % Length of central wire (m)
L2 = 25e-3; % Length of end wires (m)
B_bias = 50e-4; % Strength of bias field in z direction (T)

% How many points to evaluate the field at in each spatial dimension
Nx = 30;
Nz = 30;

% Length of each line segment
dL = 0.1e-3;
N1 = floor(L1 / dL);
N2 = floor(L2 / dL);
N = N1 + 2*N2;

%% Setup wire positions
% x, y, and z coordinates of the wire
wx = zeros(1, N);
wy = zeros(1, N);
wz = zeros(1, N); % All wires in z=0 plane

% First end wire - parallel to y axis
wx(1:N2) = -(L1 / 2);
wy(1:N2) = linspace(0, L2, N2);

% Central wire - parallel to x axis
wx((N2+1):(N1+N2)) = linspace(-(L1/2), L1/2, N1);
wy((N2+1):(N1+N2)) = 0.0;

% Second end wire - parallel to y axis
wx((N1+N2+1):N) = (L1 / 2);
wy((N1+N2+1):N) = linspace(-L2, 0, N2);

%% Setup arrays for axes and the field
% Points in space where we will evaluate the field
x = linspace(-(L1), L1, Nx);
z = linspace(0, 4e-3, Nz);
[X, Z] = meshgrid(x, z);
Y = zeros(Nx, Nz);

% x, y, and z components of the total field at each point
Bx = zeros(Nx, Nz);
By = zeros(Nx, Nz);
Bz = zeros(Nx, Nz);

%% Calculate the magnetic field
% Loop through all points in space
for i = 1:Nx
    for k = 1:Nz
        % Loop through all line segments
        for n = 1:N-1
            % Vector of length dl pointing in direction of current flow
            dl = [wx(n+1)-wx(n) wy(n+1)-wy(n) wz(n+1)-wz(n)];
            % Mid-point of the line segment
            midpoint = 0.5 .* [wx(n)+wx(n+1) wy(n)+wy(n+1) wz(n)+wz(n+1)];
            % Position vector pointing from midpoint of line segment to
            % the point where we're calculating the field
            r = [X(i,k) 0.0 Z(i,k)] - midpoint;
            r_hat = (1/norm(r)) .* r; % Unit vector used in cross product

            % Biot-Savart law
            dB = cross(dl, r_hat);
            dB = dB .* ((mu_0*I)/(4*pi*norm(r)^2));
            % Add the contribution from this line segment to the total
            % B-field at this point
            Bx(i,k) = Bx(i,k) + dB(1);
            By(i,k) = By(i,k) + dB(2);
            Bz(i,k) = Bz(i,k) + dB(3);
        end
    end
end

% Add bias field
Bz = Bz + B_bias;

%% Plot the results
figure();
hold on;
% Plot the magnetic field
quiver3(X,Y,Z, Bx,By,Bz, 'AutoScaleFactor', 5);
xlabel('x');
ylabel('y');
zlabel('z');
% Show the positions of the wires
line(wx(1:N2), wy(1:N2), wz(1:N2), 'color', 'r', 'linewidth', 2);
line(wx((N2+1):(N1+N2)), wy((N2+1):(N1+N2)), wz((N2+1):(N1+N2)), 'color', 'r', 'linewidth', 2);
line(wx((N1+N2+1):N), wy((N1+N2+1):N), wz((N1+N2+1):N), 'color', 'r', 'linewidth', 2);
hold off;

%% Plot individual components of B
% figure();
% surf(X,Z,Bx);
% xlabel('x');
% ylabel('z');
% zlabel('B_x');
% 
% figure();
% surf(X,Z,By);
% xlabel('x');
% ylabel('z');
% zlabel('B_y');
% 
% figure();
% surf(X,Z,Bz);
% xlabel('x');
% ylabel('z');
% zlabel('B_z');
B = sqrt(Bx.^2 + By.^2 + Bx.^2);
figure();
surf(X,Z,B);
xlabel('x');
ylabel('z');
zlabel('B');

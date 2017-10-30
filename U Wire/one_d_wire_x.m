% Calculate the B field of a 1d wire using the finite element method
clear all;

%% Define constants and parameters
mu_0 = pi*4e-7;
% Parameters taken from [Reichel et al]
I = 2; % Current through wire (A)
L1 = 250e-6; % Length of central wire (m)
L2 = L1 * 5; % Length of end wires (m)
B_bias = 162e-4; % Strength of bias field (T)
width = 60e-6; % Width of wire

% How many points to evaluate the field at in each spatial dimension
Nx = 600;

% Length of each line segment
dL = 1e-6;
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
wy((N1+N2+1):N) = linspace(0, L2, N2);

%% Setup arrays for axes and the field
% Points in space where we will evaluate the field
X = linspace(-300e-6, 300e-6, Nx);
Y = zeros(1, Nx);
Z = ones(1, Nx) .* 20e-6;

% x, y, and z components of the total field at each point
Bx = zeros(1, Nx);
By = zeros(1, Nx);
Bz = zeros(1, Nx);

%% Calculate the magnetic field
% Loop through all points in space
for i = 1:Nx
    % Loop through all line segments
    for n = 1:N-1
        % Vector of length dl pointing in direction of current flow
        dl = [wx(n+1)-wx(n) wy(n+1)-wy(n) wz(n+1)-wz(n)];
        % Mid-point of the line segment
        midpoint = 0.5 .* [wx(n)+wx(n+1) wy(n)+wy(n+1) wz(n)+wz(n+1)];
        % Position vector pointing from midpoint of line segment to
        % the point where we're calculating the field
        r = [X(i) Y(i) Z(i)] - midpoint;
        r_hat = (1/norm(r)) .* r; % Unit vector used in cross product

        % Biot-Savart law
        dB = cross(dl, r_hat);
        dB = dB .* ((mu_0*I)/(4*pi*norm(r)^2));
        % Add the contribution from this line segment to the total
        % B-field at this point
        Bx(i) = Bx(i) + dB(1);
        By(i) = By(i) + dB(2);
        Bz(i) = Bz(i) + dB(3);
    end
end

Bz = Bz + B_bias;

%% Plot |B|
figure();
hold on;
plot(X, Bz);
plot(X, By);
plot(X, Bx);
xlabel('x (m)');
ylabel('B (T)');
title('Z-wire at y=0 and z=20\mu m');
line(wx((N2+1):(N1+N2)), zeros(1, N1), zeros(1, N1), 'color', 'r', 'linewidth', 1);
legend({'Bx', 'By', 'Bz', 'Wire'});
hold off;
figure();
plot(X, sqrt(Bx.^2 + By.^2 + Bz.^2));
line(wx((N2+1):(N1+N2)), zeros(1, N1), zeros(1, N1), 'color', 'r', 'linewidth', 1);
xlabel('x (m)');
ylabel('|B| (T)');
title('Z-wire at y=0 and z=20\mu m');

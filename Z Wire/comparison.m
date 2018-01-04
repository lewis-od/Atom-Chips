% Calculate the B field of a 1d wire using the finite element method
% Inspired by https://uk.mathworks.com/matlabcentral/fileexchange/47368-3d-magnetic-field-computation-of-a-straight-wire-of-finite-length-using-biot-savart-s-law
clear all;

%% Define constants and parameters
mu_0 = pi*4e-7;
% Parameters taken from [Haase et al]
I = 15; % Current through wire (A)
B_bias = 20e-4; % Strength of bias field in z direction (T)

% How many points to evaluate the field at in each spatial dimension
Nx = 20;
Nz = 20;

% Length of each line segment
dL = 0.01;
N = floor(10/dL);

%% Setup wire positions
% x, y, and z coordinates of the wire
% "Infinitely" long wire along y axis
wx = zeros(1, N);
wy = linspace(5, -5, N);
wz = zeros(1, N);

%% Setup arrays for axes and the field
% Points in space where we will evaluate the field
x = linspace(-5e-3, 5e-3, Nx);
Y = 0;
z = linspace(-5e-3, 5e-3, Nz);
[X, Z] = meshgrid(x, z);

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
            r = [X(i,k) Y Z(i,k)] - midpoint;
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

%% Calculate field from analytic expression

r = sqrt(X.^2 + Y.^2+ Z.^2);
B = (mu_0*I)./(2*pi.*r);

theta = atan2(Z, X);
Bx_analytic = -B.*sin(theta);
Bz_analytic = B.*cos(theta);

%% Plot results

X = X./1e-3;
Z = Z./1e-3;

scale = 1.7;
figure('Position', [680, 558, 560*2*scale 420*scale]);
subplot(1,2,1);
hold on;
quiver(X, Z, Bx_analytic, Bz_analytic);
quiver(X, Z, Bx, Bz);
hold off;
xlabel('x [mm]', 'FontSize', 18);
ylabel('z [mm]', 'FontSize', 18);
legend({'B_{analytic}', 'B_{computational}'}, 'FontSize', 16);
xlim([min(min(X)) max(max(X))]);
ylim([min(min(Z)) max(max(Z))]);

subplot(1,2,2);
hold on;
quiver(X, Z, Bx_analytic, Bz_analytic);
quiver(X, Z, Bx, Bz);
hold off;
xlabel('x [mm]', 'FontSize', 18);
ylabel('z [mm]', 'FontSize', 18);
legend({'B_{analytic}', 'B_{computational}'}, 'FontSize', 16);
xlim([min(min(X)) max(max(X))]./2);
ylim([min(min(Z)) max(max(Z))]./2);

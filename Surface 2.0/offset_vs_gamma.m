clear all;

V0 = 1.6e-3 * 640; % Voltage difference of V0 across wire
n = 3.3e15; % Mean electron density of 2DEG [m^-2]
mu = 140; % Mobility of 2DEG [m^2 V^-2 s^-1]
add_noise = false;

sigma = n*mu*1.6e-19; % Conductivity of 2DEG [S m^-1]

%% Calculate current density
res = 5e6; % Points per unit for interpolation
[x, y, phi] = calc_potential(V0, res);

h = 1/res;
[Ex, Ey] = gradient(phi, h);
Ex = -Ex;
Ey = -Ey;

% Remove NaNs
Ex(isnan(Ex)) = 0.0;
Ey(isnan(Ey)) = 0.0;

if add_noise
    phi_noise = calc_noise(x, n, 80e-9);
    [Ex_noise, Ey_noise] = gradient(phi_noise, h)
    Ex_noise = -Ex_noise;
    Ey_noise = -Ey_noise;
    
    Ex = Ex + Ex_noise;
    Ey = Ey + Ey_noise;
end

Jx = sigma.*Ex;
Jy = sigma.*Ey;

%% Calculate Bs and z0

offset_factors = linspace(0.0, 0.5, 10);
bias_factor = 0.9;

% Size of current elements [m]
dx = 2e-6;
dy = 2e-6;

% Field at surface of conductor
[Bx0, By0, Bz0] = calc_field(x, y, Jx, Jy, dx, dy, 0, 0, 0.5e-12);
Bs = sqrt(Bx0.^2 + By0.^2 + Bz0.^2);

B_bias = bias_factor*Bs;

% Calculate B for r=(0, 0, z
zq = linspace(0, 15e-6, 301);
[Bx3, By3, Bz3] = calc_field(x, y, Jx, Jy, dx, dy, 0, 0, zq);
Bx3 = Bx3 + B_bias; % No need to add offset - doesn't affect z0
B3 = sqrt(Bx3.^2 + By3.^2 + Bz3.^2);

% Find z0
[o, z0_ind] = min(B3);
z0 = zq(z0_ind);

%% Calculate omega and gamma for varying B_offset
N = length(offset_factors);
gammas = zeros(1, N);
omegas = zeros(1, N);
omega_hats = zeros(1, N);
z0s = zeros(1, N);
r2s = zeros(3, N);
for i = 1:length(offset_factors)
   offset_factor = offset_factors(i);
   B_offset = offset_factor * Bs;
   [omega, omega_hat, gamma] = calc_frequency(x, y, Jx, Jy, B_bias, B_offset, z0);
   
   gammas(i) = gamma;
   omegas(i) = omega;
   omega_hats(i) = omega_hat;
   z0s(i) = z0;
end

%% Plot results

figure();
semilogy(offset_factors, gammas);
xlabel('B_{offset}/B_s', 'FontSize', 18);
ylabel('\Gamma_{MF} [s^{-1}]', 'FontSize', 18);

figure();
hold on;
plot(offset_factors, omegas);
plot(offset_factors, omega_hats);
xlabel('B_{offset}/B_s', 'FontSize', 18);
ylabel('\omega [Hz]', 'FontSize', 18);
legend({'$\omega$', '$\hat{\omega}$'}, 'FontSize', 16, 'Interpreter', 'Latex');

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

%% Calculate trap frequencies, trap distances, and loss-rates

offset_factors = linspace(0, 0.2, 10);
bias_factor = 0.9;

N = length(offset_factors);
gammas = zeros(1, N);
omegas = zeros(1, N);
z0s = zeros(1, N);
r2s = zeros(3, N);
for i = 1:length(offset_factors)
   offset_factor = offset_factors(i);
   [omega, omega_hat, gamma, z0, Bs] = calc_frequency(x, y, Jx, Jy, bias_factor, offset_factor);
   
   gammas(i) = gamma;
   omegas(i) = omega;
   z0s(i) = z0;
end

%%

figure();
semilogy(offset_factors, gammas);
xlabel('B_{offset}/B_s', 'FontSize', 18);
ylabel('\Gamma_{MF} [s^{-1}]', 'FontSize', 18);
clear all;

V0 = 1.6e-2 * 640; % Voltage difference of V0 across wire
n = 3.3e15; % Mean electron density of 2DEG [m^-2]
mu = 140; % Mobility of 2DEG [m^2 V^-2 s^-1]
add_noise = false;

sigma = n*mu*1.6e-19; % Conductivity of 2DEG [S m^-1]
sigma = sigma.*1e6;

%% Calculate current density
res = 5e6; % Points per unit for interpolation
[x, y, phi] = calc_potential(V0, res);

[Ex, Ey] = gradient(phi);
Ex = -Ex;
Ey = -Ey;

% Remove NaNs
Ex(isnan(Ex)) = 0.0;
Ey(isnan(Ey)) = 0.0;

if add_noise
    phi_noise = calc_noise(x, n, 80e-9);
    Ex_noise = -Ex_noise;
    Ey_noise = -Ey_noise;
    
    Ex = Ex + Ex_noise;
    Ey = Ey + Ey_noise;
end

Jx = sigma.*Ex;
Jy = sigma.*Ey;

%%
offset_factors = linspace(0, 0.5, 20);
bias_factor = 0.9;

N = length(offset_factors);
gammas = zeros(1, N);
gamma_errors = zeros(1, N);
omegas = zeros(1, N);
omega_errors = zeros(1, N);
r2s = zeros(3, N);
for i = 1:length(offset_factors)
   offset_factor = offset_factors(i);
   [omega, omega_err, r_squared, omega_hat, gamma, gamma_err] = calc_frequency(x, y, Jx, Jy, bias_factor, offset_factor);
   
   gammas(i) = gamma;
   gamma_errors(i) = gamma_err;
   omegas(i) = norm(omega);
   omega_error = sqrt(sum(omega_err.^2 ./ omega.^2))/norm(omega);
   omega_errors(i) = omega_error;
   r2s(:,i) = r_squared;
   
%    disp(r_squared);
end

%%

figure();
semilogy(offset_factors, gammas);
% ax = axes;
% errorbar(offset_factors, gammas, gamma_errors);
% ax.YScale = 'log';
xlabel('B_{offset}/B_s', 'FontSize', 18);
ylabel('\Gamma_{MF} [s^{-1}]', 'FontSize', 18);
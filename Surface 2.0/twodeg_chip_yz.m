%% Calculate the field for a z shaped wire on a 2DEG
% Based on the paper "Atom chips with free-standing two-dimensional 
% electron gases: advantages and challenges" by G. A. Sinuco-Leon, P.
% Kruger, and T.M. Fromhold
clear all;

%% Parameters and constants
V0 = 1.6e-3 * 640; % Voltage difference of V0 across wire
n = 3.3e15; % Mean electron density of 2DEG [m^-2]
mu = 140; % Mobility of 2DEG [m^2 V^-2 s^-1]
mu_0 = 4e-7 * pi; % Permeability of free space
B_bias_factor = 0.9; % B_bias = B_bias_factor * Bs
B_offset_factor = 0.2; % B_offset = B_offset_factor * Bs

x_fixed = 0; % x position to evaluate field at [m]

add_noise = false;

% Size of current elements to consider [m]
dx = 2e-6;
dy = 2e-6;

sigma = n*mu*1.6e-19; % Conductivity of 2DEG [S m^-1]

%% Calculate and plot potential
res = 5e6; % Points per unit for interpolation
[x, y, phi] = calc_potential(V0, res);

% Plot potential
% surf(x, y, phi, 'EdgeColor', 'none');
% view(2);
% xlim([min(min(x)), max(max(x))]);
% ylim([min(min(y)), max(max(y))]);
% colormap cool;

%% Calculate current desnity
h = 1/res;

[Ex, Ey] = gradient(phi, h);
Ex = -Ex;
Ey = -Ey;

% Remove NaNs
Ex(isnan(Ex)) = 0.0;
Ey(isnan(Ey)) = 0.0;

%% Add potential fluctiations
if add_noise
    phi_noise = calc_noise(x, n, 80e-9);
    [Ex_noise, Ey_noise] = gradient(phi_noise, h);
    Ex_noise = -Ex_noise;
    Ey_noise = -Ey_noise;
    
    Ex = Ex + Ex_noise;
    Ey = Ey + Ey_noise;
end

Jx = sigma.*Ex;
Jy = sigma.*Ey;

%% Apply finite element method to calculate magnetic field
% Points to calculate B at
resolution = 501;
% yq = linspace(-150e-6, 150e-6, resolution);
yq = 0;
zq = linspace(0, 30e-6, resolution);
% [yq, zq] = meshgrid(yq, zq);

[Bx, By, Bz] = calc_field(x, y, Jx, Jy, dx, dy, x_fixed, yq, zq);

[Bx0, By0, Bz0] = calc_field(x, y, Jx, Jy, dx, dy, 0, 0, 1e-11);
Bs = sqrt(Bx0.^2 + By0.^2 + Bz0.^2);

B_bias = B_bias_factor*Bs;
B_offset = B_offset_factor*Bs;

%% Add bias/offset and plot results

By = By + B_offset;
Bx = Bx + B_bias;

B = sqrt(Bx.^2 + By.^2 + Bz.^2);
% Plot results
% figure();
% surf(yq, zq, B, 'EdgeColor', 'none', 'FaceColor', 'interp');
% xlabel('y [m]', 'FontSize', 18);
% ylabel('z [m]', 'FontSize', 18);
% c = colorbar();
% c.Label.String = '|B| [T]';
% c.Label.FontSize = 18;
% colormap jet;
% view(2);
% 
% dims = size(yq);
% y0 = ceil(dims(2)/2);
% 
% figure();
plot(zq/1e-6, B);
xlabel('z [\mum]', 'FontSize', 20);
ylabel('|{\bfB}| [T]', 'FontSize', 20);

%% Plot how z0 varies with B_bias

B_bias_factor = linspace(0.5, 0.95, 5);
figure();
z0 = zeros(1, length(B_bias_factor));
subplot(1, 2, 1);
hold on;
for i = 1:length(B_bias_factor)
    Bx_plot = Bx + B_bias_factor(i)*Bs;
    B = sqrt(Bx_plot.^2 + By.^2 + Bz.^2);
    
    [t, ind] = min(B);
    z0(i) = zq(ind);

    plot(zq, B);
end
xlabel('z', 'FontSize', 18);
ylabel('|B|', 'FontSize', 18);
legend(split(num2str(B_bias_factor)).');

subplot(1, 2, 2);
plot(B_bias_factor, z0);
xlabel('B_{bias}/B_s', 'FontSize', 18);
ylabel('z_0 [m]', 'FontSize', 18);

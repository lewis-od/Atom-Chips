%% Calculate the field for a z shaped wire on a 2DEG
% Based on the paper "Atom chips with free-standing two-dimensional 
% electron gases: advantages and challenges" by G. A. Sinuco-Leon, P.
% Kruger, and T.M. Fromhold
clear all;

%% Parameters and constants
V0 = (5e5/9.743508582048404)*200; % Voltage difference of V0 across wire
E0 = 0; % Electric field applied across 2DEG [V m^-1]
n = 3.3e15; % Mean electron density of 2DEG [m^-2]
mu = 140; % Mobility of 2DEG [m^2 V^-2 s^-1]
mu_0 = 4e-7 * pi; % Permeability of free space
B_bias = 1.54e-4; % Bias field strength [T]
B_offset = 0.2e-4; % Offset field strength [T]

x_fixed = 0; % z position to evaluate field at [um]

add_noise = false;

% Size of current elements to consider [um]
dx = 2;
dy = 2;

sigma = n*mu*1.6e-19; % Conductivity of 2DEG [S m^-1]

%% Calculate and plot potential
res = 5; % Points per unit for interpolation
[x, y, phi] = calc_potential(V0, res);

% Scale down all variables
x = x.*1e-6;
y = y.*1e-6;
x_fixed = x_fixed.*1e-6;
dx = dx.*1e-6;
dy = dy.*1e-6;
res = res.*1e-6;

% Plot potential
% surf(x, y, phi, 'EdgeColor', 'none');
% view(2);
% xlim([min(min(x)), max(max(x))]);
% ylim([min(min(y)), max(max(y))]);
% colormap cool;

%% Calculate current desnity
[Ex, Ey] = gradient(phi);
Ex = -Ex;
Ey = -Ey;

% Remove NaNs
Ex(isnan(Ex)) = 0.0;
Ey(isnan(Ey)) = 0.0;

%% Add potential fluctiations
if add_noise
    q = 1.6e-19; % Charge on an electron [C]
    epsilon_0 = 8.85e-12; % Permittivity of free space [F m^-1]
    epsilon = 12.9; % Relative permittivity of GaAs
    k_s = 2.1e8; % Screening wave vector [m^-1]
    d = 52.9e-9; % Distance of 2DEG from donor layer [m]
    grain_size = 80e-9; % Typical grain size of noise [m]

    scale_factor = grain_size/(max(max(x)) - min(min(x)));
    scale_factor = scale_factor * res * 1e6;

    % Add white noise to donor density distribution
    n_dist = ones(size(x)/scale_factor).*n;
    noise = (imnoise(n_dist, 'gaussian', 0, 0.1) - 0.5)*2;
    n_dist = noise.*n + n_dist;
    n_dist = imresize(n_dist, scale_factor, 'nearest');

    delta_n = ones(size(x)).*n - n_dist;

    % Apply Thomas-Fermi screening model
    delta_n_hat = fft2(delta_n);
    delta_n_hat = fftshift(delta_n_hat); % Shift k=0 to centre

    dims = size(x);
    kx = (0:1:(dims(1)-1)).*(dims(1)./res);
    ky = (0:1:(dims(2)-1)).*(dims(2)./res);
    % Ensure k=0 is at centre
    kx = kx - (max(kx)/2);
    ky = ky - (max(ky)/2);
    [kx, ky] = meshgrid(kx, ky);
    k = sqrt(kx.^2 + ky.^2);

    % Integrand of Thomas-Fermi screening equation
    phi_noise_hat = (exp(-k.*d) .* delta_n_hat)./(k + k_s);

    phi_noise = ifft2(phi_noise_hat);
    phi_noise = abs(phi_noise);

    % Scale appropriately
    phi_noise = phi_noise .* (q^2)/(4*pi*epsilon*epsilon_0); % In [J]
    phi_noise = phi_noise ./ q; % In [eV]
    phi_noise = phi_noise ./ (1e-3); % In [meV]

    figure();
    subplot(1, 2, 1);
    surf(x, y, n_dist, 'EdgeColor', 'none');
    colorbar();
    view(2);
    title('Ionised Donor Desnity');

    subplot(1, 2, 2);
    surf(x, y, phi_noise, 'EdgeColor', 'none');
    colorbar();
    view(2);
    title('Screened Potential');
    
    [Ex_noise, Ey_noise] = gradient(phi_noise);
    Ex_noise = -Ex_noise;
    Ey_noise = -Ey_noise;
    
    Ex = Ex + Ex_noise;
    Ey = Ey + Ey_noise;
end

Ex = Ex + E0;

Jx = sigma.*Ex;
Jy = sigma.*Ey;

%% Apply finite element method to calculate magnetic field
% Points to calculate B at
resolution = 101;
yq = linspace(-150e-6, 150e-6, resolution);
zq = linspace(0, 20e-6, resolution);
[yq, zq] = meshgrid(yq, zq);

[Bx, By, Bz] = calc_field(x, y, Jx, Jy, dx, dy, x_fixed, yq, zq);

%% Add bias/offset and plot results

By = By + B_offset;

B = sqrt(Bx.^2 + By.^2 + Bz.^2);
% Plot results
figure();
surf(yq, zq, B, 'EdgeColor', 'none', 'FaceColor', 'interp');
xlabel('y', 'FontSize', 18);
ylabel('z', 'FontSize', 18);
c = colorbar();
c.Label.String = '|B|';
c.Label.FontSize = 18;
colormap jet;
view(2);

dims = size(yq);
y0 = ceil(dims(2)/2);

B_bias = (5:1:10).*1e-5;
z0 = zeros(1, length(B_bias));
figure();
subplot(1, 2, 1);
hold on;
for i = 1:length(B_bias)
    Bx_plot = Bx + B_bias(i);
    B = sqrt(Bx_plot.^2 + By.^2 + Bz.^2);
    
    [t, ind] = min(B(:, y0));
    z0(i) = zq(ind, y0);

    plot(zq(:, y0), B(:, y0));
end
xlabel('z');
ylabel('|B|');
legend(split(num2str(B_bias)).');

subplot(1, 2, 2);
plot(B_bias, z0);
xlabel('B_{bias}');
ylabel('z_0');

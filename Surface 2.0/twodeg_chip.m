%% Calculate the field for a z shaped wire on a 2DEG
% Based on the paper "Atom chips with free-standing two-dimensional 
% electron gases: advantages and challenges" by G. A. Sinuco-Leon, P.
% Kruger, and T.M. Fromhold
clear all;

%% Parameters and constants
V0 = 2; % Voltage difference of 2*V0 across wire
E0 = 0; % Electric field applied across 2DEG [V m^-1]
n = 3.3e15; % Mean electron density of 2DEG [m^-2]
mu = 140; % Mobility of 2DEG [m^2 V^-2 s^-1]
mu_0 = 4e-7 * pi; % Permeability of free space
B_bias = 0.081 * 3.4e-10 * 5; % Bias field strength [T]
B_offset = 0.081 * 3.4e-11; % Offset field strength [T]
z = 0.5; % z position to evaluate field at [um]

add_noise = true;

% Size of current elements to consider
dx = 0.05;
dy = 0.05;

sigma = n*mu*1.6e-19; % Conductivity of 2DEG [S m^-1]

%% Calculate and plot potential
res = 100; % Points per unit for interpolation
[x, y, phi] = calc_potential(V0, res);

% Scale down all variables
x = x.*1e-6;
y = y.*1e-6;
z = z.*1e-6;
dx = dx.*1e-6;
dy = dy.*1e-6;
res = res.*1e-6;

% Plot potential
surf(x, y, phi, 'EdgeColor', 'none');
view(2);
xlim([min(min(x)), max(max(x))]);
ylim([min(min(y)), max(max(y))]);
colormap cool;

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
xq = linspace(-5e-6, 5e-6, resolution);
yq = linspace(-5e-6, 5e-6, resolution);
[xq, yq] = meshgrid(xq, yq);

[Bx, By, Bz] = calc_field(x, y, z, Jx, Jy, dx, dy, xq, yq);

%% Add bias/offset and plot results

Bx = Bx + B_bias;
By = By + B_offset;

B = sqrt(Bx.^2 + By.^2 + Bz.^2);

% Plot results
figure();
surf(xq, yq, B, 'EdgeColor', 'none', 'FaceColor', 'interp');
xlabel('x', 'FontSize', 18);
ylabel('y', 'FontSize', 18);
c = colorbar();
c.Label.String = '|B|';
c.Label.FontSize = 18;
colormap jet;
view(2);

dims = size(xq);
x0 = ceil(dims(2)/2);
y0 = ceil(dims(1)/2);

figure();
hold on;
plot(xq(y0, :), B(y0, :));
plot(yq(:, x0), B(:, x0));
xlabel('x', 'FontSize', 18);
ylabel('|B|', 'FontSize', 18);
legend({'B(x,0)', 'B(0,y)'}, 'FontSize', 16);

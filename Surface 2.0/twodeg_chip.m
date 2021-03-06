%% Calculate the field for a z shaped wire on a 2DEG
% Based on the paper "Atom chips with free-standing two-dimensional 
% electron gases: advantages and challenges" by G. A. Sinuco-Leon, P.
% Kruger, and T.M. Fromhold
clear all;

%% Parameters and constants
V0 = 1.6e-3 * 640; % 1.6 mV/um
n = 3.3e15; % Mean electron density of 2DEG [m^-2]
mu = 140; % Mobility of 2DEG [m^2 V^-2 s^-1]
mu_0 = 4e-7 * pi; % Permeability of free space
B_bias_factor = 0.9; % B_bias = B_bias_factor * Bs
B_offset_factor = 0.2; % B_offset = B_offset_factor * Bs

z = 1.56e-06; % z position to evaluate field at [m]

add_noise = false;

% Size of current elements to consider [um]
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
%     phi_noise = phi_noise ./ (1.6e-19); % In [eV]
%     phi_noise = phi_noise ./ (1e-3); % In [meV]
%     surf(x, y, phi_noise, 'EdgeColor', 'none');
%     view(2);
%     colorbar();
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
resolution = 161;
xq = linspace(-150e-6, 150e-6, resolution);
yq = linspace(-150e-6, 150e-6, resolution);
[xq, yq] = meshgrid(xq, yq);

[Bx, By, Bz] = calc_field(x, y, Jx, Jy, dx, dy, xq, yq, z);

[Bx0, By0, Bz0] = calc_field(x, y, Jx, Jy, dx, dy, 0, 0, 1e-11);
Bs = sqrt(Bx0.^2 + By0.^2 + Bz0.^2);

B_bias = B_bias_factor*Bs;
B_offset = B_offset_factor*Bs;

%% Add bias/offset and plot results
Bx = Bx + B_bias;
By = By + B_offset;

B = sqrt(Bx.^2 + By.^2 + Bz.^2);

%% Plot results
figure();
surf(xq/1e-6, yq/1e-6, B, 'EdgeColor', 'none', 'FaceColor', 'interp');
xlabel('x [\mum]', 'FontSize', 20);
ylabel('y [\mum]', 'FontSize', 20);
c = colorbar();
c.Label.String = '|{\bfB}| [T]';
c.Label.FontSize = 20;
colormap jet;
view(2);

dims = size(xq);
x0 = ceil(dims(2)/2);
y0 = ceil(dims(1)/2);

figure();
hold on;
plot(xq(y0, :)/1e-6, B(y0, :));
plot(yq(:, x0)/1e-6, B(:, x0));
xlabel('x/y [\mum]', 'FontSize', 20);
ylabel('|{\bfB}| [T]', 'FontSize', 20);
legend({'|{\bfB}(x,0)|', '|{\bfB}(0,y)|'}, 'FontSize', 20);

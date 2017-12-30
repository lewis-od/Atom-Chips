%% Calculate the trap frequency for a z shaped wire on a 2DEG
% Based on the paper "Atom chips with free-standing two-dimensional 
% electron gases: advantages and challenges" by G. A. Sinuco-Leon, P.
% Kruger, and T.M. Fromhold
clear all;

%% Parameters and constants
V0 = 1.6e-3 * 640; % Voltage difference of V0 across wire
n = 3.3e15; % Mean electron density of 2DEG [m^-2]
mu = 140; % Mobility of 2DEG [m^2 V^-2 s^-1]
mu_0 = 4e-7 * pi; % Permeability of free space
B_bias = 10.5e-6; % Bias field strength [T]
B_offset = 1e-6; % Offset field strength [T]

x_fixed = 0; % z position to evaluate field at [m]

add_noise = false;

% Size of current elements to consider [m]
dx = 2e-6;
dy = 2e-6;

sigma = n*mu*1.6e-19; % Conductivity of 2DEG [S m^-1]
sigma = sigma.*1e6;

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
[Ex, Ey] = gradient(phi);
Ex = -Ex;
Ey = -Ey;

% Remove NaNs
Ex(isnan(Ex)) = 0.0;
Ey(isnan(Ey)) = 0.0;

%% Add potential fluctiations
if add_noise
    phi_noise = calc_noise(x, n, 80e-9);
    Ex_noise = -Ex_noise;
    Ey_noise = -Ey_noise;
    
    Ex = Ex + Ex_noise;
    Ey = Ey + Ey_noise;
end

Jx = sigma.*Ex;
Jy = sigma.*Ey;

%% Apply finite element method to calculate magnetic field
% Points to calculate B at
resolution = 201;
yq = linspace(-150e-6, 150e-6, resolution);
% yq = 0;
zq = linspace(0, 30e-6, resolution);
[yq, zq] = meshgrid(yq, zq);

[Bx, By, Bz] = calc_field(x, y, Jx, Jy, dx, dy, x_fixed, yq, zq);

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

dims = size(yq);
y0 = ceil(dims(2)/2);

figure();
plot(zq(:, y0), B(:, y0));
xlabel('z', 'FontSize', 18);
ylabel('|B|', 'FontSize', 18);

B_along_z = B(:, y0);
zs = zq(:, y0);

[o, z0_ind] = min(B_along_z);
z0 = zs(z0_ind);

B_along_y = B(z0_ind-1, :);
ys = yq(z0_ind-1, :);

%% Calculate trap frequency
% Trap freq. for an Rubidium-87 atom in the |F=2, m_F=2> state
mF = 2;
gF = 1/2; % Lande g-factor (theoretical value)
mu_B = 9.274e-24; % Bohr magneton [J/T]
m = 87 * 1.66054e-27; % Mass of atoms [kg]

Vz = mu_B .* B_along_z;
z = zq(:, y0);
[z, Vz] = prepareCurveData(z, Vz);

Vz_fit = @(omega, x) (m*omega.^2)/2 * (x-z0).^2 + mF*gF*mu_B*By(z0_ind, y0);
ft = fittype(Vz_fit);
fo = fitoptions('Method', 'NonLinearLeastSquares');
% fo.Exclude = z > 3*z0;
fo.Exclude = z > 2*z0;
fo.StartPoint = 5e3;
fo.TolFun = 1e-80;
fo.MaxIter = 8000;
fo.MaxFunEvals = 10000;
fo.Lower = -Inf;
fo.Upper = Inf;
fo.DiffMaxChange = 1e3;
fo.DiffMinChange = 1e-5;
res_z = fit(z, Vz, ft, fo);
figure();
subplot(1,2,1);
plot(res_z, z, Vz);
xlim([0 2*z0]);

omega_z = res_z.omega;
bounds = confint(res_z);
omega_z_err = (bounds(2) - bounds(1))/2;

Vy = mu_B .* B_along_y;
y = yq(z0_ind, :);

[o, ind] = max(Vy);
max_y = y(ind);
y_lim = -0.5*max_y;

[y, Vy] = prepareCurveData(y, Vy);

Vy_fit = @(omega, x) (m*omega.^2)/2 * x.^2 + mF*gF*mu_B*By(z0_ind, y0);
ft = fittype(Vy_fit);
fo = fitoptions('Method', 'NonLinearLeastSquares');
fo.Exclude = or(y < -y_lim, y > y_lim);
fo.StartPoint = 5e3;
fo.TolFun = 1e-80;
fo.MaxIter = 8000;
fo.MaxFunEvals = 10000;
fo.Lower = -Inf;
fo.Upper = Inf;
fo.DiffMaxChange = 1e3;
fo.DiffMinChange = 1e-5;
res_y = fit(y, Vy, ft, fo);
subplot(1, 2, 2);
plot(res_y, y, Vy);
xlim([-y_lim, y_lim]);

omega_y = res_y.omega;
bounds = confint(res_y);
omega_y_err = (bounds(2) - bounds(1))/2;

disp(sprintf("omega_z = %.4g +/- %.4gHz", omega_z, omega_z_err));
disp(sprintf("omega_y = %.4g +/- %.4gHz", omega_y, omega_y_err));


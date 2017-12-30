%% Calculate the trap frequency for a z shaped wire on a 2DEG
% Based on the paper "Atom chips with free-standing two-dimensional 
% electron gases: advantages and challenges" by G. A. Sinuco-Leon, P.
% Kruger, and T.M. Fromhold
clear all;

%% Parameters and constants
V0 = 1.6e-2 * 640; % Voltage difference of V0 across wire
n = 3.3e15; % Mean electron density of 2DEG [m^-2]
mu = 140; % Mobility of 2DEG [m^2 V^-2 s^-1]
mu_0 = 4e-7 * pi; % Permeability of free space
B_bias_factor = 0.9; % B_bias = B_bias_factor * Bs
B_offset_factor = 0.1; % B_offset = B_offset_factor * Bs

x_fixed = 0; % z position to evaluate field at [m]

add_noise = false;
should_plot = false;

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

[Bx1, By1, Bz1] = calc_field(x, y, Jx, Jy, dx, dy, x_fixed, yq, zq);

[Bx0, By0, Bz0] = calc_field(x, y, Jx, Jy, dx, dy, 0, 0, 0.5e-12);
Bs = sqrt(Bx0.^2 + By0.^2 + Bz0.^2);

B_bias = B_bias_factor*Bs;
B_offset = B_offset_factor*Bs;

%% Add bias/offset and plot results

By1 = By1 + B_offset;
Bx1 = Bx1 + B_bias;

B1 = sqrt(Bx1.^2 + By1.^2 + Bz1.^2);

% Plot results
% figure();
% surf(yq, zq, B1, 'EdgeColor', 'none', 'FaceColor', 'interp');
% xlabel('y [m]', 'FontSize', 18);
% ylabel('z [m]', 'FontSize', 18);
% c = colorbar();
% c.Label.String = '|B| [T]';
% c.Label.FontSize = 18;
% colormap jet;
% view(2);

dims = size(yq);
y0 = ceil(dims(2)/2);

if should_plot
    figure();
    plot(zq(:, y0), B1(:, y0));
    xlabel('z', 'FontSize', 18);
    ylabel('|B|', 'FontSize', 18);
end

B_along_z = B1(:, y0);
zs = zq(:, y0);

[o, z0_ind] = min(B_along_z);
z0 = zs(z0_ind);

disp(sprintf("z_0 = %.4g", z0));

B_along_y = B1(z0_ind, :);

xq = linspace(-150e-6, 150e-6, resolution);
[Bx2, By2, Bz2] = calc_field(x, y, Jx, Jy, dx, dy, xq, 0, z0);
Bx2 = Bx2(1, :);
By2 = By2(1, :);
Bz2 = Bz2(1, :);
By2 = By2 + B_offset;
Bx2 = Bx2 + B_bias;
B_along_x = sqrt(Bx2.^2 + By2.^2 + Bz2.^2);

%% Calculate trap frequency
% Trap freq. for an Rubidium-87 atom in the |F=2, m_F=2> state
mF = 2;
gF = 1/2; % Lande g-factor (theoretical value)
mu_B = 9.274e-24; % Bohr magneton [J/T]
m = 87 * 1.66054e-27; % Mass of atoms [kg]

Vz = mu_B .* B_along_z;
z = zq(:, y0);
[z, Vz] = prepareCurveData(z, Vz);

Vz_fit = @(omega, x) (m*omega.^2)/2 * (x-z0).^2 + mF*gF*mu_B*By1(z0_ind, y0);
ft = fittype(Vz_fit);
fo = fitoptions('Method', 'NonLinearLeastSquares');
% fo.Exclude = z > 2*z0;
fo.Exclude = or(z < 0.5*z0, z > 1.5*z0);
fo.StartPoint = 5e3;
fo.TolFun = 1e-80;
fo.MaxIter = 8000;
fo.MaxFunEvals = 10000;
fo.Lower = -Inf;
fo.Upper = Inf;
fo.DiffMaxChange = 1e3;
fo.DiffMinChange = 1e-5;
[res_z, gof_z] = fit(z, Vz, ft, fo);
if should_plot
    figure();
    subplot(1,3,1);
    hold on;
    plot(res_z, z, Vz);
    plot(0.5.*[z0 z0], [min(Vz) max(Vz)], 'g');
    plot(1.5.*[z0 z0], [min(Vz) max(Vz)], 'g');
    hold off;
    xlabel('z [m]', 'FontSize', 18);
    ylabel('V(0, z) [J]', 'FontSize', 18);
    legend({'Calculated Potential', 'Curve Fit'}, 'FontSize', 16);
end

omega_z = res_z.omega;
bounds = confint(res_z);
omega_z_err = (bounds(2) - bounds(1))/2;

Vy = mu_B .* B_along_y;
y = yq(z0_ind, :);

inds = find(imregionalmax(Vy));
if length(inds) > 2
    max_y = y(inds(3));
else
    max_y = y(inds(2));
end
y_lim = 0.5*max_y;

[y, Vy] = prepareCurveData(y, Vy);

Vy_fit = @(omega, x) (m*omega.^2)/2 * x.^2 + mF*gF*mu_B*By1(z0_ind, y0);
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
[res_y, gof_y] = fit(y, Vy, ft, fo);
if should_plot
    subplot(1, 3, 2);
    hold on;
    plot(res_y, y, Vy);
    plot(y_lim.*[-1 -1], [0 max(Vy)], 'g');
    plot(y_lim.*[1 1], [0 max(Vy)], 'g');
    hold off;
    xlabel('y [m]', 'FontSize', 18);
    ylabel('V(y, z_0) [J]', 'FontSize', 18);
    legend({'Calculated Potential', 'Curve Fit'}, 'FontSize', 16);
end

omega_y = res_y.omega;
bounds = confint(res_y);
omega_y_err = (bounds(2) - bounds(1))/2;

disp(sprintf("omega_z = %.4g +/- %.4gHz", omega_z, omega_z_err));
disp(sprintf("omega_y = %.4g +/- %.4gHz", omega_y, omega_y_err));

%% Calculate omega x

Vx = mu_B .* B_along_x;
x = xq;

[o, ind] = max(gradient(Vx));
x_lim = x(ind)*1.5;

[x, Vx] = prepareCurveData(x, Vx);

ft = fittype(Vy_fit);
fo = fitoptions('Method', 'NonLinearLeastSquares');
fo.Exclude = or(x < -x_lim, x > x_lim);
fo.StartPoint = 5e3;
fo.TolFun = 1e-80;
fo.MaxIter = 8000;
fo.MaxFunEvals = 10000;
fo.Lower = -Inf;
fo.Upper = Inf;
fo.DiffMaxChange = 1e3;
fo.DiffMinChange = 1e-5;
[res_x, gof_x] = fit(x, Vx, ft, fo);
if should_plot
    subplot(1, 3, 3);
    hold on;
    plot(res_x, x, Vx);
    plot(x_lim.*[-1 -1], [0, max(Vx)], 'g');
    plot(x_lim.*[1 1], [0, max(Vx)], 'g');
    hold off;
    xlabel('x [m]', 'FontSize', 18);
    ylabel('V(x, z_0) [J]', 'FontSize', 18);
    legend({'Calculated Potential', 'Curve Fit'}, 'FontSize', 16);
end

omega_x = res_x.omega;
bounds = confint(res_x);
omega_x_err = (bounds(2) - bounds(1))/2;

disp(sprintf("omega_x = %.4g +/- %.4gHz", omega_x, omega_x_err));

%% Calculate loss rate due to spin flips
h_bar = 1.054571e-34;
q = 1.6e-19;

z = zq(:, y0);
y = yq(z0_ind, :);

hz = (max(z) - min(z))/length(z);
hy = (max(y) - min(y))/length(y);

B_perp = sqrt(Bx1.^2 + (By1-By1(z0_ind, y0)).^2 + Bz1.^2);

B_perp_z = B_perp(:, y0);
dz_B_perp = gradient(abs(B_perp_z), hz);
dz_B_perp = dz_B_perp(z0_ind);

B_perp_y = B_perp(z0_ind, :);
dy_B_perp = gradient(abs(B_perp_y), hy);
dy_B_perp = dy_B_perp(y0);

% omega_z = sqrt(mF*gF*mu_B/(m*B_offset)) * abs(dz_B_perp);
% omega_y = sqrt(mF*gF*mu_B/(m*B_offset)) * abs(dy_B_perp);

omega = sqrt(omega_x.^2 + omega_y.^2 + omega_z.^2);

% omega_hat = (((mF*gF*mu_B)/(h_bar*By1(z0_ind, y0))) * omega_z.^2)^(1/3);
omega_hat = ((((mF*gF*mu_B)^2)/(h_bar*m))*dz_B_perp^2)^(1/3);

gamma_mf = (pi/(2*sqrt(q)))*omega*exp(-(omega_hat^3)/(omega^3));

disp(sprintf("omega = %.4g", omega));
disp(sprintf("omega_hat = %.4g", omega_hat));

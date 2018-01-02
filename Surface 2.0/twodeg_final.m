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
B_bias_factor = 0.9; % B_bias = B_bias_factor * Bs
B_offset_factor = 0.2; % B_offset = B_offset_factor * Bs

add_noise = false;
should_plot = false;

% Size of current elements to consider [m]
dx = 2e-6;
dy = 2e-6;

sigma = n*mu*1.6e-19; % Conductivity of 2DEG [S m^-1]

%% Calculate potential
res = 5e6; % Points per unit for interpolation
[x, y, phi] = calc_potential(V0, res);

%% Calculate current density
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
resolution = 201;
xq = linspace(-150e-6, 150e-6, resolution);
yq = linspace(-150e-6, 150e-6, resolution);
zq = linspace(0, 30e-6, resolution);

x0_ind = ceil(size(xq)/2);
y0_ind = ceil(size(yq)/2);

% Field at surface of conductor
[Bx0, By0, Bz0] = calc_field(x, y, Jx, Jy, dx, dy, 0, 0, 1e-11);
Bs = sqrt(Bx0.^2 + By0.^2 + Bz0.^2);

B_bias = B_bias_factor*Bs;
B_offset = B_offset_factor*Bs;

% Calculate B for r=(0, 0, z)
[Bx3, By3, Bz3] = calc_field(x, y, Jx, Jy, dx, dy, 0, 0, zq);
By3 = By3 + B_offset;
Bx3 = Bx3 + B_bias;
B3 = sqrt(Bx3.^2 + By3.^2 + Bz3.^2);

% Find z0
[o, z0_ind] = min(B3);
z0 = zq(z0_ind);

% Calculate B for r=(0, y, z0)
[Bx2, By2, Bz2] = calc_field(x, y, Jx, Jy, dx, dy, 0, yq, z0);
By2 = By2 + B_offset;
Bx2 = Bx2 + B_bias;
B2 = sqrt(Bx2.^2 + By2.^2 + Bz2.^2);

% Calculate B for r=(x, 0, z0)
[Bx1, By1, Bz1] = calc_field(x, y, Jx, Jy, dx, dy, xq, 0, z0);
Bz1 = Bz1(1,:);
By1 = By1(1,:) + B_offset;
Bx1 = Bx1(1,:) + B_bias;
B1 = sqrt(Bx1.^2 + By1.^2 + Bz1.^2);

%% Define constants for frequency calculations
% Trap freq. for an Rubidium-87 atom in the |F=2, m_F=2> state
mF = 2;
gF = 1/2; % Lande g-factor (theoretical value)
mu_B = 9.274e-24; % Bohr magneton [J/T]
m = 87 * 1.66054e-27; % Mass of atoms [kg]

%% Calculate omega_z
Vz = mu_B .* B3;
[z, Vz] = prepareCurveData(zq, Vz);

Vz_fit = @(omega, x) (m*omega.^2)/2 * (x-z0).^2 + mF*gF*mu_B*B_offset;
ft = fittype(Vz_fit);
fo = fitoptions('Method', 'NonLinearLeastSquares');
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

disp(sprintf("omega_z = %.4g +/- %.4gHz", omega_z, omega_z_err));

%% Calculate omega_y
Vy = mu_B .* B2;

inds = find(imregionalmax(Vy));
if length(inds) > 2
    max_y = yq(inds(3));
else
    max_y = yq(inds(2));
end
y_lim = 0.5*max_y;

[y, Vy] = prepareCurveData(yq, Vy);

Vy_fit = @(omega, x) (m*omega.^2)/2 * x.^2 + mF*gF*mu_B*B_offset;
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
    figure();
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

disp(sprintf("omega_y = %.4g +/- %.4gHz", omega_y, omega_y_err));

%% Calculate omega_x

Vx = mu_B .* B1;

[o, ind] = max(gradient(Vx));
x_lim = xq(ind)*1.5;

[x, Vx] = prepareCurveData(xq, Vx);

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
    figure();
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

hz = (max(zq) - min(zq))/length(zq);
hy = (max(yq) - min(yq))/length(yq);

B_perp_z = sqrt(Bx3.^2 + (By3-B_offset).^2 + Bz3.^2);
dz_B_perp = gradient(abs(B_perp_z), hz);
dz_B_perp = dz_B_perp(z0_ind);

B_perp_y = sqrt(Bx2.^2 + (By2-B_offset).^2 + Bz3.^2);
dy_B_perp = gradient(abs(B_perp_y), hy);
dy_B_perp = dy_B_perp(y0_ind);

% omega_z = sqrt(mF*gF*mu_B/(m*B_offset)) * abs(dz_B_perp);
% omega_y = sqrt(mF*gF*mu_B/(m*B_offset)) * abs(dy_B_perp);

omega = sqrt(omega_x.^2 + omega_y.^2 + omega_z.^2);

% omega_hat = (((mF*gF*mu_B)/(h_bar*By1(z0_ind, y0))) * omega_z.^2)^(1/3);
omega_hat = ((((mF*gF*mu_B)^2)/(h_bar*m))*dz_B_perp^2)^(1/3);

gamma_mf = (pi/(2*sqrt(exp(1))))*omega*exp(-(omega_hat^3)/(omega^3));

disp(sprintf("omega = %.4g", omega));
disp(sprintf("omega_hat = %.4g", omega_hat));
disp(sprintf("loss_rate = %.4g", gamma_mf));

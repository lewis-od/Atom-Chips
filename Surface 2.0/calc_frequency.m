function [omega, omega_err, r2, omega_hat, gamma_mf, gamma_err] = calc_frequency(x, y, Jx, Jy, bias_factor, offset_factor)
%calc_frequency Calculate trap distance, trap frequency, critical
%frequency, and loss rate for the given bias and offset fields

%% Apply finite element method to calculate magnetic field

% Size of current elements to consider [m]
dx = 2e-6;
dy = 2e-6;

% Points to calculate B at
resolution = 201;
xq = linspace(-150e-6, 150e-6, resolution);
yq = linspace(-150e-6, 150e-6, resolution);
zq = linspace(0, 30e-6, resolution);

% Indicies of 0 in xq and yq arrays
x0_ind = ceil(size(xq)/2);
y0_ind = ceil(size(yq)/2);

% Field at surface of conductor
[Bx0, By0, Bz0] = calc_field(x, y, Jx, Jy, dx, dy, 0, 0, 0.5e-12);
Bs = sqrt(Bx0.^2 + By0.^2 + Bz0.^2);

B_bias = bias_factor*Bs;
B_offset = offset_factor*Bs;

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

omega_z = res_z.omega;
bounds = confint(res_z);
omega_z_err = (bounds(2) - bounds(1))/2;


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

omega_y = res_y.omega;
bounds = confint(res_y);
omega_y_err = (bounds(2) - bounds(1))/2;

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

omega_x = res_x.omega;
bounds = confint(res_x);
omega_x_err = (bounds(2) - bounds(1))/2;

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

omega = [omega_x, omega_y, omega_z];
omega_err = [omega_x_err, omega_y_err, omega_z_err];
r2 = [gof_x.rsquare, gof_y.rsquare, gof_z.rsquare];

% Critical frequency [Hz]
omega_hat = ((((mF*gF*mu_B)^2)/(h_bar*m))*dz_B_perp^2)^(1/3);

% Loss rate [s^-1]
gamma_mf = (pi/(2*sqrt(exp(1))))*norm(omega)*exp(-(omega_hat^3)/(norm(omega)^3));

% Error in norm(omega)
sigma_omega = sqrt(sum(omega.^2 .* omega_err.^2))./norm(omega);
% sigma_omega = (omega_max - omega_min) / 2;

% Derivative of gamma with respect to omega
dgamma_omega = (gamma_mf/norm(omega)) * (1 + 3*(omega_hat/norm(omega))^3);

% Error in gamma
gamma_err = dgamma_omega * sigma_omega;

end


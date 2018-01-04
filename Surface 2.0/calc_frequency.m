function [omega, omega_hat, gamma_mf, z0, Bs] = calc_frequency(x, y, Jx, Jy, bias_factor, offset_factor)
%calc_frequency Calculate trap distance, trap frequency, critical
%frequency, and loss rate for the given bias and offset fields

%% Apply finite element method to calculate magnetic field

% Size of current elements to consider [m]
dx = 2e-6;
dy = 2e-6;

% Points to calculate B at
resolution = 301;
xq = linspace(-150e-6, 150e-6, resolution);
yq = linspace(-150e-6, 150e-6, resolution);
zq = linspace(0, 15e-6, resolution);

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

% Calculate B for r=(x, 0, z0)
[Bx1, By1, Bz1] = calc_field(x, y, Jx, Jy, dx, dy, xq, 0, z0);
Bz1 = Bz1(1,:);
By1 = By1(1,:) + B_offset;
Bx1 = Bx1(1,:) + B_bias;
B1 = sqrt(Bx1.^2 + By1.^2 + Bz1.^2);

%% Define constants for frequency/loss-rate calculations
% Trap freq. for an Rubidium-87 atom in the |F=2, m_F=2> state
mF = 2;
gF = 1/2; % Lande g-factor (theoretical value)
mu_B = 9.274e-24; % Bohr magneton [J/T]
m = 87 * 1.66054e-27; % Mass of atoms [kg]
h_bar = 1.054571e-34;

%% Calculate loss rate due to spin flips

hz = (max(zq) - min(zq))/length(zq);
hy = (max(yq) - min(yq))/length(yq);

dB_dz = gradient(B3, hz);
d2B_dz2 = gradient(dB_dz, hz);

omega = sqrt((mF*gF*mu_B)/m * d2B_dz2(z0_ind));

B_perp_z = sqrt(Bx3.^2 + (By3-B_offset).^2 + Bz3.^2);
dz_B_perp = gradient(abs(B_perp_z), hz);
dz_B_perp = dz_B_perp(z0_ind);

% Critical frequency [Hz]
omega_hat = ((((mF*gF*mu_B)^2)/(h_bar*m))*dz_B_perp^2)^(1/3);

% Loss rate [s^-1]
gamma_mf = (pi/(2*sqrt(exp(1))))*norm(omega)*exp(-(omega_hat^3)/(norm(omega)^3));

end

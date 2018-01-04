function [phi_noise] = calc_noise(x, n, grain_size)
%calc_noise Calculate electrostatic noise from a random ionised donor
% distribution using the Thomas-Fermi model

q = 1.6e-19; % Charge on an electron [C]
epsilon_0 = 8.85e-12; % Permittivity of free space [F m^-1]
epsilon = 12.9; % Relative permittivity of GaAs
k_s = 2.1e8; % Screening wave vector [m^-1]
d = 52.9e-9; % Distance of 2DEG from donor layer [m]

current_size = (max(max(x)) - min(min(x)))/length(x);
res = 1/current_size; % Points per unit

scale_factor = grain_size/current_size;
% scale_factor = scale_factor * res * 1e6;

% Add white noise to donor density distribution
n_dist = ones(ceil(size(x)/scale_factor)).*n;
noise = imnoise(n_dist, 'gaussian', 0, 1) - 0.5;
n_dist = noise.*n + n_dist;
n_dist = imresize(n_dist, size(x), 'nearest');

delta_n = ones(size(x)).*n - n_dist;

% Apply Thomas-Fermi screening model
delta_n_hat = fft2(delta_n);
delta_n_hat = fftshift(delta_n_hat); % Shift k=0 to centre

dims = size(x);
kx = (0:1:(dims(1)-1)).*(res/dims(1));
ky = (0:1:(dims(2)-1)).*(res/dims(2));
% Ensure k=0 is at centre
kx = kx - (max(kx)/2);
ky = ky - (max(ky)/2);
[kx, ky] = meshgrid(kx, ky);
k = sqrt(kx.^2 + ky.^2);
k = k.*(2*pi)^1; % Correct for factor of 2pi from each FFT

% Integrand of Thomas-Fermi screening equation
phi_noise_hat = (exp(-k.*d) .* delta_n_hat)./(k + k_s);

phi_noise = ifft2(phi_noise_hat);
phi_noise = abs(phi_noise);

% Scale appropriately
% Potential energy of electron [J]
phi_noise = phi_noise .* (q^2)/(4*pi*epsilon*epsilon_0);
phi_noise = phi_noise / q; % Electrostatic potential [J/C]

end

